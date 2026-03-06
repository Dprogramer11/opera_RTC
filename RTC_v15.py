
"""
OPERA RTC-S1 云原生 SAR 时序预处理流水线 v15
=============================================
所有用户配置通过命令行参数传入，由 Run_RTC.sh 统一管理。
直接查看参数说明：python RTC_v15.py --help

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
🌊 核心工作流程 (Workflow)
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
本流水线专为云原生 COG (Cloud-Optimized GeoTIFF) 数据设计，遵循“按需读取、全内存流转、单次写盘”的极致 I/O 优化原则。整体分为五个阶段：

1. 【环境与认证初始化】
   - 解析命令行参数，配置专为弱网/长连接优化的 GDAL 环境变量（增大缓存、关闭多路复用、缩短连接超时）。
   - 通过 Earthdata .netrc 验证身份，建立持久化的 Cookie 授权会话。

2. 【统一空间网格构建】
   - 根据用户输入的 ROI 边界和分辨率（默认 0.00027° ≈ 30m），计算全局统一的仿射变换矩阵（Transform）和目标宽高（Width/Height）。
   - 根据 WKT 生成布尔掩膜（AOI Mask），用于最终边界的精准裁切。

3. 【元数据分页检索】
   - 调用 asf_search 接口，按月分段检索指定时间范围内的 OPERA RTC-S1 极化切片（Burst）。
   - 按采集日期（YYYY-MM-DD）对影像进行分组，为后续的日拼接做准备。

4. 【单日并行拉取与拼接 (Warp & Mosaic)】
   - 多线程并发处理每个日期：基于统一网格，利用 GDAL /vsicurl/ 虚拟文件系统，发起 HTTP Range 请求。
   - **核心优势**：只下载与 ROI 相交的像素块（不下载完整影像），并在内存中直接重投影（Reproject）到目标坐标系。
   - 将同一天的多个 Burst 在内存中按像素累加求均值，生成单日拼接图并落盘。

5. 【时序合成与后处理 (Composite & TV & Clip)】
   - 根据配置（10d/12d/monthly）划分时间窗口，读取单日拼接图进行中值合成（Median Composite），有效去除瞬态噪声和部分条带。
   - （可选）对合成图进行 TV (Total Variation) 全变分去斑，平滑同质区域并保留边缘。
   - （可选）应用 AOI 掩膜将多边形外像素置 NaN。
   - ★ 唯一一次最终写盘：输出带有 NaN 且清洗完毕的高质量时序 TIF，直接供 Whittaker Smoothing 等平滑算法使用。

"""


import os
import sys
import time
import threading
import warnings
import logging
import argparse
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
from collections import defaultdict
from typing import Optional, Tuple

import numpy as np
import pandas as pd
import asf_search as asf
import earthaccess
import rasterio
import rasterio.windows
import rasterio.features
from rasterio.windows import from_bounds
from rasterio.warp import transform_bounds, reproject, Resampling
from rasterio.crs import CRS
from scipy.ndimage import convolve
from scipy.ndimage import label as labeler
from skimage.restoration import denoise_tv_bregman
from shapely import wkt as shapely_wkt
from shapely.geometry import mapping
from tqdm import tqdm

warnings.filterwarnings("ignore")

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger("opera_v15")


# ═══════════════════════════════════════════════════════════════
# 命令行参数解析
# ═══════════════════════════════════════════════════════════════

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="OPERA RTC-S1 时序预处理流水线 v15（由 run_RTC.sh 调用）",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # ── 研究区域 ──────────────────────────────────────────────
    p.add_argument("--aoi", required=True,
                   help="研究区 WKT 多边形字符串 POLYGON(...)")
    p.add_argument("--roi-bounds", nargs=4, type=float,
                   metavar=("MINLON", "MINLAT", "MAXLON", "MAXLAT"),
                   default=[118.0, 27.2, 123.0, 31.3],
                   help="ROI 矩形裁切范围，用于 vsicurl 读取窗口")

    # ── 时间 ─────────────────────────────────────────────────
    p.add_argument("--start",  default="2017-01-01",
                   help="起始日期 YYYY-MM-DD")
    p.add_argument("--end",    default="2017-12-31",
                   help="结束日期 YYYY-MM-DD")

    # ── 极化 ─────────────────────────────────────────────────
    p.add_argument("--polarization", default="VH",
                   choices=["VH", "VV", "HV", "HH"],
                   help="SAR 极化方式")

    # ── 输出 ─────────────────────────────────────────────────
    p.add_argument("--base-dir", default=r"E:\ESA\sentinel_1\VH",
                   help="输出根目录，自动创建 downloads/ 和 composite/ 子目录")
    p.add_argument("--resolution", type=float, default=0.00027,
                   help="输出分辨率（度/像素），0.00027≈30m")
    p.add_argument("--output-crs", default="EPSG:4326",
                   help="输出坐标系（默认 EPSG:4326，可改为 EPSG:32650 等 UTM）")

    # ── 合成 ─────────────────────────────────────────────────
    p.add_argument("--composite", default="10d",
                   help="合成模式：10d / 12d / monthly")

    # ── TV 去斑 ───────────────────────────────────────────────
    p.add_argument("--tv-reg",    type=float, default=5.0,
                   help="TV 正则化强度（越大越平滑，推荐 3~8）")
    p.add_argument("--tv-max",    type=float, default=1.0,
                   help="有效像素线性功率上限")
    p.add_argument("--tv-min",    type=float, default=1e-7,
                   help="有效像素线性功率下限")
    p.add_argument("--tv-interp", default="bilinear",
                   choices=["bilinear", "none"],
                   help="TV 前空洞插值方法")
    p.add_argument("--tv-eps",    type=float, default=0.01,
                   help="TV 收敛阈值")

    # ── 功能开关 ──────────────────────────────────────────────
    def _bool(s: str) -> bool:
        return s.strip().lower() in ("true", "1", "yes")

    p.add_argument("--do-tv",     type=_bool, default=True,
                   help="是否做 TV 去斑（true/false）")
    p.add_argument("--do-clip",   type=_bool, default=True,
                   help="是否裁切到 AOI 多边形（true/false）")
    p.add_argument("--debug-raw", type=_bool, default=False,
                   help="是否保存合成前原始拼接图（true/false）")

    # ── 中间文件保留策略 ──────────────────────────────────────
    p.add_argument("--keep-daily", type=_bool, default=True,
                   help="合成完成后是否保留单日期 warp TIF（false=自动删除，节省磁盘）")
    p.add_argument("--keep-raw-composite", type=_bool, default=False,
                   help="是否保留合成后 TV 前的原始合成图（须同时设 --debug-raw true）")

    # ── 网络 ─────────────────────────────────────────────────
    p.add_argument("--workers",     type=int, default=2,
                   help="日期并发线程数")
    p.add_argument("--max-retries", type=int, default=5,
                   help="网络错误最大重试次数")
    p.add_argument("--retry-delay", type=int, default=3,
                   help="首次重试等待秒数（指数退避底数）")

    # ── Cookie ────────────────────────────────────────────────
    p.add_argument("--cookie",
                   default=str(Path(__file__).parent / "cookies.txt"),
                   help="Earthdata Cookie 文件路径")

    return p.parse_args()


# ═══════════════════════════════════════════════════════════════
# 全局运行时状态（main 初始化后注入）
# ═══════════════════════════════════════════════════════════════

CFG: argparse.Namespace = None          # type: ignore  # 解析后的参数
OUTPUT_GRID: dict       = {}
AOI_MASK: Optional[np.ndarray] = None
GDAL_CONFIG: dict       = {}

# OUTPUT_CRS 由 CFG.output_crs 动态提供（默认 EPSG:4326）
SEARCH_RETRY = 5
CB_THRESHOLD = 8     # 熔断触发阈值：连续失败次数
CB_COOLDOWN  = 90    # 熔断冷却时间（秒）
COOKIE_REFRESH_INTERVAL = 4 * 3600  # Cookie 刷新间隔（秒）

# 熔断器（线程安全）
_cb_lock                 = threading.Lock()
_cb_consecutive_failures = 0
_cb_tripped              = False
_cb_trip_time            = 0.0

_last_cookie_refresh: float = 0.0


# ═══════════════════════════════════════════════════════════════
# ★ v15 错误分类表（先 DATA 后 NETWORK，顺序不能颠倒）
# ═══════════════════════════════════════════════════════════════

# 数据/文件错误 → 立即跳过，重试无意义
_DATA_ERRORS = (
    "tifffilltile",                              # corrupt tile  err_no=1
    "tiffreadencoded",                           # corrupt tile  err_no=1
    "ireadblock failed",                         # corrupt tile  err_no=1
    "not recognized as a supported file format", # 文件损坏（纯净，无 vsicurl 断连前缀）
)

# 网络/瞬态错误 → 重试（含 vsicurl 断连后的假 404）
_NET_ERRORS = (
    "could not resolve",
    "curl error",
    "connection refused",
    "connection reset",
    "timeout",
    "ssl",
    "eof",
    "http error",
    "503", "502", "429",
    "response_code=0",
    # ★ v15 新增：CURL err_no=11 断连后 vsicurl 立即报的瞬态 404
    #   日志模式：err_no=11 -> err_no=4("does not exist...") 成对出现
    "does not exist in the file system",
    "not recognized as a supported dataset name",
    "failed to connect",
    "could not connect",
)


# ═══════════════════════════════════════════════════════════════
# Part 0: 统一输出网格
# ═══════════════════════════════════════════════════════════════

def compute_output_grid() -> dict:
    res                  = CFG.resolution
    crs_str              = CFG.output_crs          # ★ 从参数读取
    minx, miny, maxx, maxy = CFG.roi_bounds
    minx = np.floor(minx / res) * res
    miny = np.floor(miny / res) * res
    maxx = np.ceil(maxx  / res) * res
    maxy = np.ceil(maxy  / res) * res
    W    = int(round((maxx - minx) / res))
    H    = int(round((maxy - miny) / res))
    tr   = rasterio.transform.from_bounds(minx, miny, maxx, maxy, W, H)
    log.info("📐 统一网格: %dx%d px | %.5f°/px | CRS=%s | %.0f MB/景",
             W, H, res, crs_str, W * H * 4 / 1e6)
    return dict(crs=CRS.from_string(crs_str), transform=tr,
                width=W, height=H,
                bounds=(minx, miny, maxx, maxy), resolution=res)


# ═══════════════════════════════════════════════════════════════
# Part 0b: AOI 掩膜
# ═══════════════════════════════════════════════════════════════

def make_aoi_mask(grid: dict) -> np.ndarray:
    geom = shapely_wkt.loads(CFG.aoi)
    mask = rasterio.features.geometry_mask(
        geometries=[mapping(geom)],
        out_shape=(grid["height"], grid["width"]),
        transform=grid["transform"],
        invert=True,
    )
    pct = mask.sum() / mask.size * 100
    log.info("🔲 AOI 掩膜: %.1f%% 在研究区内 | %.0f MB", pct, mask.nbytes / 1e6)
    return mask


def apply_aoi_clip(arr: np.ndarray, mask: np.ndarray) -> None:
    """原地将 AOI 外像素置 NaN。"""
    arr[~mask] = np.nan


# ═══════════════════════════════════════════════════════════════
# Part 1: 认证 & GDAL 环境
# ═══════════════════════════════════════════════════════════════

def setup_gdal_env() -> None:
    """
    禁用 GDAL 内部重试（MAX_RETRY=0），连接超时缩至 12s，
    通过 gdal.SetConfigOption() 直接写入 C 层（Windows 必须）。
    """
    global GDAL_CONFIG
    cookie = str(Path(CFG.cookie).resolve()).replace("\\", "/")

    GDAL_CONFIG = {
        "GDAL_HTTP_COOKIEFILE":             cookie,
        "GDAL_HTTP_COOKIEJAR":              cookie,
        "GDAL_HTTP_DNS_CACHE_TIMEOUT":      "600",
        "GDAL_HTTP_MAX_CONNECTIONS":        "2",
        "GDAL_HTTP_MULTIPLEX":              "YES",
        "CPL_VSIL_CURL_USE_HEAD":           "NO",
        "GDAL_HTTP_CONNECTTIMEOUT":         "12",  # ★ 12s（原 30s）
        "GDAL_HTTP_TIMEOUT":                "90",
        "GDAL_HTTP_MAX_RETRY":              "0",   # ★ 禁用 GDAL 内部重试
        "GDAL_HTTP_RETRY_DELAY":            "0",
        "GDAL_DISABLE_READDIR_ON_OPEN":     "EMPTY_DIR",
        "CPL_VSIL_CURL_ALLOWED_EXTENSIONS": ".tif,.tiff",
        "CPL_VSIL_CURL_CACHE_SIZE":         "200000000",
        "GDAL_CACHEMAX":                    "512",
    }

    os.environ.update(GDAL_CONFIG)

    try:
        from osgeo import gdal as _gdal
        for k, v in GDAL_CONFIG.items():
            _gdal.SetConfigOption(k, v)
        log.info("🔧 GDAL 配置完成（内部重试=0，超时=12s，osgeo 直写）")
    except ImportError:
        log.warning("⚠️  osgeo.gdal 不可用，仅通过 os.environ 设置")


def _cookie_is_fresh(max_age_hours: int = 20) -> bool:
    """Cookie 文件存在、大小 > 500B 且修改时间在 max_age_hours 内。"""
    p = Path(CFG.cookie)
    if not p.exists() or p.stat().st_size < 500:
        return False
    return (time.time() - p.stat().st_mtime) / 3600 < max_age_hours


def authenticate() -> None:
    """
    认证策略（v15 修复）：
      1. Cookie 文件存在且新鲜（20h 内）→ 直接复用，不重新登录
      2. 否则（Cookie 不存在 / 已过期）→ 通过 ~/.netrc 重新登录
         登录后若 earthaccess 写出了 cookie 则配置到 GDAL；
         即使未写出 cookie 文件，earthaccess 的内存 session 仍可用。
    """
    global _last_cookie_refresh
    p = Path(CFG.cookie)

    # ── 路径 1：Cookie 有效，直接复用 ─────────────────────────────
    if _cookie_is_fresh(20):
        age_h = (time.time() - p.stat().st_mtime) / 3600
        log.info("🔐 Cookie 有效（%.0fh 内），跳过登录", age_h)
        setup_gdal_env()
        _last_cookie_refresh = time.time()
        return

    # ── 路径 2：通过 .netrc 登录 ──────────────────────────────────
    if not p.exists():
        log.info("🔐 Cookie 文件不存在，通过 ~/.netrc 登录...")
    else:
        log.info("🔐 Cookie 已过期，通过 ~/.netrc 重新登录...")

    last_exc: Exception = RuntimeError("未知错误")
    for attempt in range(3):
        try:
            auth = earthaccess.login(strategy="netrc")
            if not auth:
                raise RuntimeError(
                    "earthaccess.login 返回 False\n"
                    "  请确认 ~/.netrc 包含 urs.earthdata.nasa.gov 的凭据，格式：\n"
                    "  machine urs.earthdata.nasa.gov login <USER> password <PASS>"
                )
            # 登录成功：检查 cookie 是否被写出
            if p.exists() and p.stat().st_size > 100:
                log.info("✅ 登录成功 | Cookie %.1f KB（已写入磁盘）",
                         p.stat().st_size / 1024)
            else:
                log.info("✅ 登录成功（earthaccess 内存 session，未生成 cookie 文件）")
            setup_gdal_env()
            _last_cookie_refresh = time.time()
            return
        except Exception as e:
            last_exc = e
            es = str(e).lower()
            if any(k in es for k in ["ssl", "eof", "connection",
                                      "timeout", "max retries"]):
                wait = 2 ** (attempt + 1)
                log.warning("⚠️  网络抖动，%ds 后重试（%d/3）", wait, attempt + 1)
                time.sleep(wait)
                continue
            # 非网络错误（如凭据错误）直接抛出
            raise
    raise RuntimeError(f"Earthdata .netrc 登录失败（3次均失败）：{last_exc}")


def maybe_refresh_cookie() -> None:
    """每 4h 检查一次认证状态，必要时重新登录。"""
    global _last_cookie_refresh
    now = time.time()
    if now - _last_cookie_refresh < COOKIE_REFRESH_INTERVAL:
        return
    p = Path(CFG.cookie)
    if p.exists() and _cookie_is_fresh(3):
        # Cookie 仍有效，只需刷新 GDAL 配置
        setup_gdal_env()
        _last_cookie_refresh = now
        log.info("🔄 GDAL 环境已刷新（Cookie 仍有效）")
    else:
        log.warning("⚠️  Cookie 不存在或即将过期，后台重新登录...")
        try:
            authenticate()
        except Exception as e:
            log.warning("⚠️  后台登录失败，继续使用当前 session：%s", e)


# ═══════════════════════════════════════════════════════════════
# Part 1b: 熔断器
# ═══════════════════════════════════════════════════════════════

def _cb_record_failure() -> None:
    global _cb_consecutive_failures, _cb_tripped, _cb_trip_time
    with _cb_lock:
        _cb_consecutive_failures += 1
        if _cb_consecutive_failures >= CB_THRESHOLD and not _cb_tripped:
            _cb_tripped   = True
            _cb_trip_time = time.time()
            log.warning("⚡ 熔断器触发：连续 %d 次网络失败，暂停 %ds...",
                        _cb_consecutive_failures, CB_COOLDOWN)


def _cb_record_success() -> None:
    global _cb_consecutive_failures, _cb_tripped
    with _cb_lock:
        _cb_consecutive_failures = 0
        if _cb_tripped:
            _cb_tripped = False
            log.info("✅ 熔断器恢复")


def _cb_wait_if_tripped() -> None:
    global _cb_tripped, _cb_consecutive_failures
    while True:
        with _cb_lock:
            if not _cb_tripped:
                return
            elapsed = time.time() - _cb_trip_time
            if elapsed >= CB_COOLDOWN:
                _cb_tripped              = False
                _cb_consecutive_failures = 0
                log.info("⚡ 熔断器冷却完成，恢复请求")
                return
            remaining = CB_COOLDOWN - elapsed
        time.sleep(min(5.0, remaining))


# ═══════════════════════════════════════════════════════════════
# Part 2: 搜索 + 覆盖率诊断
# ═══════════════════════════════════════════════════════════════

def _get_url_tif(row: pd.Series, pol: str) -> str:
    urls  = ([row.url] if pd.notna(row.url) else []) + (row.additionalUrls or [])
    valid = [u for u in urls if (
        ("_VH.tif" in u or "_HV.tif" in u) if pol == "crosspol"
        else ("_VV.tif" in u or "_HH.tif" in u)
    )]
    return valid[0] if valid else ""


def diagnose_coverage(date_groups: dict) -> None:
    log.info("=" * 60)
    log.info("📊 轨道覆盖诊断")
    log.info("=" * 60)
    track_dates: dict = defaultdict(list)
    all_dates = sorted(date_groups.keys())
    for date_str, bursts in date_groups.items():
        for b in bursts:
            if t := b.get("track"):
                track_dates[t].append(date_str)
    for track in sorted(track_dates.keys()):
        dates = sorted(set(track_dates[track]))
        log.info("  轨道 T%03d: %d 次 -> %s...",
                 track, len(dates), ", ".join(dates[:4]))
    if len(all_dates) >= 2:
        intervals = [(pd.to_datetime(all_dates[i]) -
                      pd.to_datetime(all_dates[i - 1])).days
                     for i in range(1, len(all_dates))]
        log.info("  重访中位间隔 %d 天 | 最大 %d 天",
                 int(np.median(intervals)), max(intervals))
    y = pd.to_datetime(CFG.start).year
    hint = ("月度合成" if y <= 2016 or y > 2021 else "12天或月度合成")
    log.info("  推荐合成周期（%d 年）：%s", y, hint)
    log.info("=" * 60)


def search_and_group_by_date() -> dict:
    log.info("🔍 搜索 OPERA RTC-S1 [%s ~ %s]...", CFG.start, CFG.end)

    seg_starts = pd.date_range(CFG.start, CFG.end, freq="MS")
    segs = []
    for i, s0 in enumerate(seg_starts):
        s1 = (seg_starts[i + 1] - pd.Timedelta(days=1)
              if i + 1 < len(seg_starts)
              else pd.to_datetime(CFG.end))
        segs.append((s0.strftime("%Y-%m-%d"), s1.strftime("%Y-%m-%d")))

    log.info("  分段数：%d（每月一段）", len(segs))
    all_res = asf.ASFSearchResults([])

    for idx, (s, e) in enumerate(segs, 1):
        n_seg = 0
        for attempt in range(SEARCH_RETRY):
            try:
                for page in asf.search_generator(
                    dataset=asf.DATASET.OPERA_S1,
                    processingLevel="RTC",
                    start=s, end=e,
                    intersectsWith=CFG.aoi,
                ):
                    all_res.extend(page)
                    n_seg += len(page)
                log.info("  段%d/%d [%s~%s]：%d 景",
                         idx, len(segs), s, e, n_seg)
                break
            except Exception as ex:
                if attempt < SEARCH_RETRY - 1:
                    wait = 5 * (attempt + 1)
                    log.warning("  ⚠️  段%d 搜索失败(%d/%d)，%ds 后重试：%s",
                                idx, attempt + 1, SEARCH_RETRY, wait, ex)
                    time.sleep(wait)
                else:
                    log.warning("  ⚠️  段%d [%s~%s] 最终失败，跳过：%s",
                                idx, s, e, ex)

    log.info("✅ 共找到 %d 景", len(all_res))
    if not all_res:
        return {}

    pol_key = "crosspol" if CFG.polarization in ("VH", "HV") else "copol"
    groups:  dict = defaultdict(list)

    for r in all_res:
        prop = r.properties
        if CFG.polarization not in prop.get("polarization", []):
            continue
        row = pd.Series({"url":            prop.get("url", ""),
                         "additionalUrls": prop.get("additionalUrls", [])})
        url = _get_url_tif(row, pol_key)
        if not url:
            continue
        acq_dt   = pd.to_datetime(prop.get("startTime"))
        date_key = acq_dt.strftime("%Y-%m-%d")
        groups[date_key].append({
            "opera_id":     prop.get("sceneName", ""),
            "url":          url,
            "acq_datetime": acq_dt,
            "track":        prop.get("pathNumber"),
        })

    n_d = len(groups)
    n_t = sum(len(v) for v in groups.values())
    log.info("📅 %d 日期 | %d 景 | 均 %.1f 景/日",
             n_d, n_t, n_t / max(n_d, 1))
    return dict(sorted(groups.items()))


# ═══════════════════════════════════════════════════════════════
# Part 3: TV 全变分去斑
# ═══════════════════════════════════════════════════════════════

def _ext_mask(img: np.ndarray) -> np.ndarray:
    is_nan = np.isnan(img)
    if not is_nan.any():
        return np.zeros_like(img, dtype=np.uint8)
    labels, _ = labeler(is_nan)
    edge = np.zeros_like(is_nan, dtype=bool)
    edge[0, :] = edge[-1, :] = edge[:, 0] = edge[:, -1] = True
    ext_ids = np.unique(labels[edge & is_nan]).tolist()
    return np.isin(labels, ext_ids).astype(np.uint8)


def _fill_bilinear(arr: np.ndarray, n: int = 10) -> np.ndarray:
    f, m = arr.copy(), np.isnan(arr)
    if not m.any():
        return f
    k = np.array([[1, 1, 1], [1, 0, 1], [1, 1, 1]], dtype=float)
    for _ in range(n):
        if not m.any():
            break
        t = np.where(m, 0., f)
        w = convolve((~m).astype(float), k, mode="constant", cval=0)
        s = convolve(t,                  k, mode="constant", cval=0)
        with np.errstate(invalid="ignore"):
            interp = np.where(w > 0, s / w, np.nan)
        upd = m & (w > 0)
        f[upd] = interp[upd]
        m = np.isnan(f)
    f[_ext_mask(arr).astype(bool)] = np.nan
    return f


def despeckle_tv(composite: np.ndarray) -> np.ndarray:
    log.info("  🎨 TV 去斑 (reg=%.0f, interp=%s, eps=%.3f)...",
             CFG.tv_reg, CFG.tv_interp, CFG.tv_eps)
    t0  = time.time()
    X_c = np.clip(composite, CFG.tv_min, CFG.tv_max)

    if CFG.tv_interp == "bilinear":
        X_c = _fill_bilinear(X_c)

    nodata = np.isnan(X_c)
    X_db   = 10 * np.log10(X_c, out=np.full(X_c.shape, np.nan), where=~nodata)
    X_db[nodata] = -23

    reg   = max(CFG.tv_reg - 2, 1.0)
    X_dsp = denoise_tv_bregman(X_db, weight=1.0 / reg,
                               isotropic=True, eps=CFG.tv_eps)
    X_out          = np.power(10, X_dsp / 10.0).astype(np.float32)
    X_out[nodata]  = np.nan
    result         = np.clip(X_out, CFG.tv_min, CFG.tv_max)
    log.info("  ✅ TV 完成 (%.1fs)", time.time() - t0)
    return result


# ═══════════════════════════════════════════════════════════════
# Part 4: 单 burst 读取 + warp
# ═══════════════════════════════════════════════════════════════

def _classify_error(err_str: str) -> str:
    """
    返回 'data'（立即跳过）、'network'（重试）或 'unknown'（跳过）。

    ★ _DATA_ERRORS 必须先判：
    "Read or write failed...TIFFReadEncodedTile" 在 v14 因顺序错误
    被归为网络错误重试了 3 次，v15 通过先检查 data 关键词解决。
    """
    s = err_str.lower()
    for pat in _DATA_ERRORS:
        if pat in s:
            return "data"
    for pat in _NET_ERRORS:
        if pat in s:
            return "network"
    return "unknown"


def read_and_warp_burst(url: str, grid: dict) -> Optional[np.ndarray]:
    vurl = url if url.startswith("/vsicurl/") else f"/vsicurl/{url}"
    roi  = tuple(CFG.roi_bounds)
    last_err: Optional[Exception] = None

    maybe_refresh_cookie()

    for attempt in range(CFG.max_retries):
        _cb_wait_if_tripped()

        try:
            with rasterio.Env():
                with rasterio.open(vurl) as src:
                    b   = transform_bounds(CRS.from_epsg(4326), src.crs, *roi)
                    win = from_bounds(*b, transform=src.transform).intersection(
                        rasterio.windows.Window(
                            0, 0, int(src.width), int(src.height))
                    ).round_offsets().round_lengths()

                    if win.width <= 0 or win.height <= 0:
                        return None   # ROI 无交集，不重试

                    arr    = src.read(1, window=win).astype(np.float32)
                    win_tr = src.window_transform(win)
                    crs    = src.crs
                    nd     = src.nodata if src.nodata is not None else 0.0

            arr[arr == nd] = np.nan
            arr[arr <= 0]  = np.nan
            if np.sum(np.isfinite(arr)) < 100:
                del arr
                return None

            dst = np.full((grid["height"], grid["width"]),
                          np.nan, dtype=np.float32)
            reproject(
                source=arr, destination=dst,
                src_transform=win_tr, src_crs=crs, src_nodata=np.nan,
                dst_transform=grid["transform"], dst_crs=grid["crs"],
                dst_nodata=np.nan, resampling=Resampling.bilinear,
            )
            del arr
            _cb_record_success()
            return dst

        except Exception as e:
            last_err  = e
            err_class = _classify_error(str(e))

            if err_class == "data":
                log.warning("  ⚠️  数据错误(跳过)：%s", str(e)[:120])
                return None

            elif err_class == "network":
                _cb_record_failure()
                if attempt < CFG.max_retries - 1:
                    wait = CFG.retry_delay * (2 ** attempt)  # 3/6/12/24/48s
                    log.warning("  ⚠️  网络错误(%d/%d)，%ds 后重试：%s",
                                attempt + 1, CFG.max_retries,
                                wait, str(e)[:80])
                    time.sleep(wait)
                    continue
                log.warning("  ❌ 网络重试 %d 次耗尽，跳过：%s",
                            CFG.max_retries, str(e)[:80])
                return None

            else:
                log.warning("  ⚠️  未知错误(跳过)：%s", str(e)[:120])
                return None

    log.warning("  ❌ 重试耗尽，跳过：%s",
                str(last_err)[:80] if last_err else "")
    return None


# ═══════════════════════════════════════════════════════════════
# Part 5: 单日期 warp + mosaic
# ═══════════════════════════════════════════════════════════════

def process_one_date(date_str: str,
                     bursts: list,
                     grid: dict,
                     output_dir: Path) -> Optional[Tuple[Path, float]]:
    out = output_dir / f"{date_str}_{CFG.polarization}.tif"
    if out.exists():
        with rasterio.open(out) as s:
            arr = s.read(1)
        return out, float(np.isfinite(arr).sum() / arr.size * 100)

    H, W = grid["height"], grid["width"]
    acc  = np.zeros((H, W), dtype=np.float64)
    cnt  = np.zeros((H, W), dtype=np.uint16)
    n_ok = 0
    t0   = time.time()

    for b in bursts:
        dst = read_and_warp_burst(b["url"], grid)
        if dst is None:
            continue
        v = np.isfinite(dst)
        acc[v] += dst[v]
        cnt[v] += 1
        del dst, v
        n_ok += 1
        if n_ok % 10 == 0:
            log.info("  [%s] burst %d/%d (%.1fs/景)",
                     date_str, n_ok, len(bursts),
                     (time.time() - t0) / max(n_ok, 1))

    if n_ok == 0:
        del acc, cnt
        return None

    with np.errstate(divide="ignore", invalid="ignore"):
        mosaic = np.where(cnt > 0, (acc / cnt).astype(np.float32), np.nan)
    del acc, cnt
    vp = float(np.isfinite(mosaic).sum() / mosaic.size * 100)

    out.parent.mkdir(parents=True, exist_ok=True)
    _write_tif(mosaic, out, grid)
    del mosaic

    log.info("  ✅ [%s] %d burst | 覆盖 %.1f%% | %.1fs",
             date_str, n_ok, vp, time.time() - t0)
    return out, vp


# ═══════════════════════════════════════════════════════════════
# Part 6: 合成 + TV + AOI 裁切（全内存，单次写盘）
# ═══════════════════════════════════════════════════════════════

def _composite_block(files: list, grid: dict,
                     chunk_cols: int = 512) -> np.ndarray:
    H, W   = grid["height"], grid["width"]
    result = np.full((H, W), np.nan, dtype=np.float32)
    for c0 in range(0, W, chunk_cols):
        c1   = min(c0 + chunk_cols, W)
        win  = rasterio.windows.Window(c0, 0, c1 - c0, H)
        arrs = [rasterio.open(fp).read(1, window=win).astype(np.float32)
                for fp in files]
        stack = np.stack(arrs, axis=0)
        result[:, c0:c1] = np.nanmedian(stack, axis=0).astype(np.float32)
        del stack, arrs
    return result


def _write_tif(arr: np.ndarray, path: Path, grid: dict) -> None:
    H, W = grid["height"], grid["width"]
    with rasterio.open(
        path, "w", driver="GTiff", dtype="float32",
        count=1, crs=grid["crs"], transform=grid["transform"],
        width=W, height=H, nodata=np.nan,
        compress="lzw", tiled=True, blockxsize=512, blockysize=512,
    ) as f:
        f.write(arr[np.newaxis])


def _build_windows(mf: pd.DataFrame):
    """根据 --composite 参数构建时间窗口列表。"""
    mode = CFG.composite
    if mode == "monthly":
        mf["ym"] = mf.dt.dt.to_period("M")
        return (
            [(str(ym)[:7], mf[mf.ym == ym].path.tolist())
             for ym, _ in mf.groupby("ym")],
            "monthly",
        )
    step   = int(mode.replace("d", "")) if mode.endswith("d") else 10
    start  = mf.dt.min()
    end    = mf.dt.max() + pd.Timedelta(days=step)
    wdates = pd.date_range(start, end, freq=f"{step}D")
    windows = []
    for w0 in wdates[:-1]:
        w1    = w0 + pd.Timedelta(days=step)
        files = mf[(mf.dt >= w0) & (mf.dt < w1)].path.tolist()
        windows.append((w0.strftime("%Y-%m-%d"), files))
    return windows, mode


def composite_and_despeckle(manifest: pd.DataFrame,
                             grid: dict,
                             output_dir: Path,
                             aoi_mask: Optional[np.ndarray]) -> pd.DataFrame:
    log.info("📅 合成+TV+裁切（mode=%s）...", CFG.composite)
    output_dir.mkdir(parents=True, exist_ok=True)

    mf = manifest.copy()
    mf["dt"] = pd.to_datetime(mf.date)
    windows, suffix = _build_windows(mf)

    records = []
    H, W    = grid["height"], grid["width"]
    pol     = CFG.polarization

    for date_str, files in tqdm(windows, desc="合成+TV+裁切", unit="窗口"):
        tag = f"{date_str}_{pol}_{suffix}"

        # 决定输出文件名
        parts = []
        if CFG.do_tv:   parts.append("TV")
        if CFG.do_clip and aoi_mask is not None: parts.append("clip")
        suffix_out = ("_" + "_".join(parts)) if parts else ""
        final_out  = output_dir / f"{tag}{suffix_out}.tif"
        raw_out    = output_dir / f"{tag}_raw.tif"

        # 已存在则跳过
        if final_out.exists():
            with rasterio.open(final_out) as s:
                arr = s.read(1)
            vp = float(np.isfinite(arr).sum() / arr.size * 100)
            log.info("  [%s] 已存在，跳过（覆盖 %.1f%%）", date_str, vp)
            records.append({"date": date_str, "path": str(final_out),
                            "n_inputs": len(files), "valid_pct": round(vp, 1)})
            continue

        # 无输入景
        if not files:
            empty = np.full((H, W), np.nan, dtype=np.float32)
            _write_tif(empty, final_out, grid)
            del empty
            records.append({"date": date_str, "path": str(final_out),
                            "n_inputs": 0, "valid_pct": 0.0})
            continue

        # 合成
        log.info("  [%s] 合成 %d 景...", date_str, len(files))
        composite = _composite_block(files, grid)
        vp_now    = float(np.isfinite(composite).sum() / composite.size * 100)
        log.info("  [%s] 合成覆盖=%.1f%%", date_str, vp_now)

        if CFG.debug_raw:
            _write_tif(composite, raw_out, grid)

        # TV 去斑
        if CFG.do_tv:
            composite = despeckle_tv(composite)
            vp_now    = float(np.isfinite(composite).sum() / composite.size * 100)

        # AOI 裁切
        if CFG.do_clip and aoi_mask is not None:
            apply_aoi_clip(composite, aoi_mask)
            vp_final = float(np.isfinite(composite).sum() / composite.size * 100)
            log.info("  [%s] 裁切后覆盖=%.1f%%", date_str, vp_final)
        else:
            vp_final = vp_now

        _write_tif(composite, final_out, grid)
        del composite

        # 清理 raw 合成图（debug_raw=true 且 keep_raw_composite=false 时删除）
        if CFG.debug_raw and not CFG.keep_raw_composite and raw_out.exists():
            raw_out.unlink()
            log.info("  🧹 已删除原始合成图：%s", raw_out.name)

        log.info("  ✅ [%s] -> %s (覆盖=%.1f%%)",
                 date_str, final_out.name, vp_final)
        records.append({"date": date_str, "path": str(final_out),
                        "n_inputs": len(files), "valid_pct": round(vp_final, 1)})

    df = pd.DataFrame(records)
    log.info("✅ 合成完成：%d 时间步 | 平均覆盖=%.1f%%",
             len(df), df.valid_pct.mean() if len(df) > 0 else 0.0)
    return df


# ═══════════════════════════════════════════════════════════════
# Part 7: 主流水线
# ═══════════════════════════════════════════════════════════════

def run_pipeline() -> pd.DataFrame:
    global OUTPUT_GRID, AOI_MASK

    base_dir      = Path(CFG.base_dir)
    output_dir    = base_dir / "downloads"
    composite_dir = output_dir / "composite"

    t_total = time.time()
    log.info("=" * 72)
    log.info("🚀 OPERA RTC-S1 时序预处理流水线 v15")
    log.info("   时间：%s ~ %s  |  极化：%s", CFG.start, CFG.end, CFG.polarization)
    log.info("   合成：%-8s  |  TV：%s (reg=%.0f)  |  AOI裁切：%s",
             CFG.composite,
             "✅" if CFG.do_tv   else "❌", CFG.tv_reg,
             "✅" if CFG.do_clip else "❌")
    log.info("   分辨率：%.5f°  |  并发：%d  |  重试：%d×（初始%ds）",
             CFG.resolution, CFG.workers, CFG.max_retries, CFG.retry_delay)
    log.info("   输出：%s", base_dir)
    log.info("=" * 72)

    # 认证
    authenticate()

    # 网格 & 掩膜
    OUTPUT_GRID = compute_output_grid()
    AOI_MASK    = make_aoi_mask(OUTPUT_GRID) if CFG.do_clip else None

    # 搜索
    date_groups = search_and_group_by_date()
    if not date_groups:
        log.error("未找到任何景，退出")
        return pd.DataFrame()

    diagnose_coverage(date_groups)
    n_dates  = len(date_groups)
    n_scenes = sum(len(v) for v in date_groups.values())
    log.info("\n📋 %d 日期 | %d 景 | %d 并发线程", n_dates, n_scenes, CFG.workers)

    # 连通性验证（尝试前 5 个 burst）
    log.info("🔎 验证数据连通性...")
    valid_ok = False
    for fb in date_groups[next(iter(date_groups))][:5]:
        log.info("   尝试：%s", fb["opera_id"])
        try:
            dst = read_and_warp_burst(fb["url"], OUTPUT_GRID)
            if dst is not None:
                log.info("✅ 验证通过 | shape=%s | 非NaN=%.1f%%",
                         dst.shape, np.isfinite(dst).sum() / dst.size * 100)
                del dst
                valid_ok = True
                break
            log.info("   无数据，尝试下一个...")
        except Exception as e:
            log.warning("   异常：%s", e)

    if not valid_ok:
        log.error("❌ 验证失败：请检查 Cookie 是否有效（删除 cookies.txt 后重新运行）")
        return pd.DataFrame()

    # 单日期并发处理
    output_dir.mkdir(parents=True, exist_ok=True)
    records: list = []
    success = failed = 0
    t_s = time.time()

    with ThreadPoolExecutor(max_workers=CFG.workers) as ex:
        futs = {
            ex.submit(process_one_date, ds, bursts, OUTPUT_GRID, output_dir): ds
            for ds, bursts in date_groups.items()
        }
        with tqdm(total=len(futs), desc="单日期warp", unit="日") as pbar:
            for fut in as_completed(futs):
                ds = futs[fut]
                try:
                    res = fut.result()
                    if res:
                        p, vp = res
                        records.append({
                            "date":      ds,
                            "path":      str(p),
                            "n_bursts":  len(date_groups[ds]),
                            "valid_pct": round(vp, 1),
                        })
                        success += 1
                    else:
                        failed += 1
                except Exception as e:
                    log.warning("[%s] 异常：%s", ds, e)
                    failed += 1
                rate = success / (time.time() - t_s) * 60 if success else 0
                pbar.set_postfix(ok=success, fail=failed, rate=f"{rate:.1f}/分")
                pbar.update(1)

    manifest = (pd.DataFrame(records)
                .sort_values("date").reset_index(drop=True))
    manifest.to_csv(output_dir / "manifest.csv", index=False)
    log.info("✅ 单日期完成：%.1f分 | 成功=%d  失败=%d",
             (time.time() - t_s) / 60, success, failed)

    if manifest.empty:
        return manifest

    # 合成 + TV + 裁切
    t_c     = time.time()
    comp_df = composite_and_despeckle(manifest, OUTPUT_GRID,
                                      composite_dir, AOI_MASK)
    comp_df.to_csv(composite_dir / "manifest_comp.csv", index=False)
    log.info("📦 合成完成：%.1f分", (time.time() - t_c) / 60)

    # ── 清理单日期中间文件 ───────────────────────────────────────
    if not CFG.keep_daily and not manifest.empty:
        log.info("🧹 清理单日期 warp 文件（keep_daily=false）...")
        n_del = n_miss = 0
        for rec in records:
            fp = Path(rec["path"])
            if fp.exists():
                fp.unlink()
                n_del += 1
            else:
                n_miss += 1
        # 如果 manifest.csv 也不需要保留，可同时删除；此处保留 csv 方便排查
        log.info("   🧹 已删除 %d 个单日期文件（%d 个已不存在）", n_del, n_miss)
    else:
        log.info("   📁 单日期文件已保留（keep_daily=true）")

    # 汇总
    total = (time.time() - t_total) / 60
    log.info("=" * 72)
    log.info("🎉 全部完成！总耗时 %.1f 分钟", total)
    log.info("   单日期均值覆盖：%.1f%%  |  合成均值覆盖：%.1f%%  |  时间步：%d",
             manifest.valid_pct.mean(),
             comp_df.valid_pct.mean() if not comp_df.empty else 0.0,
             len(comp_df))
    log.info("   坐标系：%s  |  单日期文件：%s  |  原始合成：%s",
             CFG.output_crs,
             "已保留" if CFG.keep_daily else "已删除",
             "已保留" if (CFG.debug_raw and CFG.keep_raw_composite) else "未生成/已删除")
    log.info("   输出目录：%s", base_dir)
    log.info("=" * 72)
    return manifest


# ═══════════════════════════════════════════════════════════════
# 程序入口
# ═══════════════════════════════════════════════════════════════

if __name__ == "__main__":
    CFG = parse_args()

    manifest = run_pipeline()

    if not manifest.empty:
        print("\n📋 单日期清单：")
        print(manifest.to_string(index=False))

        cp = Path(CFG.base_dir) / "downloads" / "composite" / "manifest_comp.csv"
        if cp.exists():
            df = pd.read_csv(cp)
            print(f"\n📊 合成清单（{CFG.composite}，含NaN）：")
            print(df.to_string(index=False))
            print(f"\n  时间步：{len(df)}  |  平均覆盖：{df.valid_pct.mean():.1f}%")
