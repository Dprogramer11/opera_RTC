---

# 🌊 OPERA RTC-S1 云原生 SAR 时序预处理流水线 (v15)

本流水线专为云原生 COG (Cloud-Optimized GeoTIFF) 数据设计，遵循“按需读取、全内存流转、单次写盘”的极致 I/O 优化原则。通过 GDAL 虚拟文件系统 (`/vsicurl/`)，直接在云端对 OPERA RTC-S1 极化切片（Burst）进行拉取、拼接、重投影、时序合成及全变分（TV）去斑处理。

## 🚀 亮点 (Release Notes)

针对弱网环境和跨洋大规模数据下载（长时运行）进行了深度重构与 I/O 调优：

* ① **智能故障分类重构**：优先检查数据损坏（`_DATA_ERRORS`），避免对已损坏的 TIFF 切片进行无效的指数退避重试。
* ② **网络瞬态修复**：将 `vsicurl` 的瞬态假 404（由断连引起）纳入网络错误（`_NET_ERRORS`），杜绝误判跳过有效影像。
* ③ **消除双重等待**：彻底禁用 GDAL 内部重试（`MAX_RETRY=0`），将重试控制权完全交还 Python 层，避免单次失败引发长达 90s+ 的卡顿。
* ④ **激进连接策略**：将连接超时时间从 30s 缩短至 12s，实现“快速失败、快速重试”。
* ⑤ **熔断器机制 (Circuit Breaker)**：跨线程全局监控，连续 8 次网络失败自动触发熔断，挂起所有请求 90s 以等待网络环境恢复。
* ⑥ **会话后台保活**：针对 9 小时以上的长时任务，引入自动 Cookie 每 4 小时后台静默刷新机制，告别 401/403 授权断连。

---

## 🛠️ 核心工作流程 (Workflow)

1. **环境与认证初始化**：解析命令行参数，配置专为弱网优化的 GDAL C 层环境变量（增大缓存、关闭多路复用）。通过 `earthaccess` 验证身份建立持久化会话。
2. **统一空间网格构建**：根据用户输入的 ROI 边界和分辨率（默认30m），计算全局统一的仿射变换矩阵（Transform）和 AOI 掩膜。
3. **元数据分页检索**：按月分段调用 `asf_search` 检索时间范围内的 Burst 极化切片，并按采集日期分组。
4. **单日并行拉取与拼接 (Warp & Mosaic)**：多线程并发，利用 `/vsicurl/` 发起 HTTP Range 请求。**仅拉取 ROI 相交像素**，内存中直接重投影、像素累加求均值，落盘生成单日拼接图。
5. **时序合成与后处理 (Composite & TV & Clip)**：
* 按时间窗口（10天 / 12天 / 月度）进行**中值合成 (Median Composite)**。
* 基于 Bregman 迭代的 **TV 全变分去斑**（有效保持边缘并消除瞬态噪声）。
* 应用 AOI 掩膜将边界外像素置 NaN。
* 最后一次性写盘输出，大幅减少 I/O 消耗。



---

## ⚙️ 依赖环境 (Prerequisites)

请确保您的 Python 环境中安装了以下核心依赖：

```bash
pip install numpy pandas rasterio asf_search earthaccess scipy scikit-image shapely tqdm gdal

```

*(注：需要配置有 NASA Earthdata 的 `~/.netrc` 凭证用于鉴权)*

---

## 💻 运行与配置说明

所有用户配置通过命令行参数传入，推荐由 `Run_RTC.sh` 统一管理调用。


## 📁 产出目录结构

流水线运行完毕后，会在您指定的 `--base-dir` 下生成如下结构：

```text
📁 base_dir/
 └── 📁 downloads/
      ├── 📄 manifest.csv                    # 单日期拉取结果报告
      ├── 📄 2023-01-01_VH.tif               # 单日 Warp & Mosaic 结果 (可选保留)
      ├── 📄 ...
      └── 📁 composite/
           ├── 📄 manifest_comp.csv          # 时序合成处理报告
           ├── 📄 2023-01-01_VH_10d_TV_clip.tif  # 最终输出：10天合成+TV去斑+掩膜裁剪
           ├── 📄 2023-01-11_VH_10d_TV_clip.tif
           └── 📄 ...

```

最终输出的 `*_TV_clip.tif` 文件包含了经过清洗的高质量时间序列数据。

---
opera的支持RTC产品

链接：https://www.earthdata.nasa.gov/data/catalog/asf-opera-l2-rtc-s1-v1-1
