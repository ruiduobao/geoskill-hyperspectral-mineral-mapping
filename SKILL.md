---
name: geoskill-hyperspectral-mineral-mapping
description: '用SAM光谱角制图与连续统去除从高中光谱影像识别并绘制矿物分布，输出矿物类别/置信度栅格与报告'
---

# 高光谱矿物制图 | Hyperspectral Mineral Mapping

(Fill in 2-3 paragraphs of introduction here: features, application scenarios, and core algorithm.)

## Dependencies / 依赖

```bash
pip install 'numpy' 'rasterio' 'scipy'
```

## Usage / 使用方法

### Basic Usage

```bash
python geoskill-hyperspectral-mineral-mapping.py --bbox 116.0 39.0 117.0 40.0 [other options]
```

### Example 1 (Synthetic Data, Offline)

```bash
python geoskill-hyperspectral-mineral-mapping.py --bbox 116.0 39.0 117.0 40.0 --synthetic --output-dir ./out
```

### Examples 2-5

(Add at least 4 real-world usage examples.)

## Output / 输出

| File | Format | Description |
|---|---|---|
| `result.tif` | GeoTIFF | Main output |
| `output-manifest.json` | JSON | Run manifest |

## Data Source / 数据源 / Source

(Describe the data source: free satellite data / local input / synthetic.)

## Privacy / 隐私声明 / Privacy

- Runs offline by default; `--synthetic` mode requires no network access at all.
- All processing is performed locally; no user data is uploaded.

## License / License

MIT

---

<!-- ===== 中文原文 (Chinese Original) ===== -->

---
name: geoskill-hyperspectral-mineral-mapping
description: '用SAM光谱角制图与连续统去除从高中光谱影像识别并绘制矿物分布，输出矿物类别/置信度栅格与报告'
---

# 高光谱矿物制图 | Hyperspectral Mineral Mapping

（在此填写 2-3 段中文介绍：功能、应用场景、核心算法。）

## 依赖

```bash
pip install 'numpy' 'rasterio' 'scipy'
```

## 使用方法

### 基本用法

```bash
python geoskill-hyperspectral-mineral-mapping.py --bbox 116.0 39.0 117.0 40.0 [其他参数]
```

### 示例 1（合成数据，离线）

```bash
python geoskill-hyperspectral-mineral-mapping.py --bbox 116.0 39.0 117.0 40.0 --synthetic --output-dir ./out
```

### 示例 2-5

（补充至少 4 个真实用法示例。）

## 输出

| 文件 | 格式 | 说明 |
|---|---|---|
| `result.tif` | GeoTIFF | 主产物 |
| `output-manifest.json` | JSON | 运行清单 |

## 数据源 / Source

（说明数据来源：免费卫星数据 / 本地输入 / 合成。）

## 隐私声明 / Privacy

- 默认离线运行，`--synthetic` 模式完全无网络。
- 所有处理在本地完成，不上传用户数据。

## License

MIT
