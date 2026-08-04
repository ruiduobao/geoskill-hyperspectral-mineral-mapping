# Hyperspectral Mineral Mapping (geoskill-hyperspectral-mineral-mapping)

> Identifies and maps mineral distribution from hyperspectral imagery using SAM spectral angle mapping and continuum removal, outputting mineral class/confidence rasters and a report

---

## 1. Overview

(Fill in 2-3 paragraphs of introduction here: functions, application scenarios, core algorithms.)

## 2. Features

(Fill in 2-3 paragraphs of introduction here: functions, application scenarios, core algorithms.)

## 3. Quick Start

```bash
pip install -r requirements.txt
python geoskill-hyperspectral-mineral-mapping.py --bbox 116 39 117 40 --synthetic --output-dir ./out
```

## 4. CLI Parameters

Run `python <skill>.py --help` for the full list. Common parameters:

| Parameter | Type | Description |
|---|---|---|
| `--bbox` | `float[4]` | WGS84 bounding box `min_lon min_lat max_lon max_lat` |
| `--input` | `path` | Local input file (GeoJSON/GeoTIFF/etc.) |
| `--output-dir` | `path` | Output directory (default `./output`) |
| `--synthetic` | `flag` | Use synthetic data instead of real input |
| `--quiet` | `flag` | Suppress non-essential stdout |

## 5. Input / Output

| File | Format | Description |
|---|---|---|
| `result.tif` | GeoTIFF | Primary output |
| `output-manifest.json` | JSON | Run manifest |


## 6. Technical Principle

(see SKILL.md for details)

## 7. Methodology

This skill has been methodologically reviewed. See [`REVIEW.md`](./REVIEW.md) for:

- P0/P1/P2 issue counts and verdicts
- Reproduction commands
- Known limitations and edge cases

## 8. License

MIT License. See [`LICENSE`](./LICENSE) for full text.

---

# 高光谱矿物制图（中文版）

> 本部分为完整中文文档，与上方英文部分对应。

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
