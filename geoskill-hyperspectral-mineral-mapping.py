#!/usr/bin/env python3
"""hyperspectral-mineral-mapping — 高光谱矿物制图

从高光谱影像（bands, H, W）识别并绘制矿物分布。核心方法：

- **SAM（光谱角制图）**：逐像元计算其光谱与矿物参考光谱库中每条光谱的夹角，
  最小夹角对应的矿物即分类结果。光谱角对光照/增益不敏感，适合矿物识别。
- **连续统去除（Continuum Removal）**：对光谱求上凸包作为连续统，逐波段相除，
  把反射光谱归一化到 [0,1] 连续统之下，突出吸收谷的形状与位置，增强矿物
  诊断性吸收特征（如 Al-OH ~2.20µm、Mg-OH/碳酸盐 ~2.32µm、Fe-OH ~2.34µm）。

内置两套简化矿物参考光谱库（USGS 风格 / ASTER 风格），包含高岭石、绿泥石、
方解石三种在短波红外具诊断吸收谷的矿物。

数据源：本地高光谱 GeoTIFF（短波红外多波段），或使用 ``--synthetic`` 生成空间
上分区分布三种矿物（各带特征吸收谷）的模拟高光谱立方体。

隐私声明 / Privacy：默认完全离线，不发起网络请求，所有处理本地完成。

Usage:
    python hyperspectral-mineral-mapping.py --bbox 116 39 117 40 --library usgs --output-dir ./out
    python hyperspectral-mineral-mapping.py --input cube.tif --library aster --output-dir ./out

License: MIT
"""
from __future__ import annotations

import argparse
import datetime as _dt
import json
import os
import sys
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

VERSION = "1.0.0"
SKILL_NAME = "hyperspectral-mineral-mapping"

try:
    from _geoskill_core.errors import (
        GeoSkillError, UsageError, ValidationError, ProcessError, to_exit_code,
    )
    from _geoskill_core.manifest import OutputManifest, OutputFile
except ImportError:  # pragma: no cover
    class GeoSkillError(Exception):
        def __init__(self, message: str, code: int = 7, kind: str = "EGeo", **kw):
            super().__init__(message)
            self.message, self.code, self.kind = message, code, kind

    class UsageError(GeoSkillError):
        def __init__(self, m, **k): super().__init__(m, code=2, kind="EUsage", **k)

    class ValidationError(GeoSkillError):
        def __init__(self, m, **k): super().__init__(m, code=6, kind="EValidate", **k)

    class ProcessError(GeoSkillError):
        def __init__(self, m, **k): super().__init__(m, code=7, kind="EProcess", **k)

    def to_exit_code(exc):
        return getattr(exc, "code", 7)

    OutputManifest = None
    OutputFile = None


def _utc_now() -> str:
    return _dt.datetime.now(_dt.timezone.utc).isoformat().replace("+00:00", "Z")


# ---------------------------------------------------------------------------
# 内置矿物参考光谱库（简化 USGS / ASTER 风格）
# 波长范围 0.4–2.5 µm；诊断吸收谷以高斯凹陷叠加在连续统上。
# ---------------------------------------------------------------------------
WL_MIN, WL_MAX, WL_NB = 0.40, 2.50, 43   # 43 波段，步长 0.05 µm


def _gaussian_dips(wl: np.ndarray, base: float, slope: float,
                   features: List[Tuple[float, float, float]]) -> np.ndarray:
    """连续统(base + slope*(wl-2)) 上叠加若干高斯吸收谷 (center, depth, width)。"""
    continuum = base + slope * (wl - 2.0)
    spec = continuum.copy()
    for (c, depth, width) in features:
        spec = spec - depth * np.exp(-((wl - c) / width) ** 2)
    return np.clip(spec, 0.01, 1.0)


# 各库：矿物 → 吸收谷参数 [(center_um, depth, width), ...]
MINERAL_LIBS: Dict[str, Dict[str, Dict[str, Any]]] = {
    "usgs": {
        "kaolinite": {
            "label": "高岭石",
            "features": [(2.165, 0.20, 0.025), (2.205, 0.25, 0.028)],
            "base": 0.62, "slope": 0.12,
        },
        "chlorite": {
            "label": "绿泥石",
            "features": [(2.315, 0.28, 0.032), (2.255, 0.06, 0.030)],
            "base": 0.45, "slope": 0.08,
        },
        "calcite": {
            "label": "方解石",
            "features": [(2.340, 0.22, 0.030), (1.875, 0.16, 0.045)],
            "base": 0.70, "slope": 0.10,
        },
    },
    "aster": {
        "kaolinite": {
            "label": "高岭石",
            "features": [(2.170, 0.18, 0.030), (2.210, 0.23, 0.032)],
            "base": 0.60, "slope": 0.14,
        },
        "chlorite": {
            "label": "绿泥石",
            "features": [(2.320, 0.26, 0.036), (2.260, 0.05, 0.034)],
            "base": 0.47, "slope": 0.10,
        },
        "calcite": {
            "label": "方解石",
            "features": [(2.345, 0.20, 0.034), (1.880, 0.15, 0.050)],
            "base": 0.72, "slope": 0.12,
        },
    },
}


def default_wavelengths(n: int = WL_NB) -> np.ndarray:
    return np.linspace(WL_MIN, WL_MAX, int(n), dtype=float)


def mineral_library(library: str) -> Tuple[np.ndarray, List[str], np.ndarray]:
    """返回 (wavelengths, mineral_names, library_spectra[n_minerals, n_bands])。"""
    if library not in MINERAL_LIBS:
        raise UsageError(
            f"unknown library '{library}'. Choose from: {sorted(MINERAL_LIBS)}",
            library=library,
        )
    wl = default_wavelengths()
    names: List[str] = []
    specs: List[np.ndarray] = []
    for name, params in MINERAL_LIBS[library].items():
        names.append(name)
        specs.append(_gaussian_dips(wl, params["base"], params["slope"], params["features"]))
    return wl, names, np.stack(specs, axis=0).astype(np.float32)


# ---------------------------------------------------------------------------
# 连续统去除
# ---------------------------------------------------------------------------
def _upper_hull(xs: np.ndarray, ys: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Andrew 单调链算法求上凸包（始终 ≥ 所有点的凹包络）。"""
    order = np.argsort(xs)
    xs = xs[order]
    ys = ys[order]
    hull: List[int] = []
    for i in range(len(xs)):
        while len(hull) >= 2:
            x1, y1 = xs[hull[-2]], ys[hull[-2]]
            x2, y2 = xs[hull[-1]], ys[hull[-1]]
            cross = (x2 - x1) * (ys[i] - y1) - (y2 - y1) * (xs[i] - x1)
            if cross >= 0:      # 左转或共线 → 弹出（保留右拐的上包络）
                hull.pop()
            else:
                break
        hull.append(i)
    return xs[hull], ys[hull]


def continuum_removal(spectrum: np.ndarray) -> np.ndarray:
    """对单条光谱做连续统去除：反射率 / 上凸包连续统，输出 ≤1，突出吸收谷。

    连续统取光谱的上凸包（始终在光谱之上），因此比值恒 ≤1，吸收谷处明显 <1。
    """
    spec = np.asarray(spectrum, dtype=float)
    n = spec.size
    if n < 4:
        return np.ones_like(spec)
    x = np.arange(n, dtype=float)
    hx, hy = _upper_hull(x, spec)
    continuum = np.interp(x, hx, hy)
    continuum = np.maximum(continuum, 1e-6)
    return np.clip(spec / continuum, 0.0, 1.0)


def continuum_removal_cube(cube: np.ndarray) -> np.ndarray:
    """对 (bands,H,W) 立方体逐像元连续统去除。"""
    nb, h, w = cube.shape
    out = np.empty_like(cube, dtype=np.float32)
    flat = cube.reshape(nb, -1)
    res = np.empty_like(flat, dtype=np.float32)
    for i in range(flat.shape[1]):
        res[:, i] = continuum_removal(flat[:, i])
    return res.reshape(nb, h, w)


# ---------------------------------------------------------------------------
# SAM 分类
# ---------------------------------------------------------------------------
def sam_angles(cube: np.ndarray, library_spectra: np.ndarray) -> np.ndarray:
    """逐像元对库中每条光谱的夹角 (弧度)。返回 (n_minerals, H, W)。

    NoData (NaN) 像元的夹角输出为 π/2（与所有库光谱正交，无法判定）。
    """
    nb, h, w = cube.shape
    n_min = library_spectra.shape[0]
    pix = cube.reshape(nb, -1).astype(float)            # (bands, N)
    valid_pix = np.all(np.isfinite(pix), axis=0)         # (N,)
    pix_norm = pix / (np.linalg.norm(pix, axis=0, keepdims=True) + 1e-12)
    lib = library_spectra.astype(float)
    lib_norm = lib / (np.linalg.norm(lib, axis=1, keepdims=True) + 1e-12)
    cosang = lib_norm @ pix_norm                         # (n_min, N)
    cosang = np.clip(cosang, -1.0, 1.0)
    angles = np.arccos(cosang)                           # (n_min, N)
    # NoData 像元 → 角 = π/2 (cosang = 0)
    angles[:, ~valid_pix] = np.pi / 2.0
    return angles.reshape(n_min, h, w)


def classify_minerals(
    cube: np.ndarray, library_spectra: np.ndarray,
    valid_mask: Optional[np.ndarray] = None,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """SAM 矿物分类。返回 (index_map[H,W], best_angle[H,W], confidence[H,W])。

    confidence = (second_angle - best_angle) / second_angle，越大越可靠。
    valid_mask (H, W) if provided, NoData 像素的 index_map 标为 -1。
    """
    angles = sam_angles(cube, library_spectra)           # (n_min, H, W)
    sorted_ang = np.sort(angles, axis=0)
    best = sorted_ang[0]
    second = sorted_ang[1] if angles.shape[0] > 1 else sorted_ang[0]
    index_map = np.argmin(angles, axis=0).astype(np.int32)
    confidence = np.zeros_like(best)
    np.divide(second - best, second, out=confidence, where=second > 1e-9)
    confidence = np.clip(confidence, 0.0, 1.0).astype(np.float32)
    if valid_mask is not None:
        index_map = np.where(valid_mask, index_map, -1).astype(np.int32)
        # NoData 区的 confidence 强制为 0
        confidence = np.where(valid_mask, confidence, 0.0).astype(np.float32)
    return index_map, best.astype(np.float32), confidence


# ---------------------------------------------------------------------------
# 合成数据：三种矿物空间分区分布
# ---------------------------------------------------------------------------
def generate_synthetic(
    bbox: List[float], library: str,
    width: int = 64, height: int = 64, seed: int = 42,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, Dict[str, Any]]:
    """生成 (bands,H,W) 高光谱立方体 + 波长 + 真值矿物索引图。

    竖向三等分：上=矿物0、中=矿物1、下=矿物2，各带噪声。
    """
    rng = np.random.default_rng(seed)
    wl, names, lib_specs = mineral_library(library)
    nb = wl.size
    n_min = len(names)
    cube = np.zeros((nb, height, width), dtype=np.float32)
    truth = np.zeros((height, width), dtype=np.float32)

    row_bands = height // n_min
    for k in range(n_min):
        r0 = k * row_bands
        r1 = height if k == n_min - 1 else (k + 1) * row_bands
        block = np.repeat(lib_specs[k][:, None], width, axis=1)   # (nb, width)
        for r in range(r0, r1):
            cube[:, r, :] = block
        truth[r0:r1, :] = float(k)

    cube = cube + rng.normal(0, 0.004, size=cube.shape).astype(np.float32)
    cube = np.clip(cube, 0.01, 1.0).astype(np.float32)

    info = {
        "bbox": bbox, "width": width, "height": height, "n_bands": nb,
        "wavelengths_um": wl.tolist(),
        "minerals": names,
        "row_bands": int(row_bands),
        "library": library,
    }
    return cube, wl, truth, info


# ---------------------------------------------------------------------------
# GeoTIFF I/O
# ---------------------------------------------------------------------------
def write_geotiff(path, array, bbox, nodata=-9999.0):
    import rasterio
    from rasterio.transform import from_bounds

    if array.ndim == 2:
        array = array[np.newaxis, ...]
    nb, h, w = array.shape
    transform = from_bounds(bbox[0], bbox[1], bbox[2], bbox[3], w, h)
    profile = {
        "driver": "GTiff", "height": h, "width": w, "count": nb,
        "dtype": "float32", "crs": "EPSG:4326", "transform": transform,
        "nodata": nodata, "compress": "deflate",
    }
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with rasterio.open(path, "w", **profile) as dst:
        for b in range(nb):
            dst.write(array[b].astype("float32"), b + 1)


def read_cube(path: str) -> Tuple[np.ndarray, List[float]]:
    import rasterio
    if not os.path.exists(path):
        raise UsageError(f"input raster not found: {path}", path=path)
    with rasterio.open(path) as src:
        cube = src.read().astype(np.float32)
        b = src.bounds
        bbox = [b.left, b.bottom, b.right, b.top]
    return cube, bbox


def read_cube_with_nodata(path: str) -> Tuple[np.ndarray, List[float], Optional[float]]:
    """Read multi-band cube and replace NoData pixels with NaN.

    A pixel is NoData if ANY band equals the nodata sentinel. Returns
    (cube (bands, H, W), bbox, nodata_value_or_None).
    """
    import rasterio
    if not os.path.exists(path):
        raise UsageError(f"input raster not found: {path}", path=path)
    with rasterio.open(path) as src:
        cube = src.read().astype(np.float32)
        b = src.bounds
        bbox = [b.left, b.bottom, b.right, b.top]
        nodata = src.nodata
    if nodata is not None and np.isfinite(nodata):
        bad_mask = np.any(cube == nodata, axis=0)
        cube[:, bad_mask] = np.nan
    return cube, bbox, nodata


def validate_bbox(bbox: Optional[List[float]], allow_none: bool = False) -> List[float]:
    """Validate a W,S,E,N bbox. Cross-180 / out-of-range / W>=E / S>=N -> ValidationError."""
    if bbox is None:
        if allow_none:
            return None  # type: ignore[return-value]
        raise ValidationError("bbox is required")
    if len(bbox) != 4:
        raise ValidationError(f"bbox must have 4 floats, got {len(bbox)}")
    w, s, e, n = bbox
    for v, name in zip([w, s, e, n], ["W", "S", "E", "N"]):
        if not isinstance(v, (int, float)) or not (-1e9 < v < 1e9):
            raise ValidationError(f"bbox {name}={v!r} not a finite number")
    if w == e or s == n:
        raise ValidationError(f"bbox has zero area: W={w} E={e} S={s} N={n}")
    if not (-180.0 <= w <= 180.0) or not (-180.0 <= e <= 180.0):
        raise ValidationError(f"bbox lon out of [-180,180]: W={w} E={e}")
    if not (-90.0 <= s <= 90.0) or not (-90.0 <= n <= 90.0):
        raise ValidationError(f"bbox lat out of [-90,90]: S={s} N={n}")
    if w > e:
        if not (w > 170.0 and e < -170.0):
            raise ValidationError(
                f"bbox has W>E (minLon > maxLon): W={w} E={e} — "
                f"if crossing the dateline, split into two bboxes (e.g. "
                f"[{w}, {s}, 180, {n}] and [-180, {s}, {e}, {n}])"
            )
        raise ValidationError(
            f"bbox crosses the 180° dateline (W={w} E={e}); "
            f"split into two non-wrapping bboxes ([{w}, {s}, 180, {n}] and "
            f"[-180, {s}, {e}, {n}]) and run separately"
        )
    if s > n:
        raise ValidationError(f"bbox has S>N (minLat > maxLat): S={s} N={n}")
    return [float(w), float(s), float(e), float(n)]


# ---------------------------------------------------------------------------
# Manifest
# ---------------------------------------------------------------------------
def write_manifest(output_dir, args, outputs, qa, started_at, exit_code, bbox):
    if OutputManifest is None:
        return None
    cmd = " ".join([SKILL_NAME] + sys.argv[1:])
    man = OutputManifest(
        skill=SKILL_NAME, skill_version=VERSION, command=cmd,
        started_at=started_at, finished_at=_utc_now(), exit_code=exit_code,
        inputs={
            "input": getattr(args, "input", None),
            "library": getattr(args, "library", None),
            "continuum_removal": bool(getattr(args, "continuum_removal", False)),
            "synthetic": bool(getattr(args, "synthetic", False)),
            "bbox": bbox,
        },
        outputs=[OutputFile(**o) for o in outputs],
        qa=qa,
        software={"python": sys.version.split()[0], "skill": SKILL_NAME},
    )
    path = os.path.join(output_dir, "output-manifest.json")
    with open(path, "w", encoding="utf-8") as f:
        json.dump(man.to_dict(), f, ensure_ascii=False, indent=2, default=str)
    return path


# ---------------------------------------------------------------------------
# 主流程
# ---------------------------------------------------------------------------
def process(args: argparse.Namespace) -> int:
    started_at = _utc_now()
    output_dir = args.output_dir

    # ---- 1. 参数验证 (前置：失败不创建 output_dir) ----
    bbox = list(args.bbox) if args.bbox else None
    if bbox is not None:
        bbox = validate_bbox(bbox)
    wl, names, lib_specs = mineral_library(args.library)

    # ---- 2. 数据获取 ----
    synth_info: Optional[Dict[str, Any]] = None
    truth: Optional[np.ndarray] = None
    input_nodata: Optional[float] = None
    valid_mask: Optional[np.ndarray] = None
    n_valid_input: int = 0
    n_total_input: int = 0
    if args.input and not args.synthetic:
        cube, file_bbox, input_nodata = read_cube_with_nodata(args.input)
        bbox = bbox if bbox is not None else file_bbox
        bbox = validate_bbox(bbox)
        valid_mask = np.all(np.isfinite(cube), axis=0)
        n_valid_input = int(valid_mask.sum())
        n_total_input = int(cube.shape[1] * cube.shape[2])
        if n_valid_input == 0:
            raise ValidationError(
                f"input cube has no valid (non-NoData) pixels "
                f"(nodata={input_nodata}, total={n_total_input})"
            )
        source_note = args.input
    else:
        if bbox is None:
            raise UsageError("provide --bbox (synthetic mode) or --input <raster>")
        cube, wl, truth, synth_info = generate_synthetic(bbox, args.library)
        source_note = "synthetic"
        n_valid_input = int(cube.size)
        n_total_input = int(cube.shape[1] * cube.shape[2])

    # ---- 3. 校验通过后创建 output_dir ----
    os.makedirs(output_dir, exist_ok=True)

    if cube.ndim != 3 or cube.size == 0:
        raise ValidationError("input must be a 3D (bands,H,W) hyperspectral cube")

    # 波段数与库对齐（真实输入可能波段数不同 → 重采样库到输入波段数）
    if cube.shape[0] != wl.size:
        wl_in = default_wavelengths(cube.shape[0])
        lib_resampled = np.stack([
            np.interp(wl_in, wl, lib_specs[k]) for k in range(lib_specs.shape[0])
        ], axis=0).astype(np.float32)
        wl = wl_in
    else:
        lib_resampled = lib_specs

    work_cube = cube
    if args.continuum_removal:
        work_cube = continuum_removal_cube(cube)
        lib_resampled = np.stack([
            continuum_removal(lib_resampled[k]) for k in range(lib_resampled.shape[0])
        ], axis=0).astype(np.float32)

    index_map, best_angle, confidence = classify_minerals(work_cube, lib_resampled, valid_mask=valid_mask)

    # 输出
    mineral_path = os.path.join(output_dir, "mineral_map.tif")
    write_geotiff(mineral_path, index_map.astype(np.float32), bbox, nodata=-1.0)
    conf_path = os.path.join(output_dir, "confidence.tif")
    write_geotiff(conf_path, confidence, bbox, nodata=0.0)

    # 统计 abundance 仅基于 valid 像素；class_counts 包含所有唯一值（含 NoData=-1）
    uniq, ucounts = np.unique(index_map, return_counts=True)
    counts = {int(k): int(v) for k, v in zip(uniq, ucounts)}
    if valid_mask is not None:
        idx_for_count = index_map[valid_mask]
    else:
        idx_for_count = index_map
    abundance = {}
    n_valid_for_abundance = int(idx_for_count.size)
    for k, v in counts.items():
        if k == -1:
            continue
        if k < len(names):
            abundance[names[k]] = float(v) / n_valid_for_abundance
    if -1 in counts:
        abundance["nodata"] = float(counts[-1]) / n_valid_for_abundance

    report = {
        "source": source_note, "library": args.library,
        "continuum_removal": bool(args.continuum_removal),
        "n_bands": int(cube.shape[0]), "shape": [int(cube.shape[1]), int(cube.shape[2])],
        "minerals": names,
        "class_counts": counts,
        "abundance": abundance,
        "mean_best_angle_deg": float(np.degrees(np.mean(best_angle))),
        "mean_confidence": float(np.mean(confidence)),
    }
    if truth is not None:
        accuracy = float(np.mean(index_map.astype(int) == truth.astype(int)))
        report["synthetic_overall_accuracy"] = accuracy

    report_path = os.path.join(output_dir, "mineral_report.json")
    with open(report_path, "w", encoding="utf-8") as f:
        json.dump(report, f, ensure_ascii=False, indent=2)

    qa: Dict[str, Any] = {
        "source": source_note, "library": args.library,
        "n_bands": int(cube.shape[0]),
        "mean_confidence": float(np.mean(confidence)),
        "mean_best_angle_deg": float(np.degrees(np.mean(best_angle))),
        "abundance": abundance,
        "n_valid_pixels": int(n_valid_input),
        "n_total_pixels": int(n_total_input),
        "input_nodata": input_nodata,
    }
    if truth is not None:
        qa["synthetic_overall_accuracy"] = report["synthetic_overall_accuracy"]

    outputs = [
        {"path": mineral_path, "kind": "raster", "crs_epsg": 4326,
         "bbox_wgs84": bbox, "band_count": 1, "nodata": -1.0},
        {"path": conf_path, "kind": "raster", "crs_epsg": 4326,
         "bbox_wgs84": bbox, "band_count": 1, "nodata": 0.0},
        {"path": report_path, "kind": "json"},
    ]
    man_path = write_manifest(output_dir, args, outputs, qa, started_at, 0, bbox)

    if not args.quiet:
        print(f"[{SKILL_NAME}] source: {source_note}  library: {args.library}  CR: {args.continuum_removal}")
        print(f"[{SKILL_NAME}] bands: {cube.shape[0]}  shape: {cube.shape[1:]}")
        print(f"[{SKILL_NAME}] abundance: {abundance}")
        print(f"[{SKILL_NAME}] mean confidence: {qa['mean_confidence']:.3f}  mean angle: {qa['mean_best_angle_deg']:.2f}°")
        if truth is not None:
            print(f"[{SKILL_NAME}] synthetic accuracy: {report['synthetic_overall_accuracy']:.3f}")
        print(f"[{SKILL_NAME}] output: {mineral_path}")
        if man_path:
            print(f"[{SKILL_NAME}] manifest: {man_path}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog=SKILL_NAME,
        description="Hyperspectral mineral mapping via SAM with optional continuum removal.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--bbox", nargs=4, type=float, metavar=("W", "S", "E", "N"),
                   help="geographic extent [minLon minLat maxLon maxLat]")
    p.add_argument("--input", help="input hyperspectral GeoTIFF (bands,H,W)")
    p.add_argument("--library", default="usgs", choices=sorted(MINERAL_LIBS.keys()),
                   help="reference spectral library (default: usgs)")
    p.add_argument("--continuum-removal", dest="continuum_removal", action="store_true",
                   help="apply continuum removal before SAM classification")
    p.add_argument("--synthetic", action="store_true",
                   help="generate a synthetic mineral scene (offline)")
    p.add_argument("--output-dir", default="./output", help="output directory")
    p.add_argument("--quiet", action="store_true", help="suppress console output")
    p.add_argument("--version", action="version", version=f"%(prog)s {VERSION}")
    return p


def main(argv: Optional[List[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return process(args)
    except GeoSkillError as exc:
        print(f"[{SKILL_NAME}] ERROR [{exc.kind}] {exc.message}", file=sys.stderr)
        return exc.code
    except KeyboardInterrupt:
        return 130
    except Exception as exc:  # noqa: BLE001
        print(f"[{SKILL_NAME}] ERROR {exc}", file=sys.stderr)
        return to_exit_code(exc)


if __name__ == "__main__":
    sys.exit(main())
