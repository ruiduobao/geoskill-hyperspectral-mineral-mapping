"""Core algorithm tests for hyperspectral-mineral-mapping."""
import sys
import os

import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(__file__))
from conftest import mod


class TestLibrary:
    def test_library_shape(self):
        for lib in mod.MINERAL_LIBS:
            wl, names, specs = mod.mineral_library(lib)
            assert wl.size == specs.shape[1]
            assert len(names) == specs.shape[0]
            assert np.all(specs > 0)

    def test_unknown_library_raises(self):
        with pytest.raises(mod.UsageError):
            mod.mineral_library("unknown_lib")

    def test_spectra_have_absorption(self):
        # 每种矿物在其特征吸收处反射率应低于连续统两端
        _, _, specs = mod.mineral_library("usgs")
        wl = mod.default_wavelengths()
        for k in range(specs.shape[0]):
            spec = specs[k]
            assert spec.min() < spec.max()  # 有吸收起伏


class TestContinuumRemoval:
    def test_output_le_one(self):
        wl = mod.default_wavelengths()
        spec = mod._gaussian_dips(wl, 0.6, 0.1, [(2.2, 0.2, 0.03)])
        cr = mod.continuum_removal(spec)
        assert cr.max() <= 1.0 + 1e-6
        # 吸收谷处连续统去除值应明显 < 1
        i_abs = int(np.argmin(spec))
        assert cr[i_abs] < 0.9

    def test_flat_spectrum_returns_one(self):
        spec = np.full(20, 0.5)
        cr = mod.continuum_removal(spec)
        np.testing.assert_allclose(cr, 1.0, atol=1e-3)

    def test_cube_shape(self):
        wl = mod.default_wavelengths()
        spec = mod._gaussian_dips(wl, 0.6, 0.1, [(2.2, 0.2, 0.03)])
        cube = np.repeat(spec[:, None, None], 8, axis=1)
        cube = np.repeat(cube, 8, axis=2)
        cr = mod.continuum_removal_cube(cube)
        assert cr.shape == cube.shape


class TestSAM:
    def test_self_angle_zero(self):
        wl, names, lib = mod.mineral_library("usgs")
        # 库中第 0 条光谱做单像元 → 对其自身夹角应为 0
        cube = lib[0][:, None, None]
        angles = mod.sam_angles(cube, lib)
        assert angles[0, 0, 0] == pytest.approx(0.0, abs=1e-5)
        assert angles[1, 0, 0] > 0.01

    def test_classify_selects_correct(self):
        wl, names, lib = mod.mineral_library("usgs")
        # lib 形状 (n_minerals, bands)；立方体需要 (bands, H, W)
        cube = lib.T[:, :, None]   # (bands, n_min, 1)
        idx, best, conf = mod.classify_minerals(cube, lib)
        assert list(idx[:, 0].astype(int)) == [0, 1, 2]
        assert np.all(best < 1e-4)
        assert np.all(conf >= 0)


class TestSyntheticAccuracy:
    def test_classification_matches_truth(self):
        bbox = [116.0, 39.0, 117.0, 40.0]
        for lib in mod.MINERAL_LIBS:
            cube, wl, truth, info = mod.generate_synthetic(bbox, lib, seed=3)
            _, _, lib_specs = mod.mineral_library(lib)
            idx, best, conf = mod.classify_minerals(cube, lib_specs)
            acc = float(np.mean(idx.astype(int) == truth.astype(int)))
            assert acc > 0.9, f"library {lib} accuracy {acc} too low"

    def test_continuum_removal_keeps_accuracy(self):
        bbox = [116.0, 39.0, 117.0, 40.0]
        cube, wl, truth, info = mod.generate_synthetic(bbox, "usgs", seed=5)
        _, _, lib = mod.mineral_library("usgs")
        cube_cr = mod.continuum_removal_cube(cube)
        lib_cr = np.stack([mod.continuum_removal(lib[k])
                           for k in range(lib.shape[0])], axis=0).astype(np.float32)
        idx, best, conf = mod.classify_minerals(cube_cr, lib_cr)
        acc = float(np.mean(idx.astype(int) == truth.astype(int)))
        assert acc > 0.9


class TestGeoTiff:
    def test_roundtrip(self, tmp_path):
        arr = np.random.uniform(0, 2, (16, 16)).astype(np.float32)
        path = str(tmp_path / "t.tif")
        mod.write_geotiff(path, arr, [116.0, 39.0, 117.0, 40.0])
        back, bb = mod.read_cube(path)
        assert back.shape == (1, 16, 16)
        np.testing.assert_allclose(bb, [116.0, 39.0, 117.0, 40.0], atol=1e-6)
