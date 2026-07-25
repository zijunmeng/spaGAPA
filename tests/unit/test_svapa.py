"""Unit tests for the SVAPA module (permutation Moran's I)."""

import numpy as np
import pandas as pd
import pytest

from spagapa.analysis import morans_i, gearys_c, svapa
from spagapa.analysis.svapa import _build_spatial_weights


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------
@pytest.fixture
def grid_coords():
    """10x10 regular grid of spots."""
    return np.mgrid[0:10, 0:10].reshape(2, -1).T.astype(float)


@pytest.fixture
def clustered_values(grid_coords):
    """Strong x-axis gradient -> high positive spatial autocorrelation."""
    rng = np.random.default_rng(0)
    return grid_coords[:, 0] + rng.normal(0, 0.05, 100)


@pytest.fixture
def random_values():
    """IID noise -> near-zero spatial autocorrelation."""
    rng = np.random.default_rng(1)
    return rng.normal(0.5, 0.3, 100)


# ---------------------------------------------------------------------------
# morans_i
# ---------------------------------------------------------------------------
class TestMoransI:
    def test_clustered_high(self, grid_coords, clustered_values):
        """A clear spatial gradient should yield a high, significant Moran's I."""
        I, p, n = morans_i(clustered_values, grid_coords, k=8, n_perm=200, seed=0)
        assert n == 100
        assert I > 0.5
        assert p < 0.05

    def test_random_near_zero(self, grid_coords, random_values):
        """IID noise should give Moran's I near zero and non-significant."""
        I, p, n = morans_i(random_values, grid_coords, k=8, n_perm=200, seed=0)
        assert n == 100
        assert -0.3 < I < 0.3
        assert p > 0.1

    def test_clustered_higher_than_random(self, grid_coords, clustered_values,
                                          random_values):
        """Ordering invariant: clustered > random."""
        I_c, _, _ = morans_i(clustered_values, grid_coords, k=8, n_perm=50, seed=0)
        I_r, _, _ = morans_i(random_values, grid_coords, k=8, n_perm=50, seed=0)
        assert I_c > I_r

    def test_negative_autocorrelation(self, grid_coords):
        """A checkerboard pattern should give a negative Moran's I."""
        checker = np.array([((i + j) % 2) for i in range(10) for j in range(10)],
                           dtype=float)
        I, p, n = morans_i(checker, grid_coords, k=4, n_perm=200, seed=0)
        assert n == 100
        assert I < 0

    def test_length_mismatch_raises(self, grid_coords):
        with pytest.raises(ValueError):
            morans_i(np.zeros(50), grid_coords, k=8, n_perm=0)

    def test_n_perm_zero_skips_pvalue(self, grid_coords, clustered_values):
        I, p, n = morans_i(clustered_values, grid_coords, k=8, n_perm=0)
        assert np.isfinite(I)
        assert np.isnan(p)
        assert n == 100

    def test_returns_three_tuple(self, grid_coords, clustered_values):
        out = morans_i(clustered_values, grid_coords, k=8, n_perm=10, seed=0)
        assert isinstance(out, tuple)
        assert len(out) == 3


# ---------------------------------------------------------------------------
# NaN handling
# ---------------------------------------------------------------------------
class TestNanHandling:
    def test_nan_entries_dropped(self, grid_coords, clustered_values):
        """NaNs should be subset out (n < 100) without error."""
        v = clustered_values.copy()
        v[::3] = np.nan
        I, p, n = morans_i(v, grid_coords, k=6, n_perm=50, seed=0)
        assert n < 100
        assert np.isfinite(I)
        # spatial structure is preserved -> still significant / positive
        assert I > 0.3

    def test_all_nan_returns_nan(self, grid_coords):
        v = np.full(100, np.nan)
        I, p, n = morans_i(v, grid_coords, k=8, n_perm=10, seed=0)
        assert np.isnan(I)
        assert np.isnan(p)
        assert n == 0

    def test_too_few_obs_returns_nan(self, grid_coords):
        v = np.full(100, np.nan)
        v[:2] = [0.1, 0.2]  # only 2 finite -> < 3
        I, p, n = morans_i(v, grid_coords, k=8, n_perm=10, seed=0)
        assert np.isnan(I)
        assert n == 2

    def test_constant_values_returns_nan(self, grid_coords):
        """Zero variance -> z'z = 0 -> Moran's I undefined."""
        v = np.full(100, 0.5)
        I, p, n = morans_i(v, grid_coords, k=8, n_perm=10, seed=0)
        assert np.isnan(I)


# ---------------------------------------------------------------------------
# gearys_c
# ---------------------------------------------------------------------------
class TestGearysC:
    def test_clustered_below_one(self, grid_coords, clustered_values):
        C, p, n = gearys_c(clustered_values, grid_coords, k=8, n_perm=200, seed=0)
        assert n == 100
        assert C < 1.0  # positive autocorrelation
        assert p < 0.05

    def test_random_near_one(self, grid_coords, random_values):
        C, p, n = gearys_c(random_values, grid_coords, k=8, n_perm=200, seed=0)
        assert n == 100
        assert 0.6 < C < 1.4


# ---------------------------------------------------------------------------
# svapa matrix
# ---------------------------------------------------------------------------
class TestSvapa:
    def test_ranking_spatial_gene_first(self, grid_coords, clustered_values,
                                        random_values):
        """The spatial gene should rank above the noise gene."""
        M = np.vstack([random_values, clustered_values])
        res = svapa(M, grid_coords, k=8, n_perm=100, seed=0, min_obs=5)
        assert list(res.columns) == [
            "gene", "morans_i", "morans_i_pvalue", "gearys_c",
            "gearys_c_pvalue", "n_obs", "padj",
        ]
        # sorted by morans_i desc -> clustered (Gene_1) first
        assert res.iloc[0]["gene"] == "Gene_1"
        assert res.iloc[0]["morans_i"] > res.iloc[1]["morans_i"]

    def test_dataframe_index_preserved(self, grid_coords, clustered_values,
                                       random_values):
        M = np.vstack([clustered_values, random_values])
        df = pd.DataFrame(M, index=["SPA", "NOISE"])
        res = svapa(df, grid_coords, k=8, n_perm=50, seed=0, min_obs=5)
        assert set(res["gene"]) == {"SPA", "NOISE"}
        assert res.iloc[0]["gene"] == "SPA"

    def test_min_obs_filters_sparse_genes(self, grid_coords, clustered_values):
        """Genes below min_obs observed spots should be dropped."""
        sparse = np.full(100, np.nan)
        sparse[:5] = np.linspace(0, 1, 5)  # only 5 obs
        M = np.vstack([clustered_values, sparse])
        res = svapa(M, grid_coords, k=8, n_perm=20, seed=0, min_obs=10)
        # sparse gene dropped, only the dense one remains
        assert len(res) == 1
        assert res.iloc[0]["n_obs"] == 100

    def test_fdr_correction_runs(self, grid_coords):
        rng = np.random.default_rng(0)
        M = rng.normal(0.5, 0.3, (10, 100))
        # inject one strong spatial gene
        M[0] = grid_coords[:, 0]
        res = svapa(M, grid_coords, k=8, n_perm=50, seed=0)
        assert "padj" in res.columns
        assert (res["padj"] >= 0).all()
        assert (res["padj"] <= 1.0 + 1e-9).all()

    def test_coords_mismatch_raises(self, grid_coords):
        M = np.zeros((3, 100))
        bad = grid_coords[:50]
        with pytest.raises(ValueError):
            svapa(M, bad, k=8, n_perm=10, seed=0)

    def test_sorted_descending(self, grid_coords):
        rng = np.random.default_rng(2)
        M = rng.normal(0.5, 0.3, (5, 100))
        res = svapa(M, grid_coords, k=8, n_perm=20, seed=0, min_obs=5)
        morans = res["morans_i"].values
        assert np.all(np.diff(morans[~np.isnan(morans)]) <= 1e-9)


# ---------------------------------------------------------------------------
# weights helper
# ---------------------------------------------------------------------------
class TestSpatialWeights:
    def test_shape_and_rowsums(self, grid_coords):
        W = _build_spatial_weights(grid_coords, k=8, symmetric=True)
        assert W.shape == (100, 100)
        row_sums = np.asarray(W.sum(axis=1)).ravel()
        # row-normalised -> each row sums to 1
        assert np.allclose(row_sums, 1.0)

    def test_k_clamped(self, grid_coords):
        # k larger than n-1 should be clamped without error
        W = _build_spatial_weights(grid_coords[:5], k=99, symmetric=True)
        assert W.shape == (5, 5)
