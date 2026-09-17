import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

from velocitykit.platforms.common import build_velocity_adata_from_anndata


def _adata(values, cells, genes):
    return AnnData(
        X=np.asarray(values),
        obs=pd.DataFrame(index=cells),
        var=pd.DataFrame(index=genes),
    )


def test_dual_subtraction_requires_exact_genes_and_aligns_order():
    total = _adata([[5, 7]], ["cell"], ["g1", "g2"])
    exonic = _adata([[2, 3]], ["cell"], ["g2", "g1"])

    result = build_velocity_adata_from_anndata(total, exonic)

    assert result.var_names.tolist() == ["g1", "g2"]
    np.testing.assert_array_equal(result.layers["spliced"], [[3, 2]])
    np.testing.assert_array_equal(result.layers["unspliced"], [[2, 5]])
    np.testing.assert_array_equal(
        result.layers["spliced"] + result.layers["unspliced"], total.X
    )
    assert result.uns["velocity_params"]["gene_join"] == "exact"


def test_dual_subtraction_rejects_gene_set_mismatch():
    total = _adata([[10]], ["cell"], ["total_only"])
    exonic = _adata([[5]], ["cell"], ["exonic_only"])

    with pytest.raises(ValueError, match="exactly the same gene identifiers"):
        build_velocity_adata_from_anndata(total, exonic)


def test_dual_subtraction_rejects_duplicate_gene_identifiers():
    total = _adata([[3, 2]], ["cell"], ["duplicate", "duplicate"])
    exonic = _adata([[1, 1]], ["cell"], ["duplicate", "duplicate"])

    with pytest.raises(ValueError, match="duplicate gene identifiers"):
        build_velocity_adata_from_anndata(total, exonic)


def test_dual_subtraction_requires_every_filtered_cell():
    total = _adata([[5], [7]], ["cell-a", "cell-b"], ["g1"])
    exonic = _adata([[3]], ["cell-a"], ["g1"])

    with pytest.raises(ValueError, match="missing 1 of 2 filtered"):
        build_velocity_adata_from_anndata(total, exonic)


def test_dual_subtraction_ignores_extra_raw_cells_and_preserves_total_order():
    total = _adata([[7], [5]], ["cell-b", "cell-a"], ["g1"])
    exonic = _adata([[3], [2], [99]], ["cell-a", "cell-b", "background"], ["g1"])

    result = build_velocity_adata_from_anndata(total, exonic)

    assert result.obs_names.tolist() == ["cell-b", "cell-a"]
    np.testing.assert_array_equal(result.layers["spliced"], [[2], [3]])
    np.testing.assert_array_equal(result.layers["unspliced"], [[5], [2]])
