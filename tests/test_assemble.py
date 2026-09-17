import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

from velocitykit.assemble import assemble_from_manifest


def _write_velocity(path, cell, genes, spliced, sample=None):
    obs = pd.DataFrame(index=[cell])
    if sample is not None:
        obs["sample"] = sample
    adata = AnnData(
        X=np.asarray([spliced]),
        obs=obs,
        var=pd.DataFrame(index=genes),
    )
    adata.layers["spliced"] = np.asarray([spliced])
    adata.layers["unspliced"] = np.asarray([spliced]) + 1
    adata.write_h5ad(path)


def test_assemble_exact_aligns_genes_and_adds_manifest_annotations(tmp_path):
    _write_velocity(tmp_path / "first.h5ad", "cell-a", ["g1", "g2"], [1, 2])
    _write_velocity(tmp_path / "second.h5ad", "cell-b", ["g2", "g1"], [20, 10])
    manifest = tmp_path / "samples.csv"
    pd.DataFrame(
        {
            "path": ["first.h5ad", "second.h5ad"],
            "source": ["first", "second"],
            "sample": ["E6", "E7"],
            "batch": ["B1", "B2"],
        }
    ).to_csv(manifest, index=False)

    combined = assemble_from_manifest(str(manifest))

    assert combined.shape == (2, 2)
    assert combined.var_names.tolist() == ["g1", "g2"]
    assert combined.obs_names.tolist() == ["cell-a:first", "cell-b:second"]
    assert combined.obs["sample"].astype(str).tolist() == ["E6", "E7"]
    assert combined.obs["batch"].astype(str).tolist() == ["B1", "B2"]
    np.testing.assert_array_equal(combined.X, [[1, 2], [10, 20]])
    np.testing.assert_array_equal(combined.layers["spliced"], combined.X)
    assert combined.uns["velocitykit_assembly"]["source_labels"].tolist() == [
        "first",
        "second",
    ]


def test_assemble_requires_exact_genes_by_default(tmp_path):
    _write_velocity(tmp_path / "first.h5ad", "a", ["g1", "g2"], [1, 2])
    _write_velocity(tmp_path / "second.h5ad", "b", ["g1", "g3"], [3, 4])
    manifest = tmp_path / "samples.csv"
    pd.DataFrame({"path": ["first.h5ad", "second.h5ad"]}).to_csv(
        manifest, index=False
    )

    with pytest.raises(ValueError, match="Gene set.*differs"):
        assemble_from_manifest(str(manifest))

    combined = assemble_from_manifest(str(manifest), gene_join="intersection")
    assert combined.var_names.tolist() == ["g1"]
    assert combined.shape == (2, 1)


def test_assemble_rejects_conflicting_manifest_annotations(tmp_path):
    _write_velocity(
        tmp_path / "input.h5ad", "cell", ["g1", "g2"], [1, 2], sample="E6"
    )
    manifest = tmp_path / "samples.csv"
    pd.DataFrame({"path": ["input.h5ad"], "sample": ["E7"]}).to_csv(
        manifest, index=False
    )

    with pytest.raises(ValueError, match="conflicts with existing"):
        assemble_from_manifest(str(manifest))


def test_assemble_rejects_duplicate_source_labels(tmp_path):
    _write_velocity(tmp_path / "first.h5ad", "a", ["g1"], [1])
    _write_velocity(tmp_path / "second.h5ad", "b", ["g1"], [2])
    manifest = tmp_path / "samples.csv"
    pd.DataFrame(
        {
            "path": ["first.h5ad", "second.h5ad"],
            "source": ["same", "same"],
        }
    ).to_csv(manifest, index=False)

    with pytest.raises(ValueError, match="source labels must be unique"):
        assemble_from_manifest(str(manifest))
