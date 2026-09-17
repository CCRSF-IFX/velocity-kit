import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

pytest.importorskip("scanpy")

from velocitykit.batch import hansen_combat_correct


def _velocity_adata():
    spliced = np.array(
        [
            [4, 0, 2, 1],
            [2, 1, 0, 3],
            [8, 2, 1, 0],
            [3, 4, 2, 1],
        ],
        dtype=float,
    )
    unspliced = np.array(
        [
            [1, 0, 2, 2],
            [2, 3, 0, 1],
            [2, 2, 3, 0],
            [1, 1, 2, 3],
        ],
        dtype=float,
    )
    adata = AnnData(
        X=spliced.copy(),
        obs=pd.DataFrame(
            {
                "batch": ["run1", "run1", "run2", "run2"],
                "condition": ["control", "treated", "control", "treated"],
            },
            index=[f"cell-{i}" for i in range(4)],
        ),
        var=pd.DataFrame(index=[f"gene-{i}" for i in range(4)]),
    )
    adata.layers["spliced"] = spliced
    adata.layers["unspliced"] = unspliced
    return adata


def test_hansen_combat_preserves_velocity_relationship(monkeypatch):
    adata = _velocity_adata()
    original_spliced = adata.layers["spliced"].copy()
    original_unspliced = adata.layers["unspliced"].copy()
    monkeypatch.setattr("velocitykit.batch.sc.pp.combat", lambda adata, **_: None)

    corrected = hansen_combat_correct(
        adata,
        batch_key="batch",
        preserve_keys=["condition"],
        target_sum=100,
        min_shared_counts=0,
        n_top_genes=0,
    )

    corrected_spliced = corrected.layers["spliced"]
    corrected_unspliced = corrected.layers["unspliced"]
    corrected_total = corrected_spliced + corrected_unspliced
    original_total = original_spliced + original_unspliced
    original_ratio = np.divide(
        original_spliced,
        original_total,
        out=np.zeros_like(original_spliced),
        where=original_total > 0,
    )
    corrected_ratio = np.divide(
        corrected_spliced,
        corrected_total,
        out=np.zeros_like(corrected_spliced),
        where=corrected_total > 0,
    )

    assert corrected.shape == adata.shape
    assert np.all(corrected_spliced >= 0)
    assert np.all(corrected_unspliced >= 0)
    np.testing.assert_allclose(corrected_total.sum(axis=1), 100, rtol=1e-6)
    np.testing.assert_allclose(
        corrected_ratio[original_total > 0],
        original_ratio[original_total > 0],
        atol=1e-7,
    )
    assert np.all(corrected_total[original_total == 0] == 0)
    np.testing.assert_allclose(corrected.X, np.log1p(corrected_spliced))
    assert corrected.obs["n_genes_by_counts"].tolist() == [3, 3, 3, 4]
    assert corrected.obs["total_counts"].tolist() == [7, 6, 11, 10]
    provenance = corrected.uns["velocitykit_batch_correction"]
    assert provenance["method"] == "hansen-combat"
    assert provenance["batch_key"] == "batch"
    assert provenance["preserve_keys"].tolist() == ["condition"]
    assert provenance["max_spliced_ratio_error"] < 1e-7


def test_hansen_combat_clips_negative_corrected_totals(monkeypatch):
    adata = _velocity_adata()

    def make_negative(combat_input, **_):
        combat_input.X[0, 0] = -1

    monkeypatch.setattr("velocitykit.batch.sc.pp.combat", make_negative)
    corrected = hansen_combat_correct(
        adata,
        batch_key="batch",
        min_shared_counts=0,
        n_top_genes=0,
    )

    assert corrected.layers["spliced"][0, 0] == 0
    assert corrected.layers["unspliced"][0, 0] == 0
    provenance = corrected.uns["velocitykit_batch_correction"]
    assert provenance["negative_totals_clipped"] == 1
    assert provenance["negative_observed_totals_clipped"] == 1
    assert provenance["negative_structural_zero_totals_clipped"] == 0


@pytest.mark.parametrize(
    "batch_key,preserve_keys,message",
    [
        ("missing", None, "was not found"),
        ("batch", ["batch"], "cannot also be a preserved"),
        ("batch", ["missing"], "was not found"),
    ],
)
def test_hansen_combat_validates_design(batch_key, preserve_keys, message):
    with pytest.raises(ValueError, match=message):
        hansen_combat_correct(
            _velocity_adata(),
            batch_key=batch_key,
            preserve_keys=preserve_keys,
            min_shared_counts=0,
            n_top_genes=0,
        )


def test_hansen_combat_requires_multiple_nontrivial_batches():
    one_batch = _velocity_adata()
    one_batch.obs["batch"] = "run1"
    with pytest.raises(ValueError, match="at least two"):
        hansen_combat_correct(
            one_batch, batch_key="batch", min_shared_counts=0, n_top_genes=0
        )

    singleton = _velocity_adata()
    singleton.obs["batch"] = ["run1", "run1", "run1", "run2"]
    with pytest.raises(ValueError, match="at least two cells"):
        hansen_combat_correct(
            singleton, batch_key="batch", min_shared_counts=0, n_top_genes=0
        )
