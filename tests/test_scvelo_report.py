import numpy as np
import pandas as pd
import pytest

pytest.importorskip("scanpy")
pytest.importorskip("scvelo")

from anndata import AnnData

from velocitykit.scvelo_report import (
    _plot_filename,
    attach_cell_metadata,
    calculate_input_qc_metrics,
    read_velocity_input,
    subset_cells,
)


def _adata():
    return AnnData(
        X=np.zeros((3, 2)),
        obs=pd.DataFrame(
            {
                "cell_id": ["cell-a", "cell-b", "cell-c"],
                "existing": [1, 1, 2],
            },
            index=["0", "1", "2"],
        ),
    )


@pytest.mark.parametrize("suffix", [".h5ad", ".loom"])
def test_read_velocity_input_accepts_h5ad_and_loom(tmp_path, suffix):
    adata = _adata()
    adata.layers["spliced"] = np.ones(adata.shape)
    adata.layers["unspliced"] = np.zeros(adata.shape)
    input_file = tmp_path / f"velocity{suffix}"
    if suffix == ".h5ad":
        adata.write_h5ad(input_file)
    else:
        adata.write_loom(input_file)

    loaded = read_velocity_input(str(input_file))

    assert loaded.shape == adata.shape
    assert {"spliced", "unspliced"}.issubset(loaded.layers)


def test_read_velocity_input_rejects_unsupported_or_missing_layers(tmp_path):
    unsupported = tmp_path / "velocity.txt"
    unsupported.write_text("not an AnnData file")
    with pytest.raises(ValueError, match="Unsupported input format"):
        read_velocity_input(str(unsupported))

    incomplete = tmp_path / "incomplete.h5ad"
    _adata().write_h5ad(incomplete)
    with pytest.raises(ValueError, match="missing required RNA-velocity layers"):
        read_velocity_input(str(incomplete))


def test_attach_cell_metadata_infers_unique_shared_key(tmp_path):
    metadata_file = tmp_path / "metadata.csv"
    pd.DataFrame(
        {
            "cell_id": ["cell-c", "cell-a", "cell-b", "extra-cell"],
            "sample": ["E9", "E6", "E7", "unused"],
            "existing": [2, 1, 1, 4],
        }
    ).to_csv(metadata_file, index=False)

    adata = _adata()
    details = attach_cell_metadata(adata, str(metadata_file))

    assert details["metadata_key"] == "cell_id"
    assert details["adata_key"] == "cell_id"
    assert details["matched_cells"] == 3
    assert details["extra_metadata_rows"] == 1
    assert adata.obs["sample"].astype(str).tolist() == ["E6", "E7", "E9"]


def test_attach_cell_metadata_can_join_to_obs_index(tmp_path):
    metadata_file = tmp_path / "metadata.tsv"
    pd.DataFrame(
        {"barcode": ["2", "0", "1"], "cell_type": ["C", "A", "B"]}
    ).to_csv(metadata_file, sep="\t", index=False)

    adata = _adata()
    details = attach_cell_metadata(
        adata,
        str(metadata_file),
        metadata_key="barcode",
        adata_key="_index",
    )

    assert details["adata_key"] == "_index"
    assert adata.obs["cell_type"].astype(str).tolist() == ["A", "B", "C"]


def test_attach_cell_metadata_rejects_duplicate_ids(tmp_path):
    metadata_file = tmp_path / "metadata.csv"
    pd.DataFrame(
        {"cell_id": ["cell-a", "cell-a", "cell-c"], "sample": ["A", "B", "C"]}
    ).to_csv(metadata_file, index=False)

    with pytest.raises(ValueError, match="duplicate cell identifiers"):
        attach_cell_metadata(
            _adata(), str(metadata_file), metadata_key="cell_id", adata_key="cell_id"
        )


def test_attach_cell_metadata_rejects_missing_cells(tmp_path):
    metadata_file = tmp_path / "metadata.csv"
    pd.DataFrame(
        {"cell_id": ["cell-a", "cell-c"], "sample": ["A", "C"]}
    ).to_csv(metadata_file, index=False)

    with pytest.raises(ValueError, match="missing 1 of 3 cells"):
        attach_cell_metadata(
            _adata(), str(metadata_file), metadata_key="cell_id", adata_key="cell_id"
        )


def test_attach_cell_metadata_rejects_conflicting_existing_values(tmp_path):
    metadata_file = tmp_path / "metadata.csv"
    pd.DataFrame(
        {"cell_id": ["cell-a", "cell-b", "cell-c"], "existing": [1, 99, 2]}
    ).to_csv(metadata_file, index=False)

    with pytest.raises(ValueError, match="conflicts with existing"):
        attach_cell_metadata(
            _adata(), str(metadata_file), metadata_key="cell_id", adata_key="cell_id"
        )


def test_plot_filename_is_safe():
    assert _plot_filename("cell type / broad") == "umap_cell_type_broad.png"


def test_subset_cells_retains_requested_values():
    adata = _adata()
    adata.obs["sample"] = pd.Categorical(["E6", "E7", "E6"])

    selected, details = subset_cells(adata, "sample", ["E6"])

    assert selected.n_obs == 2
    assert selected.obs_names.tolist() == ["0", "2"]
    assert details == {
        "column": "sample",
        "values": ["E6"],
        "original_cells": 3,
        "retained_cells": 2,
    }
    assert selected.uns["velocitykit_cell_subset"] == details


def test_subset_cells_accepts_multiple_values_and_removes_duplicates():
    adata = _adata()
    adata.obs["sample"] = ["E6", "E7", "E8"]

    selected, details = subset_cells(adata, "sample", ["E8", "E6", "E8"])

    assert selected.obs_names.tolist() == ["0", "2"]
    assert details["values"] == ["E8", "E6"]


def test_subset_cells_rejects_unknown_column_or_value():
    adata = _adata()
    adata.obs["sample"] = ["E6", "E7", "E8"]

    with pytest.raises(ValueError, match="Subset column 'missing'"):
        subset_cells(adata, "missing", ["E6"])
    with pytest.raises(ValueError, match="Requested values were not found"):
        subset_cells(adata, "sample", ["E9"])


def test_input_qc_metrics_are_calculated_before_gene_filtering():
    adata = AnnData(
        X=np.array(
            [
                [1, 0, 2, 0],
                [0, 3, 4, 5],
                [0, 0, 0, 6],
            ]
        )
    )

    calculate_input_qc_metrics(adata)
    expected_genes = [2, 3, 1]
    expected_counts = [3, 12, 6]

    # Simulate downstream HVG selection. Observation-level input QC is retained.
    adata = adata[:, :2].copy()
    assert adata.obs["n_genes_by_counts"].tolist() == expected_genes
    assert adata.obs["total_counts"].tolist() == expected_counts
    assert adata.uns["velocitykit_input_qc"] == {
        "matrix": "X",
        "n_genes_before_filtering": 4,
    }
