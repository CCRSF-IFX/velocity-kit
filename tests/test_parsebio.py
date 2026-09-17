import gzip
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from velocitykit.platforms.parsebio import (
    build_parse_sublibrary_adata,
    combine_parse_sublibraries,
    resolve_sublibrary_labels,
    resolve_sublibrary_locations,
    resolve_sublibrary_inputs,
)

TRANSCRIPT_COLUMNS = [
    "bc_wells",
    "genome",
    "gene",
    "gene_name",
    "count",
    "exonic",
    "rt_type",
    "bc_bcis",
    "cell_barcode",
    "polyN",
]


def _write_sublibrary(root: Path, cells, transcripts):
    process = root / "process"
    dge = root / "all-sample" / "DGE_filtered"
    process.mkdir(parents=True)
    dge.mkdir(parents=True)

    metadata = pd.DataFrame(
        {
            "bc_wells": cells,
            "sample": [f"sample_{cell}" for cell in cells],
            "tscp_count": [1] * len(cells),
        }
    )
    metadata.to_csv(dge / "cell_metadata.csv", index=False)

    table = pd.DataFrame(transcripts, columns=TRANSCRIPT_COLUMNS)
    with gzip.open(process / "tscp_assignment.csv.gz", "wt") as handle:
        table.to_csv(handle, index=False)
    return root


def _row(barcode, gene_id, gene_name, exonic, count=1):
    return [
        barcode,
        "genome",
        gene_id,
        gene_name,
        count,
        exonic,
        "R",
        "1_1_1",
        "ACGT_ACGT_ACGT",
        "AAAAAAAAAA",
    ]


def test_build_parse_sublibrary_counts_rows_and_filters_cells(tmp_path):
    sublibrary = _write_sublibrary(
        tmp_path / "sub1",
        ["01_01_01", "01_01_02"],
        [
            _row("01_01_01", "ENSG1", "G1", True, count=20),
            _row("01_01_01", "ENSG1", "G1", False, count=30),
            _row("01_01_01", "ENSG2", "G2", True),
            _row("01_01_02", "ENSG1", "G1", True),
            _row("not_a_cell", "ENSG3", "G3", True),
        ],
    )
    transcript_path, metadata_path = resolve_sublibrary_inputs(sublibrary)

    adata = build_parse_sublibrary_adata(transcript_path, metadata_path, chunk_size=2)

    assert adata.obs_names.tolist() == ["01_01_01", "01_01_02"]
    assert adata.var_names.tolist() == ["G1", "G2"]
    np.testing.assert_array_equal(adata.layers["spliced"].toarray(), [[1, 1], [1, 0]])
    np.testing.assert_array_equal(adata.layers["unspliced"].toarray(), [[1, 0], [0, 0]])
    np.testing.assert_array_equal(adata.X.toarray(), [[2, 1], [1, 0]])
    assert adata.uns["velocity_params"]["n_transcripts"] == 4


def test_gene_id_mode_and_sublibrary_combination(tmp_path):
    sub1 = _write_sublibrary(
        tmp_path / "sub1",
        ["01_01_01"],
        [_row("01_01_01", "ENSG1", "G1", True)],
    )
    sub2 = _write_sublibrary(
        tmp_path / "sub2",
        ["01_01_01", "01_01_03"],
        [
            _row("01_01_01", "ENSG1", "G1", False),
            _row("01_01_03", "ENSG3", "G3", True),
        ],
    )
    adatas = []
    metadata_paths = []
    for sublibrary in (sub1, sub2):
        transcript_path, metadata_path = resolve_sublibrary_inputs(sublibrary)
        metadata_paths.append(metadata_path)
        adatas.append(
            build_parse_sublibrary_adata(
                transcript_path,
                metadata_path,
                gene_column="gene",
                chunk_size=1,
            )
        )

    combined_root = tmp_path / "combined"
    combined_dge = combined_root / "all-sample" / "DGE_filtered"
    combined_dge.mkdir(parents=True)
    combined_metadata = combined_dge / "cell_metadata.csv"
    pd.DataFrame(
        {
            "bc_wells": ["01_01_01__s1", "01_01_01__s2", "01_01_03__s2"],
            "sample": [
                "sample_01_01_01",
                "sample_01_01_01",
                "sample_01_01_03",
            ],
            "tscp_count": [1, 1, 1],
        }
    ).to_csv(combined_metadata, index=False)
    (combined_root / "split-pipe_v1_8_2.log").write_text(
        "# sublibraries              Item_1_of_2\t/old/location/sub1\n"
        "# sublibraries              Item_2_of_2\t/old/location/sub2\n"
    )

    relocated = resolve_sublibrary_locations(None, tmp_path, combined_root)
    assert relocated == [sub1, sub2]

    labels = resolve_sublibrary_labels(
        [sub2, sub1],
        [metadata_paths[1], metadata_paths[0]],
        combined_root,
    )
    assert labels == ["2", "1"]

    combined = combine_parse_sublibraries(
        [adatas[1], adatas[0]], combined_root, labels=labels
    )

    assert combined.obs_names.tolist() == [
        "01_01_01__s1",
        "01_01_01__s2",
        "01_01_03__s2",
    ]
    assert combined.var_names.tolist() == ["ENSG1", "ENSG3"]
    assert combined.var.loc["ENSG3", "gene_name"] == "G3"
    assert combined.obs["sample"].tolist() == [
        "sample_01_01_01",
        "sample_01_01_01",
        "sample_01_01_03",
    ]
    assert combined.obs["sublibrary"].tolist() == ["1", "2", "2"]
    np.testing.assert_array_equal(
        combined.layers["spliced"].toarray(), [[1, 0], [0, 0], [0, 1]]
    )
    np.testing.assert_array_equal(
        combined.layers["unspliced"].toarray(), [[0, 0], [1, 0], [0, 0]]
    )

    # With no log available, the exact metadata match finds the same mapping.
    standalone_metadata = tmp_path / "combined_cell_metadata.csv"
    pd.read_csv(combined_metadata).to_csv(standalone_metadata, index=False)
    inferred = resolve_sublibrary_labels(
        [sub2, sub1],
        [metadata_paths[1], metadata_paths[0]],
        standalone_metadata,
    )
    assert inferred == ["2", "1"]


def test_invalid_exonic_value_is_rejected(tmp_path):
    sublibrary = _write_sublibrary(
        tmp_path / "sub1",
        ["01_01_01"],
        [_row("01_01_01", "ENSG1", "G1", "ambiguous")],
    )
    transcript_path, metadata_path = resolve_sublibrary_inputs(sublibrary)

    with pytest.raises(ValueError, match="Invalid values.*exonic"):
        build_parse_sublibrary_adata(transcript_path, metadata_path)


def test_missing_gene_symbol_falls_back_to_gene_id(tmp_path):
    sublibrary = _write_sublibrary(
        tmp_path / "sub1",
        ["01_01_01"],
        [_row("01_01_01", "ENSG1", None, True)],
    )
    transcript_path, metadata_path = resolve_sublibrary_inputs(sublibrary)

    adata = build_parse_sublibrary_adata(transcript_path, metadata_path)

    assert adata.var_names.tolist() == ["ENSG1"]
    assert adata.var.loc["ENSG1", "gene"] == "ENSG1"
    assert adata.var.loc["ENSG1", "gene_name"] == "ENSG1"
    assert int(adata.X.sum()) == 1


def test_combined_metadata_must_include_every_generated_cell(tmp_path):
    sublibrary = _write_sublibrary(
        tmp_path / "sub1",
        ["01_01_01"],
        [_row("01_01_01", "ENSG1", "G1", True)],
    )
    transcript_path, metadata_path = resolve_sublibrary_inputs(sublibrary)
    adata = build_parse_sublibrary_adata(transcript_path, metadata_path)
    incomplete = tmp_path / "cell_metadata.csv"
    pd.DataFrame({"bc_wells": ["different__s1"]}).to_csv(incomplete, index=False)

    with pytest.raises(ValueError, match="missing 1 generated cells"):
        combine_parse_sublibraries([adata], incomplete)


def test_suffix_inference_rejects_ambiguous_metadata(tmp_path):
    sub1 = _write_sublibrary(
        tmp_path / "sub1",
        ["01_01_01"],
        [_row("01_01_01", "ENSG1", "G1", True)],
    )
    sub2 = _write_sublibrary(
        tmp_path / "sub2",
        ["01_01_01"],
        [_row("01_01_01", "ENSG1", "G1", True)],
    )
    combined_metadata = tmp_path / "ambiguous.csv"
    pd.DataFrame(
        {
            "bc_wells": ["01_01_01__s1", "01_01_01__s2"],
            "sample": ["sample_01_01_01", "sample_01_01_01"],
            "tscp_count": [1, 1],
        }
    ).to_csv(combined_metadata, index=False)
    metadata_paths = [resolve_sublibrary_inputs(path)[1] for path in (sub1, sub2)]

    with pytest.raises(ValueError, match="mapping is ambiguous"):
        resolve_sublibrary_labels([sub1, sub2], metadata_paths, combined_metadata)
