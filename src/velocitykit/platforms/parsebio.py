"""Parse Biosciences Split Pipe support for velocity-kit.

Split Pipe records one row per assigned transcript in
``process/tscp_assignment.csv.gz``.  Its ``exonic`` column therefore provides
the spliced/unspliced classification directly; unlike the 10x and PIPseq
workflows, Parse data do not require a second count run or subtraction.
"""

import logging
import re
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

logger = logging.getLogger(__name__)

TRANSCRIPT_FILENAME = "tscp_assignment.csv.gz"
METADATA_RELATIVE_PATH = Path("all-sample/DGE_filtered/cell_metadata.csv")


def add_arguments(parser):
    """Add Parse Biosciences arguments to an argparse parser."""
    sublibrary_group = parser.add_mutually_exclusive_group(required=True)
    sublibrary_group.add_argument(
        "--sublibrary",
        action="append",
        help=(
            "Split Pipe sublibrary output directory. Repeat once per sublibrary. "
            "Each directory must contain process/tscp_assignment.csv.gz and "
            "all-sample/DGE_filtered/cell_metadata.csv."
        ),
    )
    sublibrary_group.add_argument(
        "--sublibraries-dir",
        help=(
            "Parent directory containing relocated sublibrary result folders. "
            "VelocityKit reads their base folder names and suffix order from "
            "the combined Split Pipe log."
        ),
    )
    parser.add_argument(
        "--combined-metadata",
        required=True,
        help=(
            "Combined Split Pipe output directory (preferred), or its "
            "cell_metadata.csv. The output directory lets VelocityKit read the "
            "combine log and reproduce Split Pipe's exact __sN mapping."
        ),
    )
    parser.add_argument(
        "--gene-column",
        choices=("gene_name", "gene"),
        default="gene_name",
        help=(
            "Transcript-assignment column to use as var_names: gene_name for "
            "symbols (default), or gene for stable gene IDs."
        ),
    )
    parser.add_argument(
        "--chunk-size",
        type=int,
        default=1_000_000,
        help="Rows read from each transcript table at a time (default: 1000000).",
    )
    parser.add_argument(
        "--out-h5ad",
        help="Output .h5ad file path (optional if --out-loom is provided).",
    )
    parser.add_argument(
        "--out-loom",
        help="Output .loom file path (optional if --out-h5ad is provided).",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="count",
        default=1,
        help="Increase verbosity level. Use -v for INFO, -vv for DEBUG.",
    )


def resolve_sublibrary_inputs(sublibrary: Path) -> Tuple[Path, Path]:
    """Resolve and validate the transcript and filtered-cell files."""
    sublibrary = Path(sublibrary)
    if not sublibrary.is_dir():
        raise NotADirectoryError(f"Parse sublibrary is not a directory: {sublibrary}")

    transcript_path = sublibrary / "process" / TRANSCRIPT_FILENAME
    metadata_path = sublibrary / METADATA_RELATIVE_PATH
    missing = [path for path in (transcript_path, metadata_path) if not path.is_file()]
    if missing:
        formatted = "\n".join(f"  - {path}" for path in missing)
        raise FileNotFoundError(
            f"Missing required Parse Split Pipe file(s):\n{formatted}"
        )
    return transcript_path, metadata_path


def resolve_combined_metadata(path: Path) -> Path:
    """Accept either a combined output directory or cell_metadata.csv itself."""
    path = Path(path)
    if path.is_dir():
        path = path / METADATA_RELATIVE_PATH
    if not path.is_file():
        raise FileNotFoundError(f"Combined cell metadata not found: {path}")
    return path


def _find_combine_log(reference: Path, metadata_path: Path) -> Optional[Path]:
    """Locate the unique Split Pipe log associated with combined metadata."""
    reference = Path(reference)
    if reference.is_dir():
        combined_root = reference
    elif metadata_path.parent.name == "DGE_filtered":
        # <combined>/all-sample/DGE_filtered/cell_metadata.csv
        try:
            combined_root = metadata_path.parents[2]
        except IndexError:
            combined_root = metadata_path.parent
    else:
        combined_root = metadata_path.parent

    logs = sorted(combined_root.glob("split-pipe_v*.log"))
    if not logs:
        return None
    if len(logs) > 1:
        formatted = "\n".join(f"  - {path}" for path in logs)
        raise ValueError(
            "Multiple Split Pipe logs were found for the combined output; "
            f"cannot choose an authoritative log:\n{formatted}"
        )
    return logs[0]


def _parse_combine_sublibraries(log_path: Path) -> Dict[Path, str]:
    """Map resolved sublibrary paths to Split Pipe's numeric suffix labels."""
    pattern = re.compile(
        r"^# sublibraries\s+Item_(\d+)_of_\d+\s+(.+?)\s*$", re.MULTILINE
    )
    entries = pattern.findall(log_path.read_text(errors="replace"))
    if not entries:
        raise ValueError(
            f"Could not find ordered '# sublibraries Item_N_of_M' entries in "
            f"Split Pipe log: {log_path}"
        )

    mapping: Dict[Path, str] = {}
    for position, raw_path in entries:
        path = Path(raw_path).expanduser().resolve()
        if path in mapping:
            raise ValueError(f"Duplicate sublibrary path in {log_path}: {path}")
        mapping[path] = position
    return mapping


def resolve_sublibrary_locations(
    sublibraries: Optional[Sequence[str]],
    sublibraries_dir: Optional[Path],
    combined_reference: Path,
) -> List[Path]:
    """Resolve explicit inputs or relocated folders named by the combine log."""
    if sublibraries:
        return [Path(path) for path in sublibraries]
    if sublibraries_dir is None:
        raise ValueError("Specify --sublibrary or --sublibraries-dir")

    root = Path(sublibraries_dir)
    if not root.is_dir():
        raise NotADirectoryError(f"--sublibraries-dir is not a directory: {root}")
    combined_metadata_path = resolve_combined_metadata(combined_reference)
    log_path = _find_combine_log(combined_reference, combined_metadata_path)
    if log_path is None:
        raise ValueError(
            "--sublibraries-dir requires the combined output directory with its "
            "split-pipe_v*.log; a standalone cell_metadata.csv is insufficient"
        )

    log_mapping = _parse_combine_sublibraries(log_path)
    ordered = sorted(log_mapping.items(), key=lambda item: int(item[1]))
    base_names = [path.name for path, _ in ordered]
    if len(set(base_names)) != len(base_names):
        raise ValueError(
            f"Sublibrary base folder names are not unique in {log_path}: {base_names}"
        )

    resolved = [root / base_name for base_name in base_names]
    missing = [path for path in resolved if not path.is_dir()]
    if missing:
        formatted = "\n".join(f"  - {path}" for path in missing)
        raise FileNotFoundError(
            "Could not find logged sublibrary base folder(s) under "
            f"{root}:\n{formatted}"
        )
    logger.info(
        "Resolved relocated Parse sublibraries from %s using %s",
        root,
        log_path,
    )
    return resolved


def _split_combined_metadata(
    combined_metadata: pd.DataFrame,
) -> Dict[str, pd.DataFrame]:
    """Partition combined metadata by suffix and restore raw barcodes."""
    extracted = combined_metadata["bc_wells"].str.extract(
        r"^(?P<base>.+)__s(?P<label>\d+)$"
    )
    invalid = extracted.isna().any(axis=1)
    if invalid.any():
        examples = combined_metadata.loc[invalid, "bc_wells"].head(5).tolist()
        raise ValueError(
            "Combined metadata contains cell barcodes without a numeric __sN "
            f"suffix (examples: {examples})"
        )

    partitions = {}
    for label in sorted(extracted["label"].unique(), key=int):
        mask = extracted["label"] == label
        partition = combined_metadata.loc[mask].copy()
        partition["bc_wells"] = extracted.loc[mask, "base"].to_numpy()
        if partition["bc_wells"].duplicated().any():
            raise ValueError(
                f"Combined metadata has duplicate raw barcodes within __s{label}"
            )
        partition = partition.set_index("bc_wells", drop=False)
        partition.index.name = None
        partitions[str(label)] = partition
    return partitions


def _metadata_matches_partition(
    local_metadata: pd.DataFrame,
    partition: pd.DataFrame,
) -> bool:
    """Compare all local metadata fields after removing a combined suffix."""
    if len(local_metadata) != len(partition):
        return False
    if set(local_metadata.index) != set(partition.index):
        return False
    if not set(local_metadata.columns).issubset(partition.columns):
        return False

    columns = list(local_metadata.columns)
    left = local_metadata.sort_index()[columns].fillna("<NA>").astype(str)
    right = partition.sort_index()[columns].fillna("<NA>").astype(str)
    return left.equals(right)


def _infer_labels_from_metadata(
    local_metadata: Sequence[pd.DataFrame],
    partitions: Dict[str, pd.DataFrame],
) -> List[str]:
    """Find the unique one-to-one local-to-combined metadata assignment."""
    candidates = []
    for metadata in local_metadata:
        matches = [
            label
            for label, partition in partitions.items()
            if _metadata_matches_partition(metadata, partition)
        ]
        candidates.append(matches)

    solutions: List[List[str]] = []

    def search(position: int, used: set, labels: List[str]) -> None:
        if len(solutions) > 1:
            return
        if position == len(candidates):
            solutions.append(labels.copy())
            return
        for label in candidates[position]:
            if label not in used:
                used.add(label)
                labels.append(label)
                search(position + 1, used, labels)
                labels.pop()
                used.remove(label)

    search(0, set(), [])
    if len(solutions) == 1:
        return solutions[0]

    details = ", ".join(
        f"input {i}: {['__s' + label for label in labels] or 'no matches'}"
        for i, labels in enumerate(candidates, start=1)
    )
    if not solutions:
        raise ValueError(
            "Could not match every sublibrary to combined cell metadata exactly; "
            f"{details}"
        )
    raise ValueError(
        "Sublibrary-to-suffix mapping is ambiguous; refusing to guess. "
        f"Candidate mappings: {details}"
    )


def resolve_sublibrary_labels(
    sublibrary_paths: Sequence[Path],
    metadata_paths: Sequence[Path],
    combined_reference: Path,
) -> List[str]:
    """Resolve authoritative suffix labels from the combine log and metadata."""
    if len(sublibrary_paths) != len(metadata_paths):
        raise ValueError("Every sublibrary must have one local metadata file")
    if not sublibrary_paths:
        raise ValueError("At least one Parse sublibrary is required")

    resolved_paths = [Path(path).expanduser().resolve() for path in sublibrary_paths]
    if len(set(resolved_paths)) != len(resolved_paths):
        raise ValueError("The same Parse sublibrary was supplied more than once")

    combined_metadata_path = resolve_combined_metadata(combined_reference)
    combined_metadata = _read_cell_metadata(combined_metadata_path)
    partitions = _split_combined_metadata(combined_metadata)
    local_metadata = [_read_cell_metadata(path) for path in metadata_paths]
    log_path = _find_combine_log(combined_reference, combined_metadata_path)

    if log_path is not None:
        log_mapping = _parse_combine_sublibraries(log_path)
        labels = None
        if all(path in log_mapping for path in resolved_paths):
            labels = [log_mapping[path] for path in resolved_paths]
        else:
            logged_by_name: Dict[str, str] = {}
            duplicate_names = set()
            for logged_path, label in log_mapping.items():
                if logged_path.name in logged_by_name:
                    duplicate_names.add(logged_path.name)
                logged_by_name[logged_path.name] = label
            if not duplicate_names and all(
                path.name in logged_by_name for path in resolved_paths
            ):
                labels = [logged_by_name[path.name] for path in resolved_paths]
                logger.info(
                    "Matched relocated Parse sublibraries to the combine log "
                    "using base folder names"
                )

        if labels is not None:
            if len(set(labels)) != len(labels):
                raise ValueError(
                    f"Split Pipe log maps multiple inputs to one suffix: {log_path}"
                )
            for path, metadata, label in zip(resolved_paths, local_metadata, labels):
                partition = partitions.get(label)
                if partition is None or not _metadata_matches_partition(
                    metadata, partition
                ):
                    raise ValueError(
                        f"Split Pipe log maps {path} to __s{label}, but its local "
                        "cell metadata does not exactly match that combined partition"
                    )
            logger.info("Resolved Parse suffixes from Split Pipe log: %s", log_path)
            return labels
        logger.warning(
            "Combine log paths or unique base folder names do not match every "
            "supplied sublibrary; "
            "inferring suffixes from exact metadata matches instead"
        )

    labels = _infer_labels_from_metadata(local_metadata, partitions)
    logger.info("Resolved Parse suffixes from exact combined-metadata matches")
    return labels


def _read_cell_metadata(path: Path) -> pd.DataFrame:
    metadata = pd.read_csv(path, dtype={"bc_wells": str})
    if "bc_wells" not in metadata.columns:
        raise ValueError(f"Missing required 'bc_wells' column in {path}")
    if metadata["bc_wells"].isna().any():
        raise ValueError(f"Null cell barcodes found in {path}")
    if metadata["bc_wells"].duplicated().any():
        duplicates = metadata.loc[metadata["bc_wells"].duplicated(), "bc_wells"].head(5)
        raise ValueError(
            f"Duplicate cell barcodes found in {path}: {duplicates.tolist()}"
        )
    metadata = metadata.set_index("bc_wells", drop=False)
    metadata.index = metadata.index.astype(str)
    metadata.index.name = None
    return metadata


def _parse_exonic(values: pd.Series, transcript_path: Path) -> pd.Series:
    """Normalize Split Pipe boolean values and reject ambiguous classifications."""
    if pd.api.types.is_bool_dtype(values.dtype):
        return values.astype(bool)

    normalized = values.astype(str).str.strip().str.lower()
    mapped = normalized.map({"true": True, "false": False, "1": True, "0": False})
    if mapped.isna().any():
        invalid = sorted(normalized[mapped.isna()].unique().tolist())[:5]
        raise ValueError(
            f"Invalid values in the 'exonic' column of {transcript_path}: {invalid}"
        )
    return mapped.astype(bool)


def _coo_parts_to_csr(
    rows: List[np.ndarray],
    columns: List[np.ndarray],
    values: List[np.ndarray],
    shape: Tuple[int, int],
) -> sp.csr_matrix:
    if not rows:
        return sp.csr_matrix(shape, dtype=np.int64)
    matrix = sp.coo_matrix(
        (
            np.concatenate(values),
            (np.concatenate(rows), np.concatenate(columns)),
        ),
        shape=shape,
        dtype=np.int64,
    )
    # Repeated cell/gene entries from different input chunks are summed here.
    return matrix.tocsr()


def build_parse_sublibrary_adata(
    transcript_path: Path,
    metadata_path: Path,
    gene_column: str = "gene_name",
    chunk_size: int = 1_000_000,
) -> ad.AnnData:
    """Build a velocity AnnData object from one Split Pipe sublibrary.

    The filtered Split Pipe metadata supplies the authoritative called-cell set.
    Every retained row of ``tscp_assignment.csv.gz`` contributes one transcript,
    matching the counting method in Parse Biosciences' scVelo tutorial. The
    table's ``count`` field is deliberately not used because it represents read
    support rather than an additional transcript count.
    """
    transcript_path = Path(transcript_path)
    metadata_path = Path(metadata_path)
    if gene_column not in {"gene", "gene_name"}:
        raise ValueError("gene_column must be either 'gene_name' or 'gene'")
    if chunk_size <= 0:
        raise ValueError("chunk_size must be greater than zero")
    if not transcript_path.is_file():
        raise FileNotFoundError(
            f"Transcript assignment file not found: {transcript_path}"
        )
    if not metadata_path.is_file():
        raise FileNotFoundError(f"Cell metadata file not found: {metadata_path}")

    obs = _read_cell_metadata(metadata_path)
    cell_to_index = {barcode: i for i, barcode in enumerate(obs.index)}
    gene_to_index = {}
    genes: List[str] = []
    gene_ids = {}
    gene_names = {}

    spliced_rows: List[np.ndarray] = []
    spliced_columns: List[np.ndarray] = []
    spliced_values: List[np.ndarray] = []
    unspliced_rows: List[np.ndarray] = []
    unspliced_columns: List[np.ndarray] = []
    unspliced_values: List[np.ndarray] = []
    retained_transcripts = 0

    required_columns = ["bc_wells", "gene", "gene_name", "exonic"]
    logger.info(
        "Reading %s in chunks of %s rows for %s called cells",
        transcript_path,
        f"{chunk_size:,}",
        f"{len(obs):,}",
    )

    try:
        chunks = pd.read_csv(
            transcript_path,
            usecols=required_columns,
            dtype={"bc_wells": str, "gene": str, "gene_name": str},
            chunksize=chunk_size,
        )
        for chunk_number, chunk in enumerate(chunks, start=1):
            chunk = chunk.loc[chunk["bc_wells"].isin(cell_to_index)].copy()
            alternate_gene_column = (
                "gene" if gene_column == "gene_name" else "gene_name"
            )
            chunk["_feature_key"] = chunk[gene_column].fillna(
                chunk[alternate_gene_column]
            )
            chunk = chunk.loc[chunk["_feature_key"].notna()].copy()
            if chunk.empty:
                continue

            chunk["_feature_key"] = chunk["_feature_key"].astype(str)
            chunk["exonic"] = _parse_exonic(chunk["exonic"], transcript_path)

            unique_features = chunk.drop_duplicates("_feature_key")
            for key, gene_id, gene_name in unique_features[
                ["_feature_key", "gene", "gene_name"]
            ].itertuples(index=False, name=None):
                if key not in gene_to_index:
                    gene_to_index[key] = len(genes)
                    genes.append(key)
                    gene_ids[key] = gene_id if pd.notna(gene_id) else key
                    gene_names[key] = gene_name if pd.notna(gene_name) else key

            chunk["_cell_index"] = chunk["bc_wells"].map(cell_to_index)
            chunk["_gene_index"] = chunk["_feature_key"].map(gene_to_index)
            grouped = (
                chunk.groupby(
                    ["_cell_index", "_gene_index", "exonic"],
                    sort=False,
                    observed=True,
                )
                .size()
                .reset_index(name="transcripts")
            )
            retained_transcripts += len(chunk)

            for is_exonic, row_parts, column_parts, value_parts in (
                (True, spliced_rows, spliced_columns, spliced_values),
                (False, unspliced_rows, unspliced_columns, unspliced_values),
            ):
                selected = grouped.loc[grouped["exonic"] == is_exonic]
                if not selected.empty:
                    row_parts.append(selected["_cell_index"].to_numpy(dtype=np.int64))
                    column_parts.append(
                        selected["_gene_index"].to_numpy(dtype=np.int64)
                    )
                    value_parts.append(selected["transcripts"].to_numpy(dtype=np.int64))

            logger.debug(
                "Processed chunk %s; retained %s called-cell transcripts so far",
                chunk_number,
                f"{retained_transcripts:,}",
            )
    except ValueError as exc:
        if "Usecols do not match columns" in str(exc):
            raise ValueError(
                f"{transcript_path} is not a supported Split Pipe transcript "
                f"assignment table; required columns are {required_columns}"
            ) from exc
        raise

    if not genes:
        raise ValueError(
            f"No transcripts matched the {len(obs):,} filtered cell barcodes in "
            f"{transcript_path}"
        )

    shape = (len(obs), len(genes))
    spliced = _coo_parts_to_csr(spliced_rows, spliced_columns, spliced_values, shape)
    unspliced = _coo_parts_to_csr(
        unspliced_rows, unspliced_columns, unspliced_values, shape
    )
    total = spliced + unspliced

    var = pd.DataFrame(
        {
            "gene": [gene_ids[key] for key in genes],
            "gene_name": [gene_names[key] for key in genes],
        },
        index=pd.Index(genes),
    )
    var.index = var.index.astype(str)
    var.index.name = None

    adata = ad.AnnData(X=total, obs=obs.copy(), var=var)
    adata.layers["spliced"] = spliced
    adata.layers["unspliced"] = unspliced
    adata.obs["barcodes"] = adata.obs_names.astype(str)
    adata.var_names_make_unique()
    adata.uns["velocity_params"] = {
        "method": "parse_transcript_assignment",
        "gene_column": gene_column,
        "n_called_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "n_transcripts": int(retained_transcripts),
        "n_spliced": int(spliced.sum()),
        "n_unspliced": int(unspliced.sum()),
    }
    logger.info(
        "Built Parse sublibrary: %s cells x %s genes; "
        "%s spliced and %s unspliced transcripts",
        f"{adata.n_obs:,}",
        f"{adata.n_vars:,}",
        f"{int(spliced.sum()):,}",
        f"{int(unspliced.sum()):,}",
    )
    return adata


def combine_parse_sublibraries(
    adatas: Sequence[ad.AnnData],
    combined_metadata_path: Optional[Path] = None,
    labels: Optional[Sequence[str]] = None,
) -> ad.AnnData:
    """Union genes and append Parse's ``__sN`` sublibrary barcode suffixes."""
    if not adatas:
        raise ValueError("At least one Parse sublibrary is required")

    if labels is None:
        labels = [str(i) for i in range(1, len(adatas) + 1)]
    else:
        labels = [str(label) for label in labels]
    if len(labels) != len(adatas):
        raise ValueError("Every Parse sublibrary must have exactly one suffix label")
    if len(set(labels)) != len(labels):
        raise ValueError("Parse sublibrary suffix labels must be unique")
    if any(not label.isdigit() for label in labels):
        raise ValueError("Parse sublibrary suffix labels must be numeric")

    combined = ad.concat(
        adatas,
        join="outer",
        label="sublibrary",
        keys=labels,
        index_unique="__s",
        merge="first",
        uns_merge="unique",
        fill_value=0,
    )
    # ``anndata.concat(merge="first")`` can leave annotations empty for genes
    # found only in later sublibraries. Rebuild var from the first occurrence
    # of every gene so gene IDs and symbols remain complete after an outer join.
    var = pd.concat([adata.var for adata in adatas], axis=0)
    var = var.loc[~var.index.duplicated(keep="first")]
    combined.var = var.reindex(combined.var_names).copy()
    combined.obs["sublibrary"] = combined.obs["sublibrary"].astype(str)

    if combined_metadata_path is not None:
        combined_metadata_path = resolve_combined_metadata(combined_metadata_path)
        combined_metadata = _read_cell_metadata(combined_metadata_path)
        missing = combined.obs_names.difference(combined_metadata.index)
        extra = combined_metadata.index.difference(combined.obs_names)
        if len(missing):
            raise ValueError(
                f"Combined metadata is missing {len(missing):,} generated cells "
                f"(examples: {missing[:5].tolist()})"
            )
        if len(extra):
            raise ValueError(
                f"Combined metadata contains {len(extra):,} cells not generated "
                f"from the supplied sublibraries (examples: {extra[:5].tolist()})"
            )
        # Match Split Pipe combine-mode cell order even when CLI inputs were
        # intentionally supplied in another order.
        combined = combined[combined_metadata.index, :].copy()
        aligned = combined_metadata.reindex(combined.obs_names)
        for column in aligned.columns:
            combined.obs[column] = aligned[column].to_numpy()

    spliced = combined.layers["spliced"]
    unspliced = combined.layers["unspliced"]
    combined.uns["velocity_params"] = {
        "method": "parse_transcript_assignment",
        "n_sublibraries": len(adatas),
        "n_called_cells": int(combined.n_obs),
        "n_genes": int(combined.n_vars),
        "n_spliced": int(spliced.sum()),
        "n_unspliced": int(unspliced.sum()),
    }
    return combined


def _write_outputs(
    adata: ad.AnnData,
    out_h5ad: Optional[str],
    out_loom: Optional[str],
) -> None:
    if not out_h5ad and not out_loom:
        raise ValueError(
            "At least one output format must be specified: --out-h5ad or --out-loom"
        )
    if out_h5ad:
        path = Path(out_h5ad)
        path.parent.mkdir(parents=True, exist_ok=True)
        logger.info("Writing H5AD to %s", path)
        adata.write_h5ad(path, compression="gzip")
    if out_loom:
        path = Path(out_loom)
        path.parent.mkdir(parents=True, exist_ok=True)
        logger.info("Writing loom to %s", path)
        adata.write_loom(path)


def run(args):
    """Run the Parse Biosciences preparation pipeline."""
    if not args.out_h5ad and not args.out_loom:
        raise ValueError(
            "At least one output format must be specified: --out-h5ad or --out-loom"
        )
    if args.chunk_size <= 0:
        raise ValueError("--chunk-size must be greater than zero")

    sublibrary_locations = resolve_sublibrary_locations(
        args.sublibrary,
        Path(args.sublibraries_dir) if args.sublibraries_dir else None,
        Path(args.combined_metadata),
    )

    resolved_inputs = []
    for sublibrary in sublibrary_locations:
        transcript_path, metadata_path = resolve_sublibrary_inputs(Path(sublibrary))
        resolved_inputs.append((Path(sublibrary), transcript_path, metadata_path))

    labels = resolve_sublibrary_labels(
        [item[0] for item in resolved_inputs],
        [item[2] for item in resolved_inputs],
        Path(args.combined_metadata),
    )

    adatas = []
    for i, (sublibrary, transcript_path, metadata_path) in enumerate(
        resolved_inputs, start=1
    ):
        logger.info("Processing Parse sublibrary %s: %s", i, sublibrary)
        adatas.append(
            build_parse_sublibrary_adata(
                transcript_path,
                metadata_path,
                gene_column=args.gene_column,
                chunk_size=args.chunk_size,
            )
        )

    combined = combine_parse_sublibraries(
        adatas,
        Path(args.combined_metadata),
        labels=labels,
    )
    _write_outputs(combined, args.out_h5ad, args.out_loom)
    logger.info(
        "Successfully built Parse velocity data: %s cells x %s genes",
        f"{combined.n_obs:,}",
        f"{combined.n_vars:,}",
    )
