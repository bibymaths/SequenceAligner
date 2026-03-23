"""Command line interface for the alignment analysis tool.

The CLI provides subcommands to perform analyses on global, local and
LCS alignment outputs individually or collectively. It also supports
comparing alignment methods and generating summary reports. The
primary entry point is the ``main()`` function defined at the end of
this module.
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any

import numpy as np
import pandas as pd

from . import (
    block_detection,
    dp_matrix,
    fasta_utils,
    file_inventory,
    path_utils,
)

from . import residue_profiles, substitution_analysis, plotting, comparison, summary

logger = logging.getLogger(__name__)


def configure_logging(log_file: Optional[Path], quiet: bool = False) -> None:
    """Configure logging for the CLI.

    Parameters
    ----------
    log_file : Path or None
        Path to a file where logs should be written. If None, logs are
        written to stderr.
    quiet : bool, optional
        If True, suppress informational messages (set level to WARNING).

    Returns
    -------
    None
        Logging is configured globally.
    """
    handlers = []
    level = logging.WARNING if quiet else logging.INFO
    formatter = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s")
    if log_file:
        log_file.parent.mkdir(parents=True, exist_ok=True)
        file_handler = logging.FileHandler(log_file)
        file_handler.setFormatter(formatter)
        handlers.append(file_handler)
    stream_handler = logging.StreamHandler()
    stream_handler.setFormatter(formatter)
    handlers.append(stream_handler)
    logging.basicConfig(level=level, handlers=handlers)


def parse_common_args(subparser: argparse.ArgumentParser) -> None:
    """Add common options to a subparser.

    Parameters
    ----------
    subparser : argparse.ArgumentParser
        The parser to which arguments are added.
    """
    subparser.add_argument(
        "--results-dir",
        type=Path,
        required=True,
        help="Directory containing alignment output files",
    )
    subparser.add_argument(
        "--outdir",
        type=Path,
        default=None,
        help="Directory to write outputs (default: results-dir)",
    )
    subparser.add_argument(
        "--prefix", type=str, default="alignment", help="Prefix for output files"
    )
    subparser.add_argument(
        "--overwrite", action="store_true", help="Overwrite existing output files"
    )
    subparser.add_argument(
        "--log-file", type=Path, default=None, help="Write logs to this file"
    )
    subparser.add_argument(
        "--blosum",
        type=str,
        default="blosum62",
        choices=["blosum62", "none"],
        help="Substitution matrix to use for similarity metrics",
    )
    subparser.add_argument(
        "--min-block-length",
        type=int,
        default=5,
        help="Minimum length of conserved block to report",
    )
    subparser.add_argument(
        "--identity-threshold",
        type=float,
        default=0.7,
        help="Threshold for classifying a block as high identity",
    )
    subparser.add_argument(
        "--similarity-threshold",
        type=float,
        default=0.8,
        help="Threshold for classifying a block as conservative",
    )
    subparser.add_argument(
        "--window",
        type=int,
        default=2,
        help="Window size for local support and gap proximity calculations",
    )
    subparser.add_argument(
        "--plot-dpi", type=int, default=150, help="Resolution (dpi) for plots"
    )
    subparser.add_argument(
        "--quiet", action="store_true", help="Suppress informational logs"
    )


def analyse_method(
    method: str,
    files: file_inventory.AlignmentFiles,
    substitution_matrix: Optional[Dict[str, Dict[str, float]]],
    outdir: Path,
    prefix: str,
    min_block_length: int,
    identity_threshold: float,
    similarity_threshold: float,
    window: int,
    plot_dpi: int,
    overwrite: bool,
) -> Tuple[Dict[str, Any], pd.DataFrame, pd.DataFrame]:
    """Perform analysis for a single alignment method.

    This function loads the necessary files for the specified method,
    computes statistics, conserved blocks, residue support and path
    metrics, writes relevant outputs to TSV/PNG files and returns
    objects for further use.

    Parameters
    ----------
    method : str
        One of ``"global"``, ``"local"``, or ``"lcs"``.
    files : AlignmentFiles
        Inventory of available files.
    substitution_matrix : dict or None
        Substitution matrix for computing similarity metrics.
    outdir : Path
        Base directory for outputs.
    prefix : str
        Prefix for output filenames.
    min_block_length, identity_threshold, similarity_threshold : see
        ``block_detection.detect_blocks_to_df``.
    window : int
        Window size for local support.
    plot_dpi : int
        DPI for plot outputs.
    overwrite : bool
        If False and outputs exist, the analysis is skipped.

    Returns
    -------
    tuple
        A tuple (method_results, residue_support_df_a, residue_support_df_b).
        ``method_results`` is a dict containing various results and
        statistics for the method; ``residue_support_df_a`` and
        ``residue_support_df_b`` are DataFrames of residue support for
        sequence A and B, respectively.
    """
    logger.info("Starting %s analysis", method)
    # Determine file paths
    if method == "global":
        aln_fasta = files.global_alignment
        dp_bin = files.global_dp_bin
        dp_txt = files.global_dp_txt
        path_file = files.global_path
        stats_file = files.global_stats
    elif method == "local":
        aln_fasta = files.local_alignment
        dp_bin = files.local_dp_bin
        dp_txt = files.local_dp_txt
        path_file = files.local_path
        stats_file = files.local_stats
    elif method == "lcs":
        aln_fasta = files.lcs_alignment
        dp_bin = files.lcs_dp_bin
        dp_txt = files.lcs_dp_txt
        path_file = files.lcs_path
        stats_file = None  # no stats for LCS
    else:
        raise ValueError(f"Unknown method: {method}")
    if aln_fasta is None:
        raise FileNotFoundError(f"Alignment FASTA missing for method {method}")
    # Output prefix
    out_prefix = f"{prefix}_{method}"
    # Prepare result container
    results: Dict[str, Any] = {}
    # Parse alignment FASTA
    seqs = fasta_utils.parse_alignment_fasta(aln_fasta)
    if len(seqs) < 2:
        raise ValueError(f"Expected two sequences in {aln_fasta}, found {len(seqs)}")
    ids = list(seqs.keys())
    seq_a_id, seq_b_id = ids[0], ids[1]
    seq_a = seqs[seq_a_id]
    seq_b = seqs[seq_b_id]
    results["sequence_ids"] = (seq_a_id, seq_b_id)
    # Build coordinate maps
    a_map, b_map = fasta_utils.build_coordinate_maps(seq_a, seq_b)
    results["a_map"] = a_map
    results["b_map"] = b_map
    # Alignment statistics
    stats = fasta_utils.compute_alignment_stats(
        seq_a, seq_b, substitution_matrix, similarity_threshold=0
    )
    results["alignment_stats"] = stats
    # Ungapped lengths
    len_a = stats["ungapped_length_a"]
    len_b = stats["ungapped_length_b"]
    # Load DP matrix
    shape = dp_matrix.infer_shape(len_a, len_b)
    # For LCS, use int dtype
    dtype = "int32" if method == "lcs" else "float64"
    try:
        dp_mat = dp_matrix.load_dp_matrix(dp_bin, dp_txt, shape, dtype=dtype)
    except Exception as exc:
        logger.warning("Failed to load DP matrix for %s: %s", method, exc)
        dp_mat = np.zeros(shape)
    results["dp_shape"] = dp_mat.shape
    results["dp_matrix"] = dp_mat
    # Load path coordinates
    path_coords: List[Tuple[int, int]] = []
    if path_file and path_file.exists():
        path_coords = path_utils.load_path(path_file)
        try:
            path_utils.validate_path_dimensions(path_coords, dp_mat.shape)
        except Exception as exc:
            logger.warning("Invalid path coordinates for %s: %s", method, exc)
    results["path_coords"] = path_coords
    # Path metrics
    path_metrics = path_utils.compute_path_metrics(path_coords)
    results["path_metrics"] = path_metrics
    # Stats metadata (for global/local only)
    if stats_file and stats_file.exists():
        try:
            with stats_file.open("r") as fh:
                stats_meta = json.load(fh)
            results["stats_metadata"] = stats_meta
        except Exception as exc:
            logger.warning("Failed to parse stats file %s: %s", stats_file, exc)
    # Conserved blocks
    blocks_df = block_detection.detect_blocks_to_df(
        seq_a,
        seq_b,
        a_map,
        b_map,
        substitution_matrix or {},
        min_block_length=min_block_length,
        identity_threshold=identity_threshold,
        similarity_threshold=similarity_threshold,
    )
    results["blocks_df"] = blocks_df
    # Residue support (for this method only) – compute for both sequences separately
    method_data_a = {
        "a_map": a_map,
        "b_map": b_map,
        "aligned_a": seq_a,
        "aligned_b": seq_b,
        "dp_matrix": dp_mat,
        "blocks": blocks_df,
    }
    method_data_b = {
        "a_map": b_map,
        "b_map": a_map,
        "aligned_a": seq_b,
        "aligned_b": seq_a,
        "dp_matrix": dp_mat.T,  # transpose for partner orientation
        "blocks": None,  # blocks detection is on seq_a orientation
    }
    # Compute residue support; note that compute_residue_support expects dict of methods -> data; we provide one
    support_df_a = residue_profiles.compute_residue_support(
        len_a, seq_a.replace("-", ""), {method: method_data_a}, window=window
    )
    support_df_b = residue_profiles.compute_residue_support(
        len_b, seq_b.replace("-", ""), {method: method_data_b}, window=window
    )
    # Write outputs
    # Prepare output dir
    outdir.mkdir(parents=True, exist_ok=True)
    # Write alignment summary
    summary_tsv = outdir / f"{out_prefix}_alignment_summary.tsv"
    if overwrite or not summary_tsv.exists():
        pd.DataFrame([stats]).to_csv(summary_tsv, sep="\t", index=False)
    # Write blocks
    blocks_tsv = outdir / f"{out_prefix}_conserved_blocks.tsv"
    if overwrite or not blocks_tsv.exists():
        blocks_df.to_csv(blocks_tsv, sep="\t", index=False)
    # Write path metrics
    path_metrics_tsv = outdir / f"{out_prefix}_path_metrics.tsv"
    if overwrite or not path_metrics_tsv.exists():
        pd.DataFrame([path_metrics]).to_csv(path_metrics_tsv, sep="\t", index=False)
    # Write residue support (two files) – columns names will have method prefix
    support_a_tsv = outdir / f"{out_prefix}_residue_support_{seq_a_id}.tsv"
    if overwrite or not support_a_tsv.exists():
        support_df_a.to_csv(support_a_tsv, sep="\t", index=False)
    support_b_tsv = outdir / f"{out_prefix}_residue_support_{seq_b_id}.tsv"
    if overwrite or not support_b_tsv.exists():
        support_df_b.to_csv(support_b_tsv, sep="\t", index=False)
    # Substitution summary
    subs_df = substitution_analysis.summarise_substitutions(
        seq_a, seq_b, substitution_matrix
    )
    subs_tsv = outdir / f"{out_prefix}_substitution_summary.tsv"
    if overwrite or not subs_tsv.exists():
        subs_df.to_csv(subs_tsv, sep="\t", index=False)
    # Indel regions (gap runs) – use path metrics to derive approximate indel region counts; for now, record path metrics
    # Plot DP heatmap
    heatmap_png = outdir / f"{out_prefix}_dp_heatmap.png"
    if overwrite or not heatmap_png.exists():
        plotting.plot_dp_heatmap(
            dp_mat, heatmap_png, title=f"{method.upper()} DP heatmap", dpi=plot_dpi
        )
    # Plot DP heatmap with path
    heatmap_path_png = outdir / f"{out_prefix}_dp_heatmap_with_path.png"
    if overwrite or not heatmap_path_png.exists():
        plotting.plot_dp_heatmap(
            dp_mat,
            heatmap_path_png,
            path_coords=path_coords,
            title=f"{method.upper()} DP with path",
            dpi=plot_dpi,
        )
    # Plot residue support for both sequences; combine both into one plot for all metrics
    support_plot_a = outdir / f"{out_prefix}_residue_support_{seq_a_id}.png"
    if overwrite or not support_plot_a.exists():
        plotting.plot_residue_support(
            support_df_a,
            [method],
            support_plot_a,
            title=f"{method.upper()} residue support – {seq_a_id}",
            dpi=plot_dpi,
        )
    support_plot_b = outdir / f"{out_prefix}_residue_support_{seq_b_id}.png"
    if overwrite or not support_plot_b.exists():
        plotting.plot_residue_support(
            support_df_b,
            [method],
            support_plot_b,
            title=f"{method.upper()} residue support – {seq_b_id}",
            dpi=plot_dpi,
        )
    return results, support_df_a, support_df_b


def compare_methods(
    support_a: Dict[str, pd.DataFrame],
    seq_a_id: str,
    outdir: Path,
    prefix: str,
    plot_dpi: int,
    overwrite: bool,
) -> Tuple[pd.DataFrame, pd.Series]:
    """Perform comparative analysis across alignment methods for one sequence.

    Parameters
    ----------
    support_a : dict[str, pandas.DataFrame]
        Mapping from method name to residue support DataFrame for sequence A.
    seq_a_id : str
        Identifier of sequence A for file naming.
    outdir : Path
        Directory for output files.
    prefix : str
        Prefix for output files.
    plot_dpi : int
        DPI for plots.
    overwrite : bool
        Overwrite existing files.

    Returns
    -------
    tuple
        (segments_df, category_series) where ``segments_df`` summarises
        contiguous category segments and ``category_series`` lists the
        category assignment for each residue.
    """
    # Build combined residue support DataFrame by merging on residue_index and residue letter
    # Start with base DataFrame from first method
    base_df = None
    for method, df in support_a.items():
        method_prefix = method
        if base_df is None:
            base_df = df[
                ["residue_index", "residue", f"{method_prefix}_participates"]
            ].copy()
        else:
            base_df = base_df.merge(
                df[["residue_index", f"{method_prefix}_participates"]],
                on="residue_index",
                how="outer",
            )
    base_df = base_df.fillna(
        {col: False for col in base_df.columns if col.endswith("_participates")}
    )
    # Rename columns to global/local/lcs participates
    base_df = base_df.rename(
        columns={
            "global_participates": "global_participates",
            "local_participates": "local_participates",
            "lcs_participates": "lcs_participates",
        }
    )
    # Assign categories
    category_series = comparison.assign_participation_categories(base_df)
    # Summarise segments
    segments_df = comparison.summarise_category_segments(category_series)
    # Write category assignments
    cat_tsv = outdir / f"{prefix}_alignment_method_comparison_categories_{seq_a_id}.tsv"
    if overwrite or not cat_tsv.exists():
        pd.DataFrame(
            {"residue_index": category_series.index, "category": category_series.values}
        ).to_csv(cat_tsv, sep="\t", index=False)
    # Write segments summary
    seg_tsv = outdir / f"{prefix}_alignment_method_comparison_{seq_a_id}.tsv"
    if overwrite or not seg_tsv.exists():
        segments_df.to_csv(seg_tsv, sep="\t", index=False)
    # Plot method comparison
    plot_path = outdir / f"{prefix}_alignment_method_comparison_{seq_a_id}.png"
    if overwrite or not plot_path.exists():
        plotting.plot_alignment_method_comparison(
            category_series,
            plot_path,
            title=f"Alignment method comparison – {seq_a_id}",
            dpi=plot_dpi,
        )
    return segments_df, category_series


def main(argv: Optional[List[str]] = None) -> int:
    """Entry point for the alignment analysis CLI.

    Parameters
    ----------
    argv : list[str], optional
        Command line arguments. If None, ``sys.argv[1:]`` is used.

    Returns
    -------
    int
        Exit code: 0 on success, non‑zero on failure.
    """
    parser = argparse.ArgumentParser(description="Protein alignment analysis tool")
    subparsers = parser.add_subparsers(
        dest="command", required=True, help="Subcommands"
    )

    # Subcommands
    for cmd in ["global", "local", "lcs", "full", "compare"]:
        sp = subparsers.add_parser(cmd, help=f"Run {cmd} analysis")
        parse_common_args(sp)

    args = parser.parse_args(argv)
    # Configure logging
    configure_logging(args.log_file, args.quiet)
    # Determine outdir
    outdir: Path = args.outdir or args.results_dir
    # Scan inventory
    try:
        files = file_inventory.scan_results_dir(args.results_dir)
    except Exception as exc:
        logger.error("Failed to scan results directory: %s", exc)
        return 1
    # Load substitution matrix
    substitution_matrix = substitution_analysis.load_substitution_matrix(args.blosum)
    # Selected command
    command = args.command
    # Prepare per‑method results and residue support dictionaries
    all_method_results: Dict[str, Dict[str, Any]] = {}
    all_support_a: Dict[str, pd.DataFrame] = {}
    all_support_b: Dict[str, pd.DataFrame] = {}
    sequence_ids: Optional[Tuple[str, str]] = None
    sequence_lengths: Optional[Tuple[int, int]] = None
    dp_shapes: Dict[str, Tuple[int, int]] = {}
    alignment_stats: Dict[str, Dict[str, float]] = {}
    blocks_top: Dict[str, List[Dict[str, Any]]] = {}
    stats_metadata: Dict[str, Any] = {}
    warnings_list: List[str] = []
    # Determine which methods to run
    methods_to_run: List[str]
    if command in {"global", "local", "lcs"}:
        methods_to_run = [command]
    elif command == "full":
        methods_to_run = [
            m for m in ["global", "local", "lcs"] if getattr(files, f"{m}_alignment")
        ]
    elif command == "compare":
        # For comparison we need at least two methods; run those available
        methods_to_run = [
            m for m in ["global", "local", "lcs"] if getattr(files, f"{m}_alignment")
        ]
        if len(methods_to_run) < 2:
            logger.error("Comparison requires at least two alignment methods available")
            return 1
    else:
        logger.error("Unknown command: %s", command)
        return 1
    # Validate files for each method
    for m in methods_to_run:
        try:
            file_inventory.validate_files(m, files)
        except Exception as exc:
            logger.error("Validation failed for %s: %s", m, exc)
            return 1
    # Analyse each method if needed (skip for compare if not full?)
    for m in methods_to_run:
        # Determine if outputs already exist and skip if not overwrite
        try:
            results, support_df_a, support_df_b = analyse_method(
                m,
                files,
                substitution_matrix,
                outdir,
                args.prefix,
                args.min_block_length,
                args.identity_threshold,
                args.similarity_threshold,
                args.window,
                args.plot_dpi,
                args.overwrite,
            )
        except Exception as exc:
            logger.error("Analysis failed for %s: %s", m, exc)
            return 1
        all_method_results[m] = results
        all_support_a[m] = support_df_a
        all_support_b[m] = support_df_b
        alignment_stats[m] = results.get("alignment_stats", {})
        dp_shapes[m] = results.get("dp_shape", ())
        if "stats_metadata" in results:
            stats_metadata[m] = results["stats_metadata"]
        # Determine sequence ids/lengths once
        if sequence_ids is None:
            sequence_ids = results["sequence_ids"]
            sequence_lengths = (
                results["alignment_stats"]["ungapped_length_a"],
                results["alignment_stats"]["ungapped_length_b"],
            )
        # Determine top blocks: sort by frac_identity descending
        blocks_df = results["blocks_df"]
        if blocks_df is not None and not blocks_df.empty:
            sorted_blocks = (
                blocks_df.sort_values(
                    by=["frac_identity", "frac_similarity"], ascending=False
                )
                .head(3)
                .to_dict(orient="records")
            )
            blocks_top[m] = sorted_blocks
        else:
            blocks_top[m] = []
    # If comparison or full, perform method comparison on sequence A and B separately
    if command in {"compare", "full"}:
        if sequence_ids is None or sequence_lengths is None:
            logger.error("Sequence information missing for comparison")
            return 1
        seq_a_id, seq_b_id = sequence_ids
        # Sequence A comparison
        seg_df_a, cat_series_a = compare_methods(
            all_support_a, seq_a_id, outdir, args.prefix, args.plot_dpi, args.overwrite
        )
        seg_df_b, cat_series_b = compare_methods(
            all_support_b, seq_b_id, outdir, args.prefix, args.plot_dpi, args.overwrite
        )
        # Summarise category counts
        category_counts = {
            seq_a_id: cat_series_a.value_counts().to_dict(),
            seq_b_id: cat_series_b.value_counts().to_dict(),
        }
        # For full: create summary JSON
        if command == "full":
            input_files = {
                "global_alignment": str(files.global_alignment)
                if files.global_alignment
                else None,
                "local_alignment": str(files.local_alignment)
                if files.local_alignment
                else None,
                "lcs_alignment": str(files.lcs_alignment)
                if files.lcs_alignment
                else None,
                "lcs": str(files.lcs) if files.lcs else None,
                "global_dp": str(files.global_dp_bin or files.global_dp_txt)
                if (files.global_dp_bin or files.global_dp_txt)
                else None,
                "local_dp": str(files.local_dp_bin or files.local_dp_txt)
                if (files.local_dp_bin or files.local_dp_txt)
                else None,
                "lcs_dp": str(files.lcs_dp_bin or files.lcs_dp_txt)
                if (files.lcs_dp_bin or files.lcs_dp_txt)
                else None,
            }
            summary_data = summary.build_summary_data(
                input_files=input_files,
                sequence_ids=sequence_ids,
                sequence_lengths=sequence_lengths,
                dp_shapes=dp_shapes,
                stats_metadata=stats_metadata,
                blocks_top=blocks_top,
                alignment_stats=alignment_stats,
                category_counts=category_counts,
                warnings=warnings_list,
                notes=[
                    "Interpretations are based solely on the provided alignment files and DP matrices.",
                    "LCS analysis captures exact matches only and may miss conservative substitutions.",
                ],
            )
            summary_json_path = outdir / f"{args.prefix}_summary.json"
            summary.generate_summary_json(summary_data, summary_json_path)
    logger.info("Analysis completed successfully")
    return 0


if __name__ == "__main__":
    sys.exit(main())
