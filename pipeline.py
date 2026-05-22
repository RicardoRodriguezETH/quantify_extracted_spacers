#!/usr/bin/env python3
"""
quantify_extracted_spacers — pypiper pipeline

Steps:
  1. extract_spacers    — cut target_seq column from counts table
  2. align_spacers      — R: match spacers to genomes, annotate genes
  3. process_alignments — R: merge counts, gene-level aggregation, CPM, plots

Usage:
    python pipeline.py --config config.yaml [pypiper options]
    python pipeline.py --config config.yaml --recover
    python pipeline.py --config config.yaml --new-start
"""

import argparse
import os
import yaml
import pypiper


def parse_args():
    parser = argparse.ArgumentParser(description="quantify_extracted_spacers pipeline")
    parser.add_argument("--config", required=True, help="Path to config.yaml")
    parser = pypiper.add_pypiper_args(parser, groups=["pypiper"])
    return parser.parse_args()


def load_config(config_path):
    """Load YAML config, resolving relative paths against the config file's directory."""
    base = os.path.dirname(os.path.abspath(config_path))
    with open(config_path) as fh:
        cfg = yaml.safe_load(fh)
    # Resolve every value that looks like a path
    for key in ("counts", "output_dir", "genome", "plasmid", "genome_gff", "plasmid_gff"):
        if not os.path.isabs(cfg[key]):
            cfg[key] = os.path.join(base, cfg[key])
    return cfg


def main():
    args   = parse_args()
    cfg    = load_config(args.config)
    scripts_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "scripts")

    out_dir     = cfg["output_dir"]
    spacer_dir  = os.path.join(out_dir, "spacer_level")
    logs_dir    = os.path.join(out_dir, "logs")
    pypiper_dir = os.path.join(out_dir, "logs", "pypiper")
    for d in (out_dir, spacer_dir, logs_dir, pypiper_dir):
        os.makedirs(d, exist_ok=True)

    pm = pypiper.PipelineManager(
        name="quantify_extracted_spacers",
        outfolder=pypiper_dir,
        args=args,
    )

    # ------------------------------------------------------------------
    # Step 1 — extract_spacers
    # ------------------------------------------------------------------
    pm.timestamp("### Step 1: Extract spacers")

    all_spacers = os.path.join(spacer_dir, "all_spacers.tsv")
    pm.run(
        f"cut -f1 {cfg['counts']} > {all_spacers}",
        target=all_spacers,
        lock_name="extract_spacers",
    )

    # ------------------------------------------------------------------
    # Step 2 — align_spacers
    # ------------------------------------------------------------------
    pm.timestamp("### Step 2: Align spacers")

    spacers_aligned = os.path.join(spacer_dir, "spacers_aligned.tsv")

    pm.run(
        f"""Rscript {scripts_dir}/align_spacers.R \\
    {all_spacers} \\
    {cfg['genome']} \\
    {cfg['plasmid']} \\
    {cfg['genome_gff']} \\
    {cfg['plasmid_gff']} \\
    {cfg['cores']} \\
    {cfg['chunk_size']} \\
    {out_dir}""",
        target=spacers_aligned,
        lock_name="align_spacers",
    )

    # ------------------------------------------------------------------
    # Step 3 — process_alignments
    # ------------------------------------------------------------------
    pm.timestamp("### Step 3: Process alignments")

    pm.run(
        f"""Rscript {scripts_dir}/process_alignments.R \\
    {spacers_aligned} \\
    {cfg['counts']} \\
    {out_dir} \\
    "" \\
    {cfg.get('cpm_max_thresh', 7)} \\
    {cfg.get('cpm_range_thresh', 0.5)}""",
        target=os.path.join(out_dir, "gene_level", "gene_counts.tsv"),
        lock_name="process_alignments",
    )

    pm.stop_pipeline()


if __name__ == "__main__":
    main()
