"""Command line entry point.

    python -m cish2 run      configs/baseline.yaml [-o key=value ...]
    python -m cish2 plan     configs/baseline.yaml
    python -m cish2 collect  results/<name>
    python -m cish2 analyze  results/<name>
"""

from __future__ import annotations

import argparse
import os
import sys

from .campaign import run_campaign
from .config import CampaignConfig
from .results import ResultStore


def _load(args) -> CampaignConfig:
    cfg = CampaignConfig.from_yaml(args.config, args.override)
    if args.name:
        cfg.name = args.name
    if args.results_root:
        cfg.results_root = args.results_root
    cfg.validate()
    return cfg


def cmd_plan(args) -> int:
    cfg = _load(args)
    print(f"campaign      : {cfg.name}")
    print(f"fingerprint   : {cfg.fingerprint}")
    print(f"run dir       : {cfg.run_dir}")
    print(f"LD arch       : {cfg.region.ld_arch}   LD for SuSiE: {cfg.method.ld_for_susie}")
    print(f"grid          : n={cfg.grid.n_individuals} ncausal={cfg.grid.ncausal}")
    print(f"                h2={cfg.grid.h2} dist={cfg.grid.effect_dist} maf={cfg.grid.maf_bin}")
    print(f"cells/region  : {cfg.grid.size}")
    print(f"regions       : {cfg.n_regions}")
    print(f"total cells   : {cfg.total_cells}")
    print(f"workers       : {cfg.compute.workers} x {cfg.compute.threads_per_worker} thread(s)")
    done = ResultStore(cfg.run_dir).completed_slots() if os.path.isdir(cfg.run_dir) else set()
    print(f"already done  : {len(done)} slots")
    return 0


def cmd_run(args) -> int:
    cfg = _load(args)
    run_campaign(cfg, resume=not args.no_resume)
    return 0


def cmd_gtex(args) -> int:
    from .gtex_campaign import load_config, run_gtex

    cfg = load_config(args.config, args.override)
    if args.name:
        cfg["name"] = args.name
    run_gtex(cfg, resume=not args.no_resume)
    return 0


def cmd_collect(args) -> int:
    store = ResultStore(args.run_dir)
    out = store.consolidate()
    df = store.load()
    print(f"{len(df):,} rows -> {out}")
    return 0


def cmd_analyze(args) -> int:
    from analysis.overview import is_simulation_run

    if not is_simulation_run(args.run_dir):
        # a real-data run has no true h2; it gets the prediction-test report
        from analysis.gtex_report import build_gtex_report

        print(build_gtex_report(args.run_dir))
        return 0

    from analysis.report import build_report

    build_report(args.run_dir, compare=args.compare)
    return 0


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(prog="cish2", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)

    def add_config_args(sp):
        sp.add_argument("config", help="campaign YAML")
        sp.add_argument("-o", "--override", action="append", default=[],
                        metavar="KEY=VALUE", help="dotted config override, repeatable")
        sp.add_argument("--name", default=None, help="override campaign name")
        sp.add_argument("--results-root", default=None, help="override results root")

    sp = sub.add_parser("plan", help="print what a config would run")
    add_config_args(sp)
    sp.set_defaults(func=cmd_plan)

    sp = sub.add_parser("run", help="run (or resume) a campaign")
    add_config_args(sp)
    sp.add_argument("--no-resume", action="store_true", help="ignore existing manifest")
    sp.set_defaults(func=cmd_run)

    sp = sub.add_parser("gtex", help="run the GTEx real-data pass")
    sp.add_argument("config", help="GTEx YAML")
    sp.add_argument("-o", "--override", action="append", default=[], metavar="KEY=VALUE")
    sp.add_argument("--name", default=None)
    sp.add_argument("--no-resume", action="store_true")
    sp.set_defaults(func=cmd_gtex)

    sp = sub.add_parser("collect", help="merge shards into cells.tsv.gz")
    sp.add_argument("run_dir")
    sp.set_defaults(func=cmd_collect)

    sp = sub.add_parser("analyze", help="build tables and figures for a run")
    sp.add_argument("run_dir")
    sp.add_argument("--compare", action="append", default=[],
                    help="additional run dir to overlay, repeatable")
    sp.set_defaults(func=cmd_analyze)

    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
