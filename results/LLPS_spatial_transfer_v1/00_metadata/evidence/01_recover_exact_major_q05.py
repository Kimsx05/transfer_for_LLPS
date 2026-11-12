#!/usr/bin/env python3
"""Serial, resumable 1000-draw recovery of exact major q05. Never calls train()."""
from __future__ import annotations

import argparse
import gc
import json
import os
import time
import traceback
from pathlib import Path

for key in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(key, "6")

import anndata as ad
import numpy as np
import pandas as pd
import torch
from cell2location.models import Cell2location

ROOT = Path("/home/data/t070721/codex_workspace/Bladder Metabolism")
OUT = ROOT / "neighborhood_v1/00_inputs_and_audit"
MANIFEST_PATH = OUT / "qc/posterior_sampling_manifest.csv"
MAPPING_PATH = OUT / "subtype_to_major_mapping.csv"
MAJOR_DIR = OUT / "major_q05"
LOG_DIR = OUT / "logs"
N_DRAWS = 1000
SEED = 20260922


def atomic_csv(df: pd.DataFrame, path: Path) -> None:
    tmp = path.with_suffix(path.suffix + ".tmp")
    df.to_csv(tmp, index=False)
    os.replace(tmp, path)


def set_manifest(section: str, **kwargs) -> None:
    m = pd.read_csv(MANIFEST_PATH)
    for col in ["status", "runtime", "factor_match", "output_path", "notes", "compute_device", "gpu_name"]:
        if col in m.columns:
            m[col] = m[col].astype("object")
    ix = m.index[m["section"].astype(str).eq(section)]
    if len(ix) != 1:
        raise RuntimeError(f"Manifest section lookup failed: {section}")
    for k, v in kwargs.items():
        if k not in m.columns:
            m[k] = ""
        m.loc[ix, k] = v
    atomic_csv(m, MANIFEST_PATH)


def validate_existing(path: Path, n_spots: int, majors: list[str]) -> bool:
    if not path.exists():
        return False
    try:
        x = pd.read_csv(path)
        expected = ["barcode", "section", "dataset"] + [f"q05__{m}" for m in majors] + [f"mean__{m}" for m in majors]
        return len(x) == n_spots and all(c in x.columns for c in expected) and not x[expected].isna().any().any()
    except Exception:
        return False


def run_section(row: pd.Series, mapping: pd.DataFrame, compute_device: str, batch_size: int) -> bool:
    section = str(row["section"]); dataset = str(row["dataset"])
    section_dir = Path(row["section_output_path"])
    adata_path = section_dir / "spatial_posterior.h5ad"
    model_path = section_dir / "model"
    output_path = MAJOR_DIR / f"{section}__exact_major_1000draw.csv.gz"
    majors = sorted(mapping["major"].unique().tolist())
    n_spots_expected = int(row["n_spots"])
    if validate_existing(output_path, n_spots_expected, majors):
        set_manifest(section, status="PASS", n_draws=N_DRAWS, output_path=str(output_path), notes="Validated existing completed output; skipped recomputation")
        print(f"SKIP_VALID {section}", flush=True)
        return True
    log_path = LOG_DIR / f"posterior_{section}.log"
    lines = []
    def log(msg):
        s = str(msg); lines.append(s); print(s, flush=True); log_path.write_text("\n".join(lines)+"\n", encoding="utf-8")
    started = time.time()
    use_cuda = compute_device == "CUDA"
    accelerator = "gpu" if use_cuda else "cpu"
    gpu_name = torch.cuda.get_device_name(0) if use_cuda and torch.cuda.is_available() else ""
    effective_batch = min(batch_size, n_spots_expected)
    set_manifest(section, status="RUNNING", n_draws=N_DRAWS, compute_device=compute_device,
                 gpu_name=gpu_name, batch_size=effective_batch,
                 notes=f"{compute_device} serial sampler started")
    try:
        torch.manual_seed(SEED)
        np.random.seed(SEED)
        torch.set_num_threads(6)
        try: torch.set_num_interop_threads(1)
        except RuntimeError: pass
        if use_cuda and not torch.cuda.is_available():
            raise RuntimeError("CUDA requested but torch.cuda.is_available() is False")
        log(f"START section={section} dataset={dataset} n_draws={N_DRAWS} device={compute_device} gpu={gpu_name!r} batch_size={effective_batch}")
        a = ad.read_h5ad(adata_path)
        model = Cell2location.load(model_path, adata=a, accelerator=accelerator, device=1)
        factors = [str(x) for x in model.factor_names_]
        expected_factors = mapping["subtype"].astype(str).tolist()
        config_factors = [str(x) for x in json.loads((section_dir / "model_config.json").read_text())["labels"]]
        factor_match = factors == config_factors and set(factors) == set(expected_factors) and len(factors) == 30
        if not factor_match:
            raise RuntimeError(f"Factor mismatch model={factors}; config={config_factors}; mapping={expected_factors}")
        log(f"MODEL_LOADED spots={a.n_obs} factors={len(factors)} factor_match=PASS")
        result = model.sample_posterior(
            adata=a, num_samples=N_DRAWS, return_sites=["w_sf"], accelerator=accelerator, device=1,
            batch_size=min(effective_batch, a.n_obs), return_samples=True,
        )
        draws = np.asarray(result["posterior_samples"]["w_sf"])
        expected_shape = (N_DRAWS, a.n_obs, len(factors))
        if draws.shape != expected_shape:
            raise RuntimeError(f"w_sf shape {draws.shape}, expected {expected_shape}")
        if not np.isfinite(draws).all() or (draws < 0).any():
            raise RuntimeError("Non-finite or negative posterior draws")
        out = pd.DataFrame({"barcode": a.obs_names.astype(str), "section": section, "dataset": dataset})
        same_draw_errors = []
        stored_mean = pd.DataFrame(a.obsm["means_cell_abundance_w_sf"], index=a.obs_names)
        stored_cols = [str(c).removeprefix("meanscell_abundance_w_sf_").replace("__SLASH__", "/") for c in stored_mean.columns]
        stored_mean.columns = stored_cols
        stored_comparison = []
        factor_to_ix = {f:i for i,f in enumerate(factors)}
        for major in majors:
            subs = mapping.loc[mapping["major"].eq(major), "subtype"].astype(str).tolist()
            ix = [factor_to_ix[s] for s in subs]
            # Promote before reduction so the two mathematically equivalent mean/sum
            # orders are tested without float32 accumulation-order artifacts.
            selected = draws[:, :, ix].astype(np.float64, copy=False)
            major_draws = selected.sum(axis=2, dtype=np.float64)
            q05 = np.quantile(major_draws, .05, axis=0)
            draw_mean = major_draws.mean(axis=0)
            algebraic_mean = selected.mean(axis=0, dtype=np.float64).sum(axis=1, dtype=np.float64)
            same_draw_errors.append(float(np.max(np.abs(draw_mean - algebraic_mean))))
            stored_sum = stored_mean[subs].sum(axis=1).to_numpy()
            stored_comparison.append(float(np.max(np.abs(draw_mean - stored_sum))))
            out[f"q05__{major}"] = q05
            out[f"mean__{major}"] = draw_mean
            del major_draws
        same_draw_error = max(same_draw_errors)
        if not np.isclose(same_draw_error, 0.0, rtol=1e-5, atol=1e-6):
            raise RuntimeError(f"Same-draw mean linearity failed: {same_draw_error}")
        out.to_csv(output_path, index=False, compression="gzip")
        runtime = time.time() - started
        stored_max = max(stored_comparison)
        notes = ("Exact same-draw linearity passed. Comparison with pre-existing posterior mean is an independent Monte Carlo-export diagnostic, "
                 f"not an algebraic identity; max_abs_difference={stored_max:.8g}.")
        set_manifest(section, status="PASS", runtime=f"{runtime:.1f}s", factor_match="PASS_EXACT_SET_AND_MODEL_ORDER",
                     mean_linearity_max_abs_error=same_draw_error, output_path=str(output_path),
                     compute_device=compute_device, gpu_name=gpu_name, batch_size=effective_batch, notes=notes)
        log(f"PASS runtime={runtime:.1f}s same_draw_mean_error={same_draw_error:.8g} stored_mean_mc_max_abs_diff={stored_max:.8g} output={output_path}")
        return True
    except Exception as exc:
        runtime = time.time() - started
        tb = traceback.format_exc(); log(f"FAIL {exc!r}\n{tb}")
        set_manifest(section, status="FAIL", runtime=f"{runtime:.1f}s", factor_match="FAIL_OR_NOT_REACHED",
                     output_path="", compute_device=compute_device, gpu_name=gpu_name,
                     batch_size=effective_batch, notes=repr(exc))
        return False
    finally:
        for name in ["draws", "result", "model", "a", "out"]:
            if name in locals(): del locals()[name]
        gc.collect()


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--section", help="Run one section only")
    ap.add_argument("--all", action="store_true", help="Run all pending/failed sections serially")
    ap.add_argument("--device", choices=["CPU", "CUDA"], default="CPU")
    ap.add_argument("--batch-size", type=int, default=512)
    ap.add_argument("--pending-only", action="store_true", help="Only run PENDING/unfinished rows; PASS rows are always preserved")
    args = ap.parse_args()
    if not args.section and not args.all:
        ap.error("Specify --section or --all")
    MAJOR_DIR.mkdir(parents=True, exist_ok=True); LOG_DIR.mkdir(parents=True, exist_ok=True)
    mapping = pd.read_csv(MAPPING_PATH)
    if len(mapping) != 30 or mapping["subtype"].nunique() != 30 or mapping["major"].isna().any():
        raise RuntimeError("Frozen mapping is not valid 30/30")
    manifest = pd.read_csv(MANIFEST_PATH)
    if args.pending_only:
        manifest = manifest.loc[~manifest["status"].astype(str).eq("PASS")]
    if args.section:
        manifest = manifest.loc[manifest["section"].astype(str).eq(args.section)]
        if len(manifest) != 1: raise RuntimeError(f"Unknown section {args.section}")
    for _, row in manifest.iterrows():
        run_section(row, mapping, args.device, args.batch_size)


if __name__ == "__main__":
    main()
