import os
import logging
import numpy as np
import pandas as pd
from typing import Dict, Optional, Iterable, Any, Tuple, List
import pyranges as pr

logging.basicConfig(level=logging.INFO)


class INTEGRATE:
    """
    Simplified INTEGRATE class.

    Key changes vs. v1
    ------------------
    * gene column  → always phenotype_id (no caqtl/sqtl gene mapping at all)
    * integrate key → (locus_range, lead_snp) instead of lead_snp alone
    * removed: caqtl_eqtl_coloc_gene_mapper, nearest-gene lookup, _ensure_gene_column
    * removed: gene_mapping output table (no mapping = no table)
    * kept:    HESS attachment, best-per-(phenotype, locus_range, lead_snp) pick
    """

    GIM_METRICS_DEFAULT: Dict[str, List[str]] = {
        "coloc":     ["PP.H3.abf", "PP.H4.abf"],
        "smr":       ["p_SMR", "p_HEIDI"],
        "fastenloc": ["GRCP"],
        "ecaviar":   ["clpp", "gwas_pip", "qtl_pip"],
    }

    # phenotype_id column name per QTL study type
    PHENO_COL_DEFAULT: Dict[str, str] = {
        "eqtl":  "phenotype_id",
        "pqtl":  "phenotype_id",
        "sqtl":  "phenotype_id",
        "caqtl": "phenotype_id",
    }

    # ------------------------------------------------------------------ init --

    def __init__(
        self,
        genecode_file: Optional[str] = None,   # kept for HESS region annotation if needed
        h3h4_threshold: float = 0.5,           # no longer used for mapping, kept for compat
        lead_col: str = "lead_snp",
        locus_col: str = "locus_range",        # NEW: column carrying the locus range
        gim_metrics: Optional[Dict[str, Iterable[str]]] = None,
        pheno_col_map: Optional[Dict[str, str]] = None,
    ):
        self.genecode_file   = genecode_file
        self.h3h4_threshold  = h3h4_threshold
        self.lead_col        = lead_col
        self.locus_col       = locus_col       # NEW
        self.GIM_metric      = gim_metrics or self.GIM_METRICS_DEFAULT
        self.phenotype_col_studytype = pheno_col_map or self.PHENO_COL_DEFAULT

        self.hess_by_gwas: Dict[str, pd.DataFrame] = {}

    # --------------------------------------------------------------- run_all --

    def run_all(
        self,
        integrate_result_dict: Dict[str, dict],
        output_dir: str,
    ) -> pd.DataFrame:
        """
        Process all configs and write one output file:
            integrate_results.tsv

        Columns (long format):
            phenotype_id | locus_range | lead_snp | study_type | context
            | gwas | population | gim | metric | value
            + HESS columns (h2, hess_region, CHR, START, END, P)

        Returns the merged long DataFrame.
        """
        os.makedirs(output_dir, exist_ok=True)

        # 1) collect HESS (once per GWAS, from gwas–eqtl cfgs)
        self._collect_hess_by_gwas(integrate_result_dict)

        # 2) iterate cfgs
        long_rows: List[dict] = []
        for cfg_key, cfg in integrate_result_dict.items():
            t1 = cfg.get("trait1", "")
            t2 = cfg.get("trait2", "")
            # skip pure QTL–QTL pairs (eqtl vs caqtl etc.)
            if t1 not in ("gwas",) and t2 not in ("gwas",):
                logging.info(f"[{cfg_key}] skipping non-GWAS pair ({t1} vs {t2})")
                continue
            logging.info(f"[{cfg_key}] processing {t1}:{cfg.get('trait1_name')} "
                         f"vs {t2}:{cfg.get('trait2_name')}")
            rows = self._process_single_cfg(cfg_key, cfg)
            long_rows.extend(rows)

        # 3) assemble long table
        merged_long = pd.DataFrame(long_rows, columns=[
            "phenotype_id", "locus_range", "lead_snp",
            "study_type", "context", "gwas", "population",
            "gim", "metric", "value",
        ])

        # 4) attach HESS
        merged_long_hess = self._attach_hess_by_gwas(merged_long)

        # 5) write
        out_path = os.path.join(output_dir, "integrate_results.tsv")
        merged_long_hess.to_csv(out_path, sep="\t", index=False)
        logging.info(f"[INTEGRATE] Written {len(merged_long_hess)} rows → {out_path}")

        return merged_long_hess

    # ---------------------------------------------------------- HESS helpers --

    def _collect_hess_by_gwas(self, integrate_result_dict: Dict[str, dict]) -> None:
        """Load HESS files once per GWAS name."""
        hessstore: Dict[str, pd.DataFrame] = {}
        for cfg_key, cfg in integrate_result_dict.items():
            if cfg.get("trait1") == "gwas" and cfg.get("trait2") == "eqtl":
                gwas_name = cfg.get("trait1_name")
                hess_path = (cfg.get("GIM_results") or {}).get("hess")
                if gwas_name and gwas_name not in hessstore and self._path_exists(hess_path):
                    df = self._read_any(hess_path)
                    if self._is_df(df):
                        hessstore[gwas_name] = df.copy()
                        logging.info(f"[HESS] loaded for GWAS '{gwas_name}': {len(df)} rows")
        self.hess_by_gwas = hessstore

    def _attach_hess_by_gwas(self, merged_long: pd.DataFrame) -> pd.DataFrame:
        if not self._is_df(merged_long):
            return merged_long

        out_all: List[pd.DataFrame] = []
        for gwas_name, sub in merged_long.groupby("gwas", dropna=False):
            hess_df = self.hess_by_gwas.get(gwas_name)
            if not self._is_df(hess_df):
                temp = sub.copy()
                for c in ["h2", "hess_region", "CHR", "START", "END", "P"]:
                    temp[c] = np.nan
                out_all.append(temp)
                continue

            hess = hess_df.copy()
            req = ["CHR", "START", "END", "NAME", "P", "h2"]
            if any(c not in hess.columns for c in req):
                logging.warning(f"[HESS] missing required columns for '{gwas_name}'")
                temp = sub.copy()
                for c in ["h2", "hess_region", "CHR", "START", "END", "P"]:
                    temp[c] = np.nan
                out_all.append(temp)
                continue

            hess["lead_snp_h"] = hess["NAME"].astype(str).str.strip()
            hess["hess_region"] = (
                hess["CHR"].astype(str).str.strip() + ":" +
                hess["START"].astype(str).str.strip() + "-" +
                hess["END"].astype(str).str.strip()
            )
            hess["_h2_num"] = pd.to_numeric(hess["h2"], errors="coerce").fillna(float("-inf"))
            hess = hess.sort_values("_h2_num", ascending=False).drop_duplicates("lead_snp_h")

            temp = sub.copy()
            temp["lead_snp_h"] = temp["lead_snp"].astype(str).str.strip()
            temp = temp.merge(
                hess[["lead_snp_h", "h2", "hess_region", "CHR", "START", "END", "P"]],
                on="lead_snp_h", how="left"
            ).drop(columns=["lead_snp_h"])
            temp["h2"] = pd.to_numeric(temp["h2"], errors="coerce").fillna(0.0)
            out_all.append(temp)

        return pd.concat(out_all, axis=0, ignore_index=True)

    # ------------------------------------------------------ per-cfg processing --

    def _process_single_cfg(self, cfg_key: str, cfg: dict) -> List[dict]:
        """
        Returns long-format rows for one cfg.
        Each row is keyed by (phenotype_id, locus_range, lead_snp).
        """
        t1, t1name = cfg.get("trait1"), cfg.get("trait1_name")
        t2, t2name = cfg.get("trait2"), cfg.get("trait2_name")
        pop        = cfg.get("population", "NA")
        results    = dict(cfg.get("GIM_results") or {})

        # determine which side is the QTL
        qtl_type = t2 if t1 == "gwas" else t1
        gwas_label = t1name if t1 == "gwas" else t2name

        if qtl_type not in self.phenotype_col_studytype:
            logging.warning(f"[{cfg_key}] unsupported QTL type '{qtl_type}'")
            return []

        pheno_col = self.phenotype_col_studytype[qtl_type]   # always "phenotype_id"

        # HESS handled globally
        results.pop("hess", None)
        # results.pop("susie", None)   # susie skipped for now

        long_rows: List[dict] = []

        for gim, path in results.items():
            if not self._path_exists(path):
                logging.warning(f"[{cfg_key}/{gim}] file not found: {path}")
                continue
            df = self._read_any(path)
            if not self._is_df(df):
                logging.warning(f"[{cfg_key}/{gim}] empty/invalid file: {path}")
                continue

            # --- verify required columns ---
            needed = [pheno_col, self.lead_col, self.locus_col]
            missing = [c for c in needed if c not in df.columns]
            if missing:
                logging.warning(f"[{cfg_key}/{gim}] missing columns {missing} — skipping")
                continue

            # --- clean & deduplicate ---
            df = df.copy()
            df[pheno_col]      = df[pheno_col].astype(str).str.strip()
            df[self.lead_col]  = df[self.lead_col].astype(str).str.strip()
            df[self.locus_col] = df[self.locus_col].astype(str).str.strip()

            # unique per (phenotype, locus_range, lead_snp)
            df = df.drop_duplicates([pheno_col, self.locus_col, self.lead_col])

            # --- pick best record per (phenotype, locus_range, lead_snp) ---
            df_best = self._pick_best_per_pheno_locus_lead(
                df, gim, pheno_col, self.locus_col, self.lead_col
            )

            metrics = [m for m in self.GIM_metric.get(gim, []) if m in df_best.columns]

            for _, r in df_best.iterrows():
                pheno  = r[pheno_col]
                locus  = r[self.locus_col]
                lead   = r[self.lead_col]

                if metrics:
                    for m in metrics:
                        long_rows.append({
                            "phenotype_id": pheno,
                            "locus_range":  locus,
                            "lead_snp":     lead,
                            "study_type":   qtl_type,
                            "context":      t2name if t1 == "gwas" else t1name,
                            "gwas":         gwas_label,
                            "population":   pop,
                            "gim":          gim,
                            "metric":       m,
                            "value":        r.get(m),
                        })
                else:
                    long_rows.append({
                        "phenotype_id": pheno,
                        "locus_range":  locus,
                        "lead_snp":     lead,
                        "study_type":   qtl_type,
                        "context":      t2name if t1 == "gwas" else t1name,
                        "gwas":         gwas_label,
                        "population":   pop,
                        "gim":          gim,
                        "metric":       None,
                        "value":        None,
                    })

        return long_rows

    # ------------------------------------------------- best-pick per key tuple --

    @staticmethod
    def _pick_best_per_pheno_locus_lead(
        df: pd.DataFrame,
        gim: str,
        pheno_col: str,
        locus_col: str,
        lead_col: str,
    ) -> pd.DataFrame:
        """
        For each (phenotype_id, locus_range, lead_snp) group keep the single
        most significant / highest-confidence row according to GIM type.
        """
        group_keys = [pheno_col, locus_col, lead_col]
        df = df.copy()

        score_col, ascending = None, False

        if gim == "smr":
            if "p_SMR" in df.columns:
                score_col, ascending = "p_SMR", True       # smaller p = better
        elif gim == "coloc":
            if "overall_H4" in df.columns:
                score_col, ascending = "overall_H4", False
            elif "PP.H3.abf" in df.columns:
                score_col, ascending = "PP.H3.abf", False
        elif gim == "fastenloc":
            if "GRCP" in df.columns:
                score_col, ascending = "GRCP", False
        elif gim == "ecaviar":
            if "clpp" in df.columns:
                score_col, ascending = "clpp", False
            elif "gwas_pip" in df.columns:
                score_col, ascending = "gwas_pip", False

        if score_col is None or score_col not in df.columns:
            return df.drop_duplicates(group_keys, keep="first")

        df["_score"] = pd.to_numeric(df[score_col], errors="coerce")
        if ascending:
            df["_score"] = df["_score"].fillna(float("inf"))
            idx = df.groupby(group_keys)["_score"].idxmin()
        else:
            df["_score"] = df["_score"].fillna(float("-inf"))
            idx = df.groupby(group_keys)["_score"].idxmax()

        return df.loc[idx].drop(columns=["_score"])

    # ------------------------------------------------------------ IO helpers --

    @staticmethod
    def _path_exists(path: Any) -> bool:
        try:
            return isinstance(path, str) and os.path.exists(path)
        except Exception:
            return False

    @staticmethod
    def _read_any(path: str) -> Optional[pd.DataFrame]:
        if path is None:
            return None
        p = str(path).lower()
        try:
            if p.endswith(".parquet"):
                return pd.read_parquet(path)
            if p.endswith((".tsv", ".txt", ".tsv.gz", ".txt.gz")):
                return pd.read_csv(path, sep="\t")
            return pd.read_csv(path)
        except Exception as e:
            logging.warning(f"[READ] failed: {path} — {e}")
            return None

    @staticmethod
    def _is_df(df: Any) -> bool:
        return isinstance(df, pd.DataFrame) and not df.empty