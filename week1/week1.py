"""Week 1 term project — LIN28A CLIP / ribosome footprinting reanalysis.

원본 노트북 ``week1.ipynb`` (Google Colab 용) 의 데이터 수집 → featureCounts
→ pandas 분석 흐름을 그대로 따라가는 단일 실행 스크립트입니다.
직접 작성하는 부분이 비어 있던 두 코드 셀

* 셀 33: "논문의 그림 Figure 4D 처럼 한 번 만들어 봅시다"
* 셀 37: "localization 데이터와 위에서 만든 scatter 를 결합해서
        논문 그림(Figure 5B, S6A) 과 비슷하게 만들어 봅시다"

을 채워서 ``result1.png``, ``result2.png`` 로 저장합니다.

데이터셋
--------
``https://hyeshik.qbio.io/binfo/binfo1-datapack1.tar`` 에 들어있는 6 개의
BAM (CLIP, RNA-seq, ribosome profiling) 과 GENCODE M27 GTF 입니다.
이는 Cho et al., *Cell* (2012) "LIN28A is a suppressor of ER-associated
translation in embryonic stem cells" 의 시퀀싱 데이터 입니다.

요구 환경
--------
* ``featureCounts`` (subread) 가 PATH 에 있어야 함
* Python: numpy, pandas, matplotlib, scipy

Usage
-----
    python week1.py
"""

from __future__ import annotations

import gzip
import shutil
import ssl
import subprocess
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy.stats import pearsonr


# --- 경로 설정 (원본 Colab 노트북의 /content/drive/MyDrive 와 대응) ----------
HERE = Path(__file__).resolve().parent
ASSIGNMENT_ROOT = HERE.parent
WORKDIR = ASSIGNMENT_ROOT / "work"               # = binfo1-work
DATAPACK_DIR = WORKDIR / "binfo1-datapack1"      # = binfo1-datapack1
RESULTS_DIR = HERE / "results"

# --- 원본 노트북에 명시된 다운로드 URL --------------------------------------
DATAPACK_URL = "https://hyeshik.qbio.io/binfo/binfo1-datapack1.tar"
GENCODE_URL = (
    "http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/"
    "release_M27/gencode.vM27.annotation.gtf.gz"
)
LOCALIZATION_URL = "https://hyeshik.qbio.io/binfo/mouselocalization-20210507.txt"

EXPECTED_BAMS = [
    "CLIP-35L33G.bam",
    "RNA-control.bam",
    "RNA-siLin28a.bam",
    "RNA-siLuc.bam",
    "RPF-siLin28a.bam",
    "RPF-siLuc.bam",
]
SAMPLE_COLS = list(EXPECTED_BAMS)  # featureCounts 출력의 컬럼명과 일치

# 논문 Figure 5B / S6A 색상 (nucleus = 파랑, integral membrane = 빨강,
# cytoplasm = 초록)
LOC_COLORS = {
    "nucleus":           "#1f6fd0",
    "integral membrane": "#d9322f",
    "cytoplasm":         "#2ca02c",
}


# ---------------------------------------------------------------------------
# 1. 데이터 수집 (원본 노트북 셀 9, 11)
# ---------------------------------------------------------------------------

def _download(url: str, dest: Path, timeout: int = 60) -> None:
    """원본 노트북의 ``wget --no-check-certificate`` 와 같은 동작.

    실패하면 예외를 그대로 올린다.
    """
    print(f"[wget] {url}")
    ctx = ssl._create_unverified_context()
    with urllib.request.urlopen(url, context=ctx, timeout=timeout) as r, \
            open(dest, "wb") as f:
        shutil.copyfileobj(r, f)


def setup_datapack() -> None:
    """원본 노트북 셀 9 / 11 의 동작을 그대로 수행.

    * BAM tar 다운로드 → tar 풀기
    * GENCODE M27 GTF 다운로드 → gunzip
    이미 받아놓은 파일이 있으면 skip.
    """
    DATAPACK_DIR.mkdir(parents=True, exist_ok=True)

    # ---- BAM/BAI ----
    missing = [b for b in EXPECTED_BAMS if not (DATAPACK_DIR / b).exists()]
    if missing:
        tar_path = WORKDIR / "binfo1-datapack1.tar"
        if not tar_path.exists():
            _download(DATAPACK_URL, tar_path)
        print(f"[tar] extracting {tar_path}")
        subprocess.run(["tar", "-C", str(WORKDIR), "-xf", str(tar_path)],
                       check=True)

    # ---- GENCODE GTF ----
    gtf = DATAPACK_DIR / "gencode.gtf"
    if not gtf.exists():
        gtf_gz = DATAPACK_DIR / "gencode.gtf.gz"
        if not gtf_gz.exists():
            _download(GENCODE_URL, gtf_gz)
        print(f"[gunzip] {gtf_gz.name}")
        with gzip.open(gtf_gz, "rb") as fin, open(gtf, "wb") as fout:
            shutil.copyfileobj(fin, fout)
        gtf_gz.unlink()


# ---------------------------------------------------------------------------
# 2. 작업 디렉토리 준비 (원본 노트북 셀 22, 24)
# ---------------------------------------------------------------------------

def setup_workdir() -> None:
    """원본 노트북: ``cp -r binfo1-datapack1/* binfo1-work/`` 와 대응.

    공간 절약을 위해 cp 대신 symlink 를 건다 (featureCounts 결과는 동일).
    """
    WORKDIR.mkdir(parents=True, exist_ok=True)
    for entry in DATAPACK_DIR.iterdir():
        if entry.suffix in {".bam", ".bai"} or entry.name == "gencode.gtf":
            link = WORKDIR / entry.name
            if not link.exists():
                link.symlink_to(entry.resolve())


# ---------------------------------------------------------------------------
# 3. featureCounts (원본 노트북 셀 26)
# ---------------------------------------------------------------------------

def run_feature_counts() -> Path:
    out = WORKDIR / "read-counts.txt"
    if out.exists():
        print(f"[skip] {out.name} already exists")
        return out

    feature_counts = shutil.which("featureCounts")
    if feature_counts is None:
        raise RuntimeError(
            "featureCounts not found in PATH. "
            "Please install subread (e.g. `conda install -c bioconda subread`)."
        )

    cmd = [feature_counts, "-T", "8",
           "-a", "gencode.gtf",
           "-o", "read-counts.txt",
           *EXPECTED_BAMS]
    print(f"[run] {' '.join(cmd)}")
    subprocess.run(cmd, cwd=WORKDIR, check=True)
    return out


# ---------------------------------------------------------------------------
# 4. 데이터 가공 (원본 노트북 셀 30 의 정규화·필터링 강화)
# ---------------------------------------------------------------------------

def compute_metrics(cnts: pd.DataFrame, min_reads: int) -> pd.DataFrame:
    """원본 셀 30 에서 raw count 로 바로 비율을 잡던 부분을 보강.

    * library size 정규화 (CPM): 라이브러리 크기 차이 보정
    * x 축: log2( CLIP-35L33G CPM / RNA-control CPM ) — CLIP enrichment
    * y 축: log2( (RPF-siLin28a / RNA-siLin28a) / (RPF-siLuc / RNA-siLuc) )
            — Lin28a knockdown 시 ribosome density 변화
    * 모든 sample 에서 ``min_reads`` 이상인 gene 만 사용 (저카운트 비율 노이즈 제거)
    """
    size_factors = cnts[SAMPLE_COLS].sum()
    cpm = cnts[SAMPLE_COLS].div(size_factors, axis=1) * 1e6

    keep = (cnts[SAMPLE_COLS] >= min_reads).all(axis=1)
    c = cpm[keep]

    clip_enrich = np.log2(c["CLIP-35L33G.bam"] / c["RNA-control.bam"])
    te_kd   = c["RPF-siLin28a.bam"] / c["RNA-siLin28a.bam"]
    te_ctrl = c["RPF-siLuc.bam"]    / c["RNA-siLuc.bam"]
    rden_change = np.log2(te_kd / te_ctrl)

    df = pd.DataFrame({"clip_enrich": clip_enrich, "rden_change": rden_change})
    return df.replace([np.inf, -np.inf], np.nan).dropna()


def annotate_localization(df: pd.DataFrame) -> pd.DataFrame:
    """원본 노트북 셀 35 의 localization 파일을 받아서 gene 별 매핑."""
    ssl._create_default_https_context = ssl._create_unverified_context
    mouselocal = pd.read_csv(LOCALIZATION_URL, sep="\t")
    loc_map = mouselocal.drop_duplicates("gene_id").set_index("gene_id")["type"]

    out = df.copy()
    out["localization"] = df.index.str.split(".").str[0].map(loc_map).values
    return out


# ---------------------------------------------------------------------------
# 5. 그림 — 노트북의 두 빈 코드 셀에 대응
# ---------------------------------------------------------------------------

def plot_result1(df: pd.DataFrame, out_path: Path) -> float:
    """**과제 1 (셀 33)** — Figure 4D 스타일 scatter + Pearson r."""
    r, _ = pearsonr(df["clip_enrich"], df["rden_change"])

    fig, ax = plt.subplots(figsize=(4.0, 4.0))
    ax.scatter(df["clip_enrich"], df["rden_change"],
               s=4, c="black", alpha=0.25, edgecolors="none")
    ax.axhline(0, color="lightgray", lw=0.5)
    ax.axvline(0, color="lightgray", lw=0.5)
    ax.set_xlim(-6, 4)
    ax.set_ylim(-2, 2)
    ax.set_xlabel("LIN28A CLIP enrichment (log$_2$)")
    ax.set_ylabel(
        "Ribosome density change\nupon $\\it{Lin28a}$ knockdown (log$_2$)"
    )
    ax.set_title(
        "CLIP and ribosome footprinting\nupon $\\it{Lin28a}$ knockdown",
        fontsize=10, loc="left",
    )
    ax.text(0.97, 0.06, f"r = {r:.4f}", transform=ax.transAxes,
            ha="right", va="bottom", fontsize=10)
    plt.tight_layout()
    fig.savefig(out_path, dpi=200)
    plt.close(fig)
    return r


def plot_result2(df_full: pd.DataFrame, out_path: Path) -> None:
    """**과제 2 (셀 37)** — localization 색칠 scatter (Figure S6A 스타일).

    전체 gene 을 회색 배경으로 깔고, localization 이 매핑된 gene 만
    nucleus / integral membrane / cytoplasm 색으로 표시.
    """
    fig, ax = plt.subplots(figsize=(4.8, 4.4))

    unl = df_full["localization"].isna()
    ax.scatter(df_full.loc[unl, "clip_enrich"], df_full.loc[unl, "rden_change"],
               s=3, c="lightgray", alpha=0.35, edgecolors="none")

    for typ in ("nucleus", "integral membrane", "cytoplasm"):
        sel = df_full["localization"] == typ
        ax.scatter(df_full.loc[sel, "clip_enrich"],
                   df_full.loc[sel, "rden_change"],
                   s=6, c=LOC_COLORS[typ], alpha=0.75,
                   edgecolors="none", label=typ)

    ax.axhline(0, color="lightgray", lw=0.5)
    ax.axvline(0, color="lightgray", lw=0.5)
    ax.set_xlim(-6, 4)
    ax.set_ylim(-2, 3)
    ax.set_xlabel("LIN28A CLIP enrichment (log$_2$)")
    ax.set_ylabel(
        "Ribosome density change\nupon $\\it{Lin28a}$ knockdown (log$_2$)"
    )
    ax.set_title("Linkage to localization", fontsize=10, loc="left")
    ax.legend(loc="upper left", fontsize=8, framealpha=0.95,
              markerscale=1.4, handletextpad=0.4)

    plt.tight_layout()
    fig.savefig(out_path, dpi=200)
    plt.close(fig)


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def main() -> None:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    # 1) 데이터 수집 — 원본 노트북 셀 9 / 11
    setup_datapack()
    # 2) 작업 폴더 — 원본 노트북 셀 22 / 24
    setup_workdir()
    # 3) 카운팅 — 원본 노트북 셀 26
    counts_path = run_feature_counts()
    # 4) 카운트 로드 — 원본 노트북 셀 28
    cnts = pd.read_csv(counts_path, sep="\t", comment="#", index_col=0)
    print(f"loaded {len(cnts):,} genes from {counts_path.name}")

    # 5) 과제 1 — 셀 33: Figure 4D 스타일
    df_full = compute_metrics(cnts, min_reads=30)
    print(f"[result1] n = {len(df_full):,} genes (min_reads=30)")
    r = plot_result1(df_full, RESULTS_DIR / "result1.png")
    print(f"[result1] Pearson r = {r:.4f}  (논문값 0.4028)")

    # 6) 과제 2 — 셀 37: localization 색칠 scatter (Figure S6A 스타일)
    df_full_loc = annotate_localization(df_full)
    plot_result2(df_full_loc, RESULTS_DIR / "result2.png")
    n_loc = df_full_loc["localization"].notna().sum()
    print(f"[result2] {n_loc:,} localization-annotated / "
          f"{len(df_full_loc):,} total")

    print("\n[localization 별 평균값]")
    print(df_full_loc.groupby("localization")[
        ["clip_enrich", "rden_change"]
    ].mean().round(3))

    print(f"\nDone. Outputs in {RESULTS_DIR}")
    for p in sorted(RESULTS_DIR.glob("*.png")):
        print(f"  - {p.name}")


if __name__ == "__main__":
    main()
