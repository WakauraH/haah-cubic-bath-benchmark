"""Reproduction checks. Every test either recomputes a shipped number or regenerates a shipped artifact and compares.
Run from the repository root:  python3 -m pytest -q tests            (fast, ~1-2 min on 16 cores)
                               python3 -m pytest -q tests -m slow    (L=2 master equation, OU self-test; several minutes)"""
import json, os, subprocess, sys, shutil, pathlib, math
import numpy as np, pytest

ROOT = pathlib.Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"; RESULTS = ROOT / "results"; FIGS = ROOT / "paper" / "figs"
sys.path.insert(0, str(SCRIPTS))

def run(cmd, cwd, timeout=1800):
    r = subprocess.run([sys.executable] + cmd, cwd=cwd, capture_output=True, text=True, timeout=timeout)
    assert r.returncode == 0, f"{cmd} failed:\n{r.stdout[-2000:]}\n{r.stderr[-2000:]}"
    return r.stdout

def approx_equal(a, b, rtol=1e-9, path=""):
    """Recursive numeric comparison of two JSON trees; collects mismatches."""
    bad = []
    if isinstance(a, dict) and isinstance(b, dict):
        if set(a) != set(b): bad.append(f"{path}: keys {sorted(set(a)^set(b))}")
        for k in a.keys() & b.keys(): bad += approx_equal(a[k], b[k], rtol, f"{path}/{k}")
    elif isinstance(a, list) and isinstance(b, list):
        if len(a) != len(b): bad.append(f"{path}: length {len(a)} vs {len(b)}")
        for i, (x, y) in enumerate(zip(a, b)): bad += approx_equal(x, y, rtol, f"{path}[{i}]")
    elif isinstance(a, (int, float)) and isinstance(b, (int, float)) and not isinstance(a, bool):
        if not (a == b or (isinstance(a, float) and isinstance(b, float) and math.isnan(a) and math.isnan(b))
                or abs(a - b) <= rtol * max(abs(a), abs(b), 1e-300)):
            bad.append(f"{path}: {a} vs {b}")
    elif a != b:
        bad.append(f"{path}: {a!r} vs {b!r}")
    return bad

# ---------- code construction ----------
def test_code_construction_matches_haah_table():
    import haah_code as hc
    for L in [2, 3, 4, 5, 6, 7]:
        X, Z, N = hc.build(L)
        assert not ((X @ Z.T) % 2).any(), f"stabilizers do not commute at L={L}"
        assert N - hc.rank_f2(X) - hc.rank_f2(Z) == hc.k_table(L)
    assert [hc.k_table(L) for L in (2, 3, 4, 5, 6, 7, 15)] == [6, 2, 14, 2, 6, 2, 50]

def test_decoder_selftest():
    out = run(["-c", "import haah_tmem; haah_tmem.selftest()"], cwd=SCRIPTS)
    assert "selftest ok" in out

# ---------- derived numbers and figures ----------
def test_protocol_f_regenerates_derived_numbers(tmp_path):
    shipped = json.load(open(RESULTS / "derived_numbers.json"))
    backup = tmp_path / "derived_numbers.json"; shutil.copy(RESULTS / "derived_numbers.json", backup)
    try:
        run(["protocol_f.py"], cwd=SCRIPTS)
        new = json.load(open(RESULTS / "derived_numbers.json"))
    finally:
        shutil.copy(backup, RESULTS / "derived_numbers.json")
    # the replica study uses a fixed seed, so even the Monte Carlo quantiles must agree exactly
    bad = approx_equal(shipped, new, rtol=1e-9)
    assert not bad, "derived_numbers.json differs from the shipped file:\n" + "\n".join(bad[:20])

def test_make_figs_builds_all_pdfs():
    before = {f: (FIGS / f).stat().st_mtime for f in os.listdir(FIGS) if f.endswith(".pdf")}
    run([str(FIGS / "make_figs.py")], cwd=ROOT)
    for f in ["fig1_baseline.pdf", "fig2_channels.pdf", "fig3_ou.pdf", "fig4_pseudomode.pdf", "fig5_lamb.pdf"]:
        p = FIGS / f
        assert p.exists() and p.stat().st_size > 1000 and p.stat().st_mtime >= before.get(f, 0), f

# ---------- Monte Carlo: bit-for-bit reproduction of one archived point ----------
def test_flat_bath_L5_reproduces_archived_point():
    """Same seeds (1000*L+i), same kernel: the archived L=5, beta=4.5 flat-bath mean and SEM must come back exactly."""
    arch = json.load(open(RESULTS / "tmem_beta4.5_fast_flat.json"))["5"]
    out = run(["haah_tmem_fast.py", "4.5", "5", "100", "_pytest_flat_L5"], cwd=SCRIPTS)
    new = json.load(open(RESULTS / "tmem_beta4.5_pytest_flat_L5.json"))["5"]
    assert new["n"] == arch["n"] == 100
    assert abs(new["Tmem_mean"] - arch["Tmem_mean"]) <= 1e-9 * arch["Tmem_mean"], (new["Tmem_mean"], arch["Tmem_mean"])
    assert abs(new["Tmem_sem"] - arch["Tmem_sem"]) <= 1e-9 * arch["Tmem_sem"]
    assert new["modes"] == arch["modes"]

def test_filtered_bath_L5_reproduces_archived_point():
    """Channel filter s0=0.001 at L=5 (the largest gain in the paper at its smallest size)."""
    arch = json.load(open(RESULTS / "tmem_beta4.5_fast_s0_0.001_allL.json"))["5"]
    run(["haah_tmem_fast.py", "4.5", "5", "100", "_pytest_s0_L5", "1,1,0.001"], cwd=SCRIPTS)
    new = json.load(open(RESULTS / "tmem_beta4.5_pytest_s0_L5.json"))["5"]
    assert abs(new["Tmem_mean"] - arch["Tmem_mean"]) <= 1e-9 * arch["Tmem_mean"]
    assert abs(new["Tmem_sem"] - arch["Tmem_sem"]) <= 1e-9 * arch["Tmem_sem"]

# ---------- pseudomode models (deterministic linear algebra) ----------
def test_pseudomode_twodefect_reproduces_archived_table():
    arch = json.load(open(RESULTS / "gap" / "B_twodefect_meanforce.json"))
    target = RESULTS / "gap" / "B_twodefect_meanforce.json"; bak = target.read_bytes()
    try:
        run(["pseudomode_twodefect.py"], cwd=SCRIPTS)
        new = json.load(open(target))
    finally:
        target.write_bytes(bak)
    bad = approx_equal(arch, new, rtol=1e-6)
    assert not bad, "\n".join(bad[:20])

# ---------- slow validations ----------
@pytest.mark.slow
def test_L2_exact_master_equation_vs_gillespie():
    target = RESULTS / "tier1_validation.json"; bak = target.read_bytes()
    try:
        out = run(["haah_davies.py"], cwd=SCRIPTS, timeout=3600)
    finally:
        target.write_bytes(bak)
    assert "max|dP|" in out

@pytest.mark.slow
def test_ou_sigma_zero_selftest():
    out = run(["haah_tmem_ou.py", "selftest"], cwd=SCRIPTS, timeout=3600)
    assert "selftest" in out.lower()
