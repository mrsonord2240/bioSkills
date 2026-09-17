#!/bin/bash
# Smoke test: confirms an LDSC install actually works before pointing it at real GWAS data.
# Exercises --h2, --rg, and --h2-cts (Finucane 2018 cell-type prioritization) against the
# tiny bundled fixtures in data/ (see data/README.md for provenance).
#
# Prerequisite: git clone https://github.com/CBIIT/ldsc.git, Python 3.9+ env per SKILL.md's
# "Tool Install Notes", and the one-line ldscore/sumstats.py patch (.loc[:,1:] -> .iloc[:,1:])
# for the --h2-cts step below.
#
# Usage:
#   LDSC_DIR=/path/to/ldsc bash smoke_test_ldsc.sh

set -euo pipefail

LDSC_DIR=${LDSC_DIR:?Set LDSC_DIR to your CBIIT/ldsc clone, e.g. LDSC_DIR=./ldsc bash smoke_test_ldsc.sh}
DATA=$(dirname "$0")/data
OUT=$(mktemp -d)

echo "=== munge_sumstats.py ==="
python "${LDSC_DIR}/munge_sumstats.py" \
    --sumstats "${DATA}/munge_sumstats_input" \
    --merge-alleles "${DATA}/merge_alleles" \
    --N 6702 --signed-sumstats OR,1 \
    --out "${OUT}/munge"
# Expect: "1 SNPs remain" and a printed Mean chi^2 / Lambda GC -- confirms munging runs.

echo "=== --h2 (total heritability) ==="
python "${LDSC_DIR}/ldsc.py" --h2 "${DATA}/sumstats/0" \
    --ref-ld "${DATA}/ldscore/oneld_onefile1" \
    --w-ld "${DATA}/ldscore/w1" \
    --out "${OUT}/h2"
grep "Total Observed scale h2" "${OUT}/h2.log"
# Expect: Total Observed scale h2 ~ 0.38 (0.04), Intercept ~2.09, Ratio ~0.12.

echo "=== --rg (cross-trait genetic correlation) ==="
python "${LDSC_DIR}/ldsc.py" --rg "${DATA}/sumstats/0,${DATA}/sumstats/1" \
    --ref-ld "${DATA}/ldscore/twold_onefile1" \
    --w-ld "${DATA}/ldscore/w1" \
    --out "${OUT}/rg"
grep "Genetic Correlation" "${OUT}/rg.log"
# Expect: Genetic Correlation ~ 0.11 (0.08), P ~ 0.15 (these two simulated traits are ~uncorrelated).

echo "=== --h2-cts (Finucane 2018 cell-type prioritization) ==="
echo "(requires the ldscore/sumstats.py .loc[:,1:] -> .iloc[:,1:] patch -- see SKILL.md)"
# test.ldcts' rows are relative to $DATA (e.g. "ldscore/twold_firstfile"), so run from there.
( cd "${DATA}" && python "${LDSC_DIR}/ldsc.py" --h2-cts "sumstats/0" \
    --ref-ld-chr "ldscore/oneld_onefile" \
    --ref-ld-chr-cts "test.ldcts" \
    --w-ld-chr "ldscore/w" \
    --out "${OUT}/cts" )
cat "${OUT}/cts.cell_type_results.txt"
# Expect a two-row table (CellTypeA, CellTypeB) with Coefficient / SE / P-value -- one row should
# land well under a Bonferroni threshold of 0.05/2 = 0.025, the other should not (these are
# synthetic annotations, not real tissues, so there is no "correct" tissue -- the point is that
# the two rows differ and neither errors out).

echo "All three LDSC entry points ran. Output in ${OUT}/"
