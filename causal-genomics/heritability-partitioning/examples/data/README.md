# Smoke-test fixtures

Everything under `sumstats/`, `ldscore/`, `merge_alleles`, and `munge_sumstats_input` is copied
verbatim from `CBIIT/ldsc`'s own unit-test suite (`test/simulate_test/`, `test/munge_test/`,
commit `1f09cf0c`) -- the tool's own CI fixtures, not real GWAS data. Two simulated 1000-SNP
Z-score panels (`sumstats/0`, `sumstats/1`) and chromosome-split (chr1/chr2), single- and
two-category LD score files. `test.ldcts` was authored for this Skill (not present upstream) to
exercise the Finucane 2018 `--h2-cts` code path, pointing at the two bundled `twold_firstfile` /
`twold_secondfile` two-category LD scores as stand-in "cell types".

`ldscore/w{1,2}.l2.ldscore.gz` are the same single regression-weight file duplicated per
chromosome -- CBIIT/ldsc ships only one `w.l2.ldscore` (not chromosome-split); duplicating it is a
smoke-test convenience, not a real ancestry-matched weights reference.

Run `bash ../smoke_test_ldsc.sh` from `examples/` to exercise `--h2`, `--rg`, and `--h2-cts`
against these fixtures and confirm your LDSC install works before pointing it at real GWAS data.
