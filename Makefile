# Reproduction targets. Run from the repository root.
PY ?= python3

check:      ## fast self-checks (code construction, decoder, derived numbers, figures, one bit-reproducible MC point)
	$(PY) -m pytest -q tests -m "not slow"

check-all:  ## also the slower validations (L=2 exact master equation, OU sigma->0 self-test)
	$(PY) -m pytest -q tests

numbers:    ## regenerate results/derived_numbers.json from results/*.json
	cd scripts && $(PY) protocol_f.py

figs:       ## regenerate paper/figs/*.pdf from results/*.json
	$(PY) paper/figs/make_figs.py

clean:
	rm -rf scripts/__pycache__ tests/__pycache__ .pytest_cache results/*_pytest*.json

.PHONY: check check-all numbers figs clean
