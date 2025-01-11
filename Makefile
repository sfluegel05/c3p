RUN = poetry run

test: doctest pytest

pytest:
	poetry run pytest tests

DOCTEST_DIR = c3p
doctest:
	find $(DOCTEST_DIR) -not -path "c3p/programs/*.py" -type f \( -name "*.rst" -o -name "*.md" -o -name "*.py" \) -print0 | xargs -0 $(RUN) python -m doctest --option ELLIPSIS --option NORMALIZE_WHITESPACE


%-doctest: %
	$(RUN) python -m doctest --option ELLIPSIS --option NORMALIZE_WHITESPACE $<

# Experiments

DATASET = results/2025/benchmark/dataset.json
VERBOSE = 0
EXPTS = claude-3-sonnet-undef gpt-4o-undef gpt-4o-hi gpt-4o-nodefs ensembl-5 o1-preview-undef deepseek-coder-undef

2025-llama-405b:
	$(RUN) c3p-learn --verbose $(VERBOSE)  -m lbl/llama-3-405b --dataset $(DATASET) -w results/2025 -x -

2025-llama:
	$(RUN) c3p-learn --verbose $(VERBOSE)  -m lbl/llama --dataset $(DATASET) -w results/2025 -x -

2025-claude-sonnet:
	$(RUN) c3p-learn  -m lbl/claude-3-sonnet --dataset $(DATASET) -w results/2025 -x -
#	$(RUN) c3p-learn  -m lbl/claude-3-sonnet --dataset $(DATASET) -w results/2025 -x -X "16beta-hydroxy steroid" -

2025-gpt-4o:
	$(RUN) c3p-learn -m gpt-4o --dataset $(DATASET) -w results/2025 -x -

2025-deepseek-coder:
	$(RUN) c3p-learn --verbose $(VERBOSE) -m deepseek-coder --dataset $(DATASET) -w results/2025 -x -X "11beta-hydroxy steroid" -

2025-o1-preview:
	$(RUN) c3p-learn -m o1-preview --dataset $(DATASET) -w results/2025 -x -

2025-gpt-4o-hi:
	$(RUN) c3p-learn  --verbose $(VERBOSE)  -e hi -m gpt-4o -f 0.9 -a 6 --dataset $(DATASET) -w results/2025 -x -

2025-gpt-4o-nodefs:
	$(RUN) c3p-learn -e nodefs -m gpt-4o --exclude-definitions --dataset $(DATASET) -w results/2025 -x -

# Create ensemble model

create-ensemble:
	$(RUN) c3p-combine --verbose 2 results/2025/* -o results/2025/ensembl-5 -m 5

eval:
	$(RUN) c3p-validate --verbose 1 -w results/2025/ensembl-5/cache -d $(DATASET)

all-validate: $(patsubst %, validate-%, $(EXPTS))
validate-%:
	$(RUN) c3p-validate --verbose 1 -w results/2025/$*/cache -d $(DATASET)

summarize-%:
	$(RUN) c3p-summarize --verbose 1 -w results/2025/$*/eval -o results/2025/$*/summary

compare:
	$(RUN) c3p-compare -o results/2025/comparison results/2025/claude-3-sonnet-undef  results/2025/gpt-4o-undef results/2025/o1-preview-undef results/2025/gpt-4o-nodefs results/2025/gpt-4o-hi results/2025/ensembl-5 results/2025/deepseek-coder-undef


compare-partial:
	$(RUN) c3p-compare -x -o results/2025/comparison-partial results/2025/claude-3-sonnet-undef  results/2025/gpt-4o-undef results/2025/o1-preview-undef results/2025/gpt-4o-nodefs results/2025/gpt-4o-hi results/2025/ensembl-5 results/2025/deepseek-coder-undef results/2025/llama-undef
