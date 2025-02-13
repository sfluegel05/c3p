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
VERBOSE = 1
EXPTS = claude-3-sonnet-undef gpt-4o-undef gpt-4o-hi gpt-4o-nodefs ensembl-5 o1-preview-undef deepseek-coder-undef
LEARN = $(RUN) c3p-learn --verbose $(VERBOSE)

2025-llama-405b:
	$(LEARN)  -m lbl/llama-3-405b --dataset $(DATASET) -w results/2025 -x -

2025-llama:
	$(LEARN)  -m lbl/llama --dataset $(DATASET) -w results/2025 -x -

2025-claude-sonnet:
	$(RUN) c3p-learn  --verbose $(VERBOSE)  -m claude-3-sonnet --dataset $(DATASET) -w results/2025 -x -

2025-claude-sonnet-use-the-force:
	$(LEARN) -e force  -m claude-3-sonnet --use-the-force --dataset $(DATASET) -w results/2025 -x -

2025-gpt-4o:
	$(LEARN) -m gpt-4o --dataset $(DATASET) -w results/2025 -x -

2025-o3-mini:
	$(LEARN) -m o3-mini --dataset $(DATASET) -w results/2025 -x -

2025-gemini2-flash:
	$(LEARN) -m gemini-2.0-flash-exp --dataset $(DATASET) -w results/2025 -x -

2025-deepseek-coder:
	$(LEARN) -m deepseek-coder --randomize-order --dataset $(DATASET) -w results/2025 -x -

2025-deepseek-coder-use-the-force:
	$(LEARN) -e force  -m deepseek-coder --use-the-force --dataset $(DATASET) -w results/2025 -x -

2025-cborg-deepthought:
	$(LEARN) -m cborg-deepthought --dataset $(DATASET) -w results/2025 -x -

2025-o1:
	$(LEARN) -m o1 --dataset $(DATASET) -w results/2025 -x -

2025-gpt-4o-hi:
	$(LEARN)  -e hi -m gpt-4o -f 0.9 -a 6 --dataset $(DATASET) -w results/2025 -x -

2025-deepseek-hi:
	$(LEARN)  -e hi -m deepseek-coder -f 0.9 -a 6 --dataset $(DATASET) -w results/2025 -x -X "carbamate ester" -

2025-deepseek-reasoner:
	$(LEARN)  -m deepseek-reasoner --dataset $(DATASET) -w results/2025 -x  -X "carbamate ester" -X "peptide antibiotic" -X "triglyceride" -

2025-together-r1:
	$(LEARN)  -m deepseek-ai/DeepSeek-R1 --dataset $(DATASET) -w results/2025 -x -

2025-deepseek-reasoner-use-the-force:
	$(LEARN) -e force  -m deepseek-reasoner --use-the-force --dataset $(DATASET) -w results/2025 -x  -X "carbamate ester" -

2025-gpt-4o-nodefs:
	$(RUN) c3p-learn -e nodefs -m gpt-4o --exclude-definitions --dataset $(DATASET) -w results/2025 -x -

2025-r1-distill:
	$(LEARN)  -m groq/deepseek-r1-distill-llama-70b --dataset $(DATASET) -w results/2025 -x -


# Create ensemble model

create-ensemble:
	$(RUN) c3p-combine --verbose 2 results/2025/* -o results/2025/ensembl-5 -m 5

create-ensemble-6:
	$(RUN) c3p-combine --verbose 2 results/2025/* -o results/2025/ensemble-6 -m 6

create-ensemble-7:
	$(RUN) c3p-combine --verbose 2 results/2025/* -o results/2025/ensemble-7 -m 7

create-ensemble-3:
	$(RUN) c3p-combine --verbose 2 results/2025/* -o results/2025/ensemble-3 -m 3

create-ensemble-4:
	$(RUN) c3p-combine --verbose 2 results/2025/* -o results/2025/ensemble-4 -m 4

eval:
	$(RUN) c3p-validate --verbose 1 -w results/2025/ensembl-5/cache -d $(DATASET)

all-validate: $(patsubst %, validate-%, $(EXPTS))
validate-%:
	$(RUN) c3p-validate --verbose 1 -w results/2025/$*/cache -d $(DATASET)

summarize-%:
	$(RUN) c3p-summarize --verbose 1 -w results/2025/$*/eval -o results/2025/$*/summary

EX_CORE = claude-3-sonnet-undef o1-undef claude-3-sonnet-force  gpt-4o-undef gpt-4o-hi gemini-2.0-flash-exp-undef o3-mini-undef ensemble-7
EX_ORIGINAL = claude-3-sonnet-undef  gpt-4o-undef o1-preview-undef gpt-4o-nodefs gpt-4o-hi ensembl-5 deepseek-coder-undef
EX_COMPLETED = $(EX_ORIGINAL) 
PARTIAL =  gpt-4o-hi o1-preview-undef deepseek-coder-undef deepseek-reasoner-undef deepseek-reasoner-force deepseek-chat-hi ensembl-7

compare:
	$(RUN) c3p-compare -o results/2025/comparison $(patsubst %,results/2025/%, $(EX_CORE))

compare-all:
	$(RUN) c3p-compare -x -o results/2025/comparison $(patsubst %,results/2025/%, $(EX_COMPLETED) $(PARTIAL) chebifier)

compare-chebifier:
	$(RUN) c3p-compare -x -o results/2025/comparison-chebifier $(patsubst %,results/2025/%, $(EX_CORE) chebifier)


compare-llama:
	$(RUN) c3p-compare -x -o results/2025/comparison-partial results/2025/claude-3-sonnet-undef  results/2025/gpt-4o-undef results/2025/o1-preview-undef results/2025/gpt-4o-nodefs results/2025/gpt-4o-hi results/2025/ensembl-5 results/2025/deepseek-coder-undef results/2025/llama-undef


compare-reasoners:
	$(RUN) c3p-compare -x -o results/2025/comparison-reasoners results/2025/o1-undef results/2025/o3-mini-undef results/2025/claude-3-sonnet-undef results/2025/deepseek-reasoner-undef 

compare-force:
	$(RUN) c3p-compare -x -o results/2025/comparison-force results/2025/claude-3-sonnet-undef results/2025/claude-3-sonnet-force

compare-reasoner-force:
	$(RUN) c3p-compare -x -o results/2025/comparison-reasoner-force results/2025/deepseek-reasoner-undef results/2025/deepseek-reasoner-force

compare-deepseek-hi:
	$(RUN) c3p-compare -x -o results/2025/comparison-reasoner-hi results/2025/deepseek-coder-undef results/2025/deepseek-coder-hi

