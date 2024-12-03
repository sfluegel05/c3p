RUN = poetry run

test: doctest pytest

pytest:
	poetry run pytest

DOCTEST_DIR = c3p
doctest:
	find $(DOCTEST_DIR) -not -path "c3p/programs/*.py" -type f \( -name "*.rst" -o -name "*.md" -o -name "*.py" \) -print0 | xargs -0 $(RUN) python -m doctest --option ELLIPSIS --option NORMALIZE_WHITESPACE


%-doctest: %
	$(RUN) python -m doctest --option ELLIPSIS --option NORMALIZE_WHITESPACE $<

