PIP=pip3

.PHONY: docs  # necessary so it doesn't look for 'docs/makefile html'

init:
	$(PIP) install pipenv --upgrade
	pipenv install --dev --skip-lock


pylint:
	pipenv run pylint -E pyrxn


benchmark:
	rm -rf .benchmarks/images/*svg
	python -m pytest --benchmark-autosave --benchmark-max-time=0.1 --benchmark-group-by=func


lock:
	pipenv lock
	pipenv lock -r > requirements.txt
	pipenv lock -r > requirements-dev.txt
	pipenv lock --dev -r >> requirements-dev.txt
