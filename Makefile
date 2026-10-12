.PHONY: all lint test install dev clean distclean html serve

PYTHON ?= python

all: ;

lint:
	q2lint
	flake8

test: all
	py.test

install: all
	$(PYTHON) -m pip install -v .

dev: all
	pip install -e .

# Documentation (q2doc + MyST, the gut-to-soil-tutorial layout). Run after
# `q2doc refresh-cache` in an environment where the plugin is installed.
html:
	cp environment-files/q2-confluence-qiime2-dev.yml docs/_static/environment.yml
	cd docs && Q2DOC_BASE_COMMAND='qiime' jupyter book build --html

serve:
	npx serve docs/_build/html/ -p 4000

clean: distclean
	rm -rf docs/_build/html/

distclean: ;
