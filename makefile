# Define the package name, module name, environment name, and Python version
PACKAGE_NAME = AquaSEBS
MODULE_NAME = $(subst -,_,$(PACKAGE_NAME))
ENVIRONMENT_NAME = $(PACKAGE_NAME)
PYTHON_VERSION = $(if $(PYTHON),$(PYTHON),3.14)

# Allow passing extra flags via command line: make <target> ARGS="-k test_name"
ARGS ?=

# Allow passing VERSION via make set-version VERSION=x.y.z or make set-version x.y.z
ifeq (set-version,$(firstword $(MAKECMDGOALS)))
  VERSION_ARG := $(wordlist 2,$(words $(MAKECMDGOALS)),$(MAKECMDGOALS))
  $(eval $(VERSION_ARG):;@:)
endif
VERSION ?= $(VERSION_ARG)

# Include a Makefile.local to override things, if found
-include Makefile.local

.DEFAULT_GOAL := help

.PHONY: help clean test wheel sdist build twine-upload dist set-version \
        install uninstall reinstall environment remove-environment \
        repomix generate-history show-issues

help: ## Show available targets
	@grep -E '^[a-zA-Z0-9_-]+:.*?## .*$$' $(MAKEFILE_LIST) | \
		awk 'BEGIN {FS = ":.*?## "}; {printf "  \033[36m%-28s\033[0m %s\n", $$1, $$2}'

##@ Cleanup & Environment
clean: ## Remove build artifacts, caches, coverage files, and logs
	rm -rf *.o *.out *.log
	rm -rf build/
	rm -rf dist/
	rm -rf *.egg-info
	rm -rf .pytest_cache
	rm -rf htmlcov/ .coverage
	rm -rf $(PACKAGE_NAME)-repomix.xml repomix-output.xml
	find . -type d -name "__pycache__" -exec rm -rf {} +
	find . -type f -name "*.py[co]" -delete

environment: ## Create a new Conda/Mamba environment (override with PYTHON=3.12)
	mamba create -y -n $(ENVIRONMENT_NAME) -c conda-forge python=$(PYTHON_VERSION)

remove-environment: ## Remove the Conda/Mamba environment
	mamba env remove -y -n $(ENVIRONMENT_NAME)

install: ## Install the package in editable mode with development dependencies
	pip install -e .[dev]

uninstall: ## Uninstall the package
	pip uninstall -y $(PACKAGE_NAME)

reinstall: ## Reinstall the package (uninstall and then install)
	$(MAKE) uninstall
	$(MAKE) install

##@ Testing
test: ## Run Python unit tests using pytest
	pytest -vv -p no:warnings $(ARGS)

##@ Packaging & Release
set-version: ## Set package version across pyproject.toml, CITATION.cff, and README.md (usage: make set-version VERSION=1.3.2)
	@test -n "$(VERSION)" || (echo "Error: VERSION is required (e.g., make set-version VERSION=1.3.2)" && exit 1)
	@python -c 'import re, sys; from datetime import date; from pathlib import Path; \
		v = sys.argv[1]; \
		today = date.today().isoformat(); \
		p = Path("pyproject.toml"); \
		p.write_text(re.sub(r"(?m)^(version\s*=\s*\")[^\"]+(\")", rf"\g<1>{v}\g<2>", p.read_text())); \
		c = Path("CITATION.cff"); \
		ct = re.sub(r"(?m)^(\s*version:\s*)\S+", rf"\g<1>{v}", c.read_text()); \
		ct = re.sub(r"(?m)^(date-released:\s*)\S+", rf"\g<1>{today}", ct); \
		c.write_text(ct); \
		r = Path("README.md"); \
		rt = r.read_text(); \
		rt = re.sub(r"(version=\{)[^}]+(\})", rf"\g<1>{v}\g<2>", rt); \
		rt = re.sub(r"(Version `)[^`]+(`)", rf"\g<1>{v}\g<2>", rt); \
		rt = re.sub(r"(Implemented Code \(`v)[^`]+(`\))", rf"\g<1>{v}\g<2>", rt); \
		r.write_text(rt); \
		print(f"Updated version to {v} (date-released: {today}) in pyproject.toml, CITATION.cff, and README.md")' "$(VERSION)"

wheel: ## Build the Python wheel
	python -m build --wheel

sdist: ## Build the Python source distribution
	python -m build --sdist

build: ## Build the Python wheel and source distribution
	python -m build

twine-upload: ## Upload the built package to PyPI using Twine
	twine upload dist/*

dist: ## Clean, build, and upload the package to PyPI
	$(MAKE) clean
	$(MAKE) build
	$(MAKE) twine-upload

##@ Utilities
repomix: ## Pack repository context for LLM review
	npx repomix -o $(PACKAGE_NAME)-repomix.xml --ignore "**/*.tif,**/*.jpeg,**/*.aux.xml,**/*.csv,**/*.pdf,**/*.svg,**/*.ipynb,*-repomix*.xml,repomix-output.xml,release_diffs.patch" --copy

generate-history: ## Generate release commit history patch
	fish -c ' \
		git fetch --tags --quiet; \
		set first_commit (git rev-list --max-parents=0 HEAD | head -n 1); \
		set releases (gh release list --limit 1000 --json tagName,createdAt --jq "sort_by(.createdAt) | .[].tagName"); \
		begin; \
			echo "========================================"; \
			echo "REPOSITORY RELEASE TIMELINE"; \
			echo "========================================"; \
			gh release list --limit 1000; \
			echo ""; \
			if test (count $$releases) -ge 1; \
				set first_rel $$releases[1]; \
				echo "========================================"; \
				echo "INITIAL COMMIT ($$first_commit) -> FIRST RELEASE ($$first_rel)"; \
				echo "========================================"; \
				git log --reverse --stat -p $$first_rel -- "*.py" "pyproject.toml"; \
				echo ""; \
			end; \
			if test (count $$releases) -ge 2; \
				for i in (seq 2 (count $$releases)); \
					set prev_idx (math $$i - 1); \
					set prev $$releases[$$prev_idx]; \
					set curr $$releases[$$i]; \
					echo "========================================"; \
					echo "RELEASE COMMIT HISTORY: $$prev -> $$curr"; \
					echo "========================================"; \
					git log --reverse --stat -p $$prev..$$curr -- "*.py" "pyproject.toml"; \
					echo ""; \
				end; \
			end; \
			if test (count $$releases) -ge 1; \
				set latest_rel $$releases[-1]; \
				echo "========================================"; \
				echo "LATEST RELEASE ($$latest_rel) -> CURRENT WORKING TREE (HEAD + UNCOMMITTED)"; \
				echo "========================================"; \
				git log --reverse --stat -p $$latest_rel..HEAD -- "*.py" "pyproject.toml"; \
				echo ""; \
				git status -s -- "*.py" "pyproject.toml"; \
				echo ""; \
				git diff HEAD -- "*.py" "pyproject.toml"; \
			else; \
				echo "========================================"; \
				echo "INITIAL COMMIT ($$first_commit) -> CURRENT WORKING TREE"; \
				echo "========================================"; \
				git log --reverse --stat -p -- "*.py" "pyproject.toml"; \
				echo ""; \
				git status -s -- "*.py" "pyproject.toml"; \
				echo ""; \
				git diff HEAD -- "*.py" "pyproject.toml"; \
			end; \
		end > release_diffs.patch; \
	'

show-issues: ## List open GitHub issues
	gh issue list --state open --json number,title,body --template '{{range .}}## #{{.number}}: {{.title}}{{"\n\n"}}{{.body}}{{"\n\n---\n\n"}}{{end}}'