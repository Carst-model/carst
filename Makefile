.PHONY: test lint clean diamond
PROJECT_DIR := $(shell pwd -P)
FIREDRAKE_VENV_FULL := $(PROJECT_DIR)/firedrake
SPUD_CHECK := firedrake/usr/include/spud

FIREDRAKE_ACTIVATION := firedrake/bin/activate
LIBSPUD_CONFIG_FLAGS := --prefix=$(FIREDRAKE_VENV_FULL)/usr

all: diamond_default.rng $(SPUD_CHECK) firedrake
	{ \
		set -e; \
		source $(FIREDRAKE_ACTIVATION); \
		diamond -s diamond_defaut.rng; \
	}

test: firedrake carst basic_tests.py scripts
	@echo "Running simulation"
	{ \
		source $(FIREDRAKE_ACTIVATION); \
		python3 basic_tests.py; \
	}
	@echo "Running visualisation"
	{ \
		source $(FIREDRAKE_ACTIVATION); \
		python3 scripts/stratiMesh.py; \
	}

lint:
	@echo "Linting carst codebase"
	@flake8 carst
	@echo "Linting carst test suite"
	@flake8 test
	@echo "Linting carst examples"
	@flake8 examples
	@echo "Linting carst scripts"
	@flake8 scripts

clean:
	@echo "Removing venv..."
	-rm -r firedrake firedrake-install

diamond_default.rng: $(SPUD_CHECK) diamond_defaut.rnc
	@echo "Building schema..."
	{ \
		set -e; \
		source $(FIREDRAKE_ACTIVATION); \
		spud-preprocess diamond_defaut.rnc; \
	}

planning_diagram.png:
	@echo "Building diagram..."
	ditaa planning_diagram.ditaa

firedrake:
	@echo "Building lastest firedrake version"
	curl -O https://raw.githubusercontent.com/firedrakeproject/firedrake/master/scripts/firedrake-install
	python3 firedrake-install --no-package-manager


$(SPUD_CHECK): firedrake spud
	@echo "Building spud..."
	{ \
		set -e; \
		cd spud; \
		./configure $(LIBSPUD_CONFIG_FLAGS); \
		make; \
		make install; \
	}
