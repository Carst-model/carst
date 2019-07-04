.PHONY: test lint clean

# Paths for setup
PROJECT_DIR := $(shell pwd -P)
FIREDRAKE_VENV_FULL := $(PROJECT_DIR)/firedrake

# Paths to alias long targets
DIAMOND_BIN := firedrake/usr/bin/diamond
DXDIFF_BIN := firedrake/usr/bin/dxdiff

# For sourcing the shell file to activate the firedrake venv
FIREDRAKE_ACTIVATION := firedrake/bin/activate

# Temporary test
all: diamond_default.rng $(DIAMOND_BIN) firedrake
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

diamond_default.rng: $(DIAMOND_BIN) diamond_default.rnc
	@echo "Building schema..."
	{ \
		set -e; \
		source $(FIREDRAKE_ACTIVATION); \
		spud-preprocess diamond_default.rnc; \
	}

planning_diagram.png:
	@echo "Building diagram..."
	ditaa planning_diagram.ditaa

firedrake:
	@echo "Building lastest firedrake version"
	curl -O https://raw.githubusercontent.com/firedrakeproject/firedrake/master/scripts/firedrake-install
	python3 firedrake-install --no-package-manager


$(DIAMOND_BIN): firedrake spud $(DXDIFF_BIN)
	@echo "Installing diamond into venv..."
	{ \
		set -e; \
		source $(FIREDRAKE_ACTIVATION); \
		cd spud/diamond; \
		pip install .; \
	}

$(DXDIFF_BIN): spud firedrake
	@echo "Installing dxdiff into venv..."
	{ \
		set -e; \
		source $(FIREDRAKE_ACTIVATION); \
		cd spud/dxdiff; \
		pip install .; \
	}
