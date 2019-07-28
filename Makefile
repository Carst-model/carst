.PHONY: test lint clean diamond

# Paths for setup
PROJECT_DIR := $(shell pwd -P)
FIREDRAKE_VENV_FULL := $(PROJECT_DIR)/firedrake

# Paths to alias long targets
DIAMOND_PIP := firedrake/lib/python3.7/site-packages/diamond
DXDIFF_PIP := firedrake/lib/python3.7/site-packages/dxdiff
SPUD_BASE := firedrake/share/spud/spud_base.rnc

# For sourcing the shell file to activate the firedrake venv
FIREDRAKE_ACTIVATION := firedrake/bin/activate

# Temporary test
all: diamond_default.rng $(DIAMOND_PIP) firedrake
	{ \
		set -e; \
		source $(FIREDRAKE_ACTIVATION); \
		diamond -s diamond_default.rng; \
	}

test: firedrake carst basic_tests.py scripts
	@echo "Running simulation"
	{ \
		set -e; \
		source $(FIREDRAKE_ACTIVATION); \
		python3 basic_tests.py; \
	}
	@echo "Running visualisation"
	{ \
		set -e; \
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
	-rm -rf firedrake firedrake-install

diamond_default.rng: $(DIAMOND_PIP) diamond_default.rnc
	@echo "Building schema..."
	{ \
		set -e; \
		source $(FIREDRAKE_ACTIVATION); \
		./firedrake/usr/bin/spud-preprocess diamond_default.rnc; \
	}

planning_diagram.png:
	@echo "Building diagram..."
	ditaa planning_diagram.ditaa

firedrake:
	@echo "Building lastest firedrake version"
	curl -O https://raw.githubusercontent.com/firedrakeproject/firedrake/master/scripts/firedrake-install
	python3 firedrake-install --no-package-manager

# Diamond alias for easy terminal use
diamond: $(DIAMOND_PIP)

$(DIAMOND_PIP): firedrake spud/diamond $(DXDIFF_PIP) $(SPUD_BASE)
	@echo "Installing diamond into venv..."
	{ \
		set -e; \
		source $(FIREDRAKE_ACTIVATION); \
		cd spud/diamond; \
		pip install -r requirements.txt; \
		pip install .; \
	}

$(SPUD_BASE): spud/schema firedrake
	@echo "Installing base spud..."
	mkdir -p firedrake/share/spud
	cp spud/schema/spud_base.* firedrake/share/spud

$(DXDIFF_PIP): spud/dxdiff firedrake
	@echo "Installing dxdiff into venv..."
	{ \
		set -e; \
		source $(FIREDRAKE_ACTIVATION); \
		cd spud/dxdiff; \
		pip install .; \
	}
