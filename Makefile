.PHONY: test lint clean plan

firedrake:
	@echo "Building lastest firedrake version"
	curl -O https://raw.githubusercontent.com/firedrakeproject/firedrake/master/scripts/firedrake-install
	python3 firedrake-install --no-package-manager

test: firedrake
	@echo "Running simulation"
	@( \
		source ./firedrake/bin/activate; \
		python3 carst/carst_it1.py; \
	)
	@echo "Running visualisation"
	@python3 scripts/stratiMesh.py; \

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
	@echo "Removing firedrake..."
	rm -r firedrake firedrake-install
	@echo "Removing versioneer manifest..."
	rm MANIFEST.in

plan:
	@echo "Removing old diagram..."
	rm planning_diagram.png
	@echo "Building diagram..."
	ditaa planning_diagram.ditaa
