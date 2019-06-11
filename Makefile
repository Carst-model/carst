.PHONY: test lint

test:
	@echo "    Running simulation"
	@( \
		source ./firedrake/bin/activate; \
		python3 carst/carst_it1.py; \
	)
	@echo "	   Running visualisation"
	@python3 scripts/stratiMesh.py; \

lint:
	@echo "    Linting carst codebase"
	@flake8 carst
	@echo "    Linting carst test suite"
	@flake8 test
	@echo "    Linting carst examples"
	@flake8 examples
	@echo "    Linting carst scripts"
	@flake8 scripts

