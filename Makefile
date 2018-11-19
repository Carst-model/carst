lint:
	@echo "    Linting carst codebase"
	@flake8 carst
	@echo "    Linting carst test suite"
	@flake8 test
	@echo "    Linting carst examples"
	@flake8 examples
	@echo "    Linting carst scripts"
	@flake8 scripts

