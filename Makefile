all: lint test integration

lint:
	# stop the build if there are Python syntax errors or undefined names
	flake8 . --count --select=E9,F63,F7,F82 --show-source --statistics
	# exit-zero treats all errors as warnings.
	flake8 . --count --exit-zero --max-complexity=20 --statistics

test:
	pytest

integration:
	pytest --override-ini='python_files=integration_*.py' --override-ini='python_functions=integration_*'
