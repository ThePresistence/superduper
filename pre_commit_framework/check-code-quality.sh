#!/bin/bash 

set -e  # If there is any error it will stop the whole bash script file :) 

black --config .black.toml . 
pylint --rcfile .pylintrc *.py 
flake8 --config .flake8.toml
mypy . 
isort . --settings .isort.cfg
# try to use ruff also If you want :) 
