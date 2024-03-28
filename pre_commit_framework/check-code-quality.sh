#!/bin/bash 

set -e  # If there is any error it will stop the whole bash script file :) 


EXIT_STATUS=0

black --config .black.toml .  || ((EXIT_STATUS++))
pylint --rcfile .pylintrc *.py  || ((EXIT_STATUS++))
flake8 --config .flake8.toml || ((EXIT_STATUS++))
mypy .  || ((EXIT_STATUS++))
isort . --settings .isort.cfg || ((EXIT_STATUS++))
# try to use ruff also If you want :) 

echo existing with status $EXIT_STATUS 
exit $EXIT_STATUS