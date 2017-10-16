#!/usr/bin/env bash

./pre-commit-hook.sh
rm -rf build
rm -rf dist
python setup.py check -r -s
python setup.py bdist_wheel
twine upload -s dist/*
