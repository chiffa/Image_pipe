#!/usr/bin/env bash

rm -rf build
rm -rf dist
python setup.py check -r -s
python setup.py bdist_wheel
twine upload -s dist/*
