#!/usr/bin/env bash

sphinx-apidoc -fo doc/source imagepipe
sphinx-build -b html doc/source /build >> doc_build.log 2>>doc_build_err.log
pandoc -o README.md inlined_readme.rst
