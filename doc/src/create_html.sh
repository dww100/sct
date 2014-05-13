#! /bin/bash

rm ../*.html
find . -name \*.md -type f -exec pandoc -t html -s --latexmathml -o ../{}.html {} \; && rename 's#\.md\.html#\.html#' ../*
