#! /bin/bash

rm ../*.html

find . -name \*.md -type f -exec pandoc -t html -s --latexmathml -o ../{}.html {} --template template.html --css stylesheets/stylesheet.css -A footer.html \; && rename 's#\.md\.html#\.html#' ../*
