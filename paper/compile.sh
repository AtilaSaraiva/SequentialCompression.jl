#!/usr/bin/env bash

ls | entr -s "latexmk -bibtex -pdf paper.tex" &
