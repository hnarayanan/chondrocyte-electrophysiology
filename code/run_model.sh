#!/usr/bin/env bash

rm -fr results/*
octave driver.m
cd results
for i in $( ls *.tex ); do
    perl -pi -e 's/results\//eps\//' ${i}
done
mkdir eps
mv *.eps eps
for i in $( ls *.tex ); do
    pdflatex ${i}
done
rm *.log *.aux
mkdir tex
mv *.tex tex
mkdir pdf
mv *.pdf pdf
cd ..
open results/pdf/*.pdf