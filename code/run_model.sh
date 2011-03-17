#!/usr/bin/env bash

rm -fr output/*
octave driver.m
cd output
for i in $( ls *.tex ); do
    perl -pi -e 's/output\//eps\//' ${i}
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
open output/pdf/*.pdf