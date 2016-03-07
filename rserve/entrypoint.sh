#!/bin/bash

cd /tmp
rm -rf /tmp/vcf*
R -e "package.skeleton(name='vcf2ld', code_files='/tmp/vcf2ld.R')"

sed -i 's/%%  ~~function to do ... ~~/  ~~function to do ... ~~/' vcf2ld/man/*
sed -i 's|%%  ~~data name/kind ... ~~|  ~~data name/kind ... ~~|' vcf2ld/man/*
sed -i 's|%%   ~~ data name/kind ... ~~|  ~~data name/kind ... ~~|' vcf2ld/man/*
sed -i 's|~~ Optionally other standard keywords, one per line, from file KEYWORDS in ~~||' vcf2ld/man/*
sed -i 's|~~ the R documentation directory ~~||' vcf2ld/man/*

R CMD build vcf2ld
R CMD INSTALL vcf2ld_1.0.tar.gz
R CMD Rserve.dbg --RS-conf /usr/share/rserve/config/rserve.conf --vanilla --no-save
