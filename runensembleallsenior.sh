#!/bin/bash


if [ $# -eq 0 ]
       then
     echo "********* Please provide a module name, e.g. cor_coxph, as an argument."
     exit
     fi
     
     
     
     for T in  janssen_pooled_partAsenior janssen_na_partAsenior janssen_la_partAsenior  janssen_pooled_partAnonsenior janssen_na_partAnonsenior janssen_la_partAnonsenior janssen_sa_partAnonsenior
     do
     export TRIAL=$T 
     echo $TRIAL
     if [[ "$1" == "cor_report" ]] 
     then
     make cor_report
     else
       make -k -C $1 all
        Rscript -e "TRIAL <- Sys.getenv('TRIAL'); print(TRIAL); bookdown::render_book(input = 'cor_threshold/report.Rmd', output_file = 'covpn_correlates_$1_$TRIAL.pdf', config_file = '_bookdown_$1.yml', output_format = bookdown::pdf_document2(toc_depth=3), quiet=TRUE)"
     fi
     
     done