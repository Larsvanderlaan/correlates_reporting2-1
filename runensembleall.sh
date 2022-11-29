#!/bin/bash


if [ $# -eq 0 ]
  then
    echo "********* Please provide a module name, e.g. cor_coxph, as an argument."
    exit
fi



for T in  janssen_sa_partA janssen_pooled_partA janssen_na_partA janssen_la_partA  
do
    export TRIAL=$T 
    echo $TRIAL
    if [[ "$1" == "cor_report" ]] 
    then
        make cor_report
    else
        make -k -C $1 all
        Rscript -e "TRIAL <- Sys.getenv('TRIAL'); print(TRIAL); bookdown::render_book(input = 'cor_threshold/report.Rmd', output_file = 'covpn_correlates_$1_$TRIAL.pdf', config_file = '_bookdown_$1.yml', output_format = bookdown::pdf_document2(toc_depth=3), quiet=TRUE)"

       # bash ./_build_chapter.sh $1
    fi

done
