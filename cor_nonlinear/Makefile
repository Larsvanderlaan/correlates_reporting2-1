## all           : clean saved results and get results to populate report
all: clean cor_nonlinear

## clean         : delete saved bootstrap results
clean:
	rm -rf output/$(TRIAL)/*

## cor_nonlinear : create output figures
cor_nonlinear: 
ifeq ($(TRIAL),$(filter $(TRIAL), moderna_mock moderna_real))
	Rscript code/cor_nonlinear.R D57
	Rscript code/cor_nonlinear.R D29
else ifeq ($(TRIAL),$(filter $(TRIAL), janssen_pooled_mock janssen_na_mock janssen_la_mock janssen_sa_mock))
	Rscript code/cor_nonlinear.R D29
	Rscript code/cor_nonlinear.R D29start1
else ifeq ($(TRIAL),$(filter $(TRIAL), janssen_pooled_realbAb janssen_pooled_realPsV janssen_pooled_realADCP janssen_na_realbAb janssen_na_realPsV janssen_na_realADCP janssen_la_realbAb janssen_la_realPsV janssen_la_realADCP janssen_sa_realbAb janssen_sa_realPsV janssen_sa_realADCP))
	Rscript code/cor_nonlinear.R D29IncludeNotMolecConfirmed
	Rscript code/cor_nonlinear.R D29IncludeNotMolecConfirmedstart1
else ifeq ($(TRIAL),$(filter $(TRIAL), hvtn705))
	Rscript code/cor_nonlinear.R D210
else ifeq ($(TRIAL),$(filter $(TRIAL), prevent19))
	Rscript code/cor_nonlinear.R D35
endif

before_script:
	echo "numPerm <- 50" >> code/params.R

## report        : create a standalone report
report:  
	Rscript -e "rmarkdown::render('cor_nonlinear_report.Rmd', output_file='cor_nonlinear_report_$(TRIAL).pdf')"

test: before_script cor_nonlinear report


# type 'make help' to show all make commands
help: Makefile
	@sed -n 's/^##//p' $<

# define PHONYS
.PHONY: help clean cor_nonlinear
