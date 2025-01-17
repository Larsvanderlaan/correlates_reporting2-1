#Sys.setenv(TRIAL = "janssen_pooled_realPsV"); COR="D29IncludeNotMolecConfirmedstart1"; Sys.setenv(VERBOSE = 1) 
#Sys.setenv(TRIAL = "prevent19"); COR="D35"; Sys.setenv(VERBOSE = 1)
renv::activate(project = here::here(".."))     
    # There is a bug on Windows that prevents renv from working properly. The following code provides a workaround:
    if (.Platform$OS.type == "windows") .libPaths(c(paste0(Sys.getenv ("R_HOME"), "/library"), .libPaths()))    
source(here::here("..", "_common.R"))

library(kyotil) # p.adj.perm, getFormattedSummary
library(marginalizedRisk)
library(tools) # toTitleCase
library(survey)
library(plotrix) # weighted.hist
library(parallel)
library(forestplot)
library(Hmisc) # wtd.quantile, cut2
library(xtable) # this is a dependency of kyotil
source(here::here("code", "params.R"))
time.start=Sys.time()
myprint(study_name)
myprint(verbose)

all.markers=paste0("Day", tpeak, assays)

# path for figures and tables etc
save.results.to = here::here("output");                                if (!dir.exists(save.results.to))  dir.create(save.results.to)
save.results.to = paste0(save.results.to, "/", attr(config,"config")); if (!dir.exists(save.results.to))  dir.create(save.results.to)
save.results.to = paste0(save.results.to, "/", COR,"/");               if (!dir.exists(save.results.to))  dir.create(save.results.to)
print(paste0("save.results.to equals ", save.results.to))

# some exploratory code
if (config$is_ows_trial) source(here::here("code", "cor_coxph_misc.R"))


# B=1e3 and numPerm=1e4 take 10 min to run with 30 CPUS for one analysis
B <-       config$num_boot_replicates 
numPerm <- config$num_perm_replicates # number permutation replicates 1e4
myprint(B)
myprint(numPerm)

# uloq censoring, done here b/c should not be done for immunogenicity reports
# note that if delta are used, delta needs to be recomputed
for (a in assays) {
  for (t in "Day"%.%tpeak ) {
    dat.mock[[t %.% a]] <- ifelse(dat.mock[[t %.% a]] > log10(uloqs[a]), log10(uloqs[a]), dat.mock[[t %.% a]])
  }
}    

dat.vac.seroneg=subset(dat.mock, Trt==1 & ph1)
dat.pla.seroneg=subset(dat.mock, Trt==0 & ph1)

# define an alias for EventIndPrimaryDxx
dat.vac.seroneg$yy=dat.vac.seroneg[[config.cor$EventIndPrimary]]
dat.pla.seroneg$yy=dat.pla.seroneg[[config.cor$EventIndPrimary]]

myprint(tfinal.tpeak)
write(tfinal.tpeak, file=paste0(save.results.to, "timepoints_cum_risk_"%.%study_name))
    
# define trichotomized markers
dat.vac.seroneg = add.trichotomized.markers (dat.vac.seroneg, tpeak, wt.col.name="wt")
marker.cutpoints=attr(dat.vac.seroneg, "marker.cutpoints")
for (a in assays) {        
    for (t in "Day"%.%tpeak) {
        q.a=marker.cutpoints[[a]][[t]]
        write(paste0(labels.axis[1,a], " [", concatList(round(q.a, 2), ", "), ")%"), file=paste0(save.results.to, "cutpoints_", t, a, "_"%.%study_name))
    }
}
    
# create verification object to be populated by the following scripts
rv=list() 
rv$marker.cutpoints=marker.cutpoints



###################################################################################################
# run PH models
###################################################################################################

#create twophase design object
design.vacc.seroneg<-twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~ph2, data=dat.vac.seroneg)

source(here::here("code", "cor_coxph_ph.R"))

# optional forest plots
if(length(config$forestplot_script)==1 & study_name!="PREVENT19") {
    tmp=here::here("code", config$forestplot_script)
    if (file.exists(tmp)) source(tmp)
    
    # unit testing 
    if (study_name == "MockCOVE") {
        tmp.1=c(sapply(rv$fr.2[-1], function (x) x[c("HR","p.value"),1])) # concatList(tmp.1, ", ")
        if (tpeak=="29") {
            tmp.2=c(2.19803e-01,3.42813e-06,4.00791e-01,1.55780e-03,2.64497e-01,2.90077e-04,2.52391e-01,3.38292e-04,3.68861e-01,6.24978e-03)
        } else if (tpeak=="57") {
            tmp.2=c(1.17284e-01,4.73761e-11,3.91017e-01,7.49144e-04,2.84943e-01,1.36601e-05,2.44480e-01,9.03454e-06,2.77729e-01,8.84152e-06)
        }
        assertthat::assert_that(
            max(abs(tmp.1-tmp.2)/abs(tmp.2))<1e-5,
            msg = "failed sanity check")    
        print("Passed sanity check")    
    }
}



###################################################################################################
# draw marginalized risk curves
###################################################################################################
    
# load ylims.cor[[1]] from D29 analyses, which is a list of two: 1 with placebo lines, 2 without placebo lines.
tmp=paste0(here::here(), paste0("/output/", attr(config,"config"), "/", COR, "/ylims.cor.", study_name, ".Rdata"))
if (file.exists(tmp)) load(tmp)
# if this does not exist, the code will find alternative ylim

source(here::here("code", "cor_coxph_marginalized_risk_no_marker.R"))
source(here::here("code", "cor_coxph_marginalized_risk_bootstrap.R"))
source(here::here("code", "cor_coxph_marginalized_risk_plotting.R"))



###################################################################################################
# save verification object rv
save(rv, file=paste0(here::here("verification"), "/", COR, ".rv."%.%study_name%.%".Rdata"))

print("cor_coxph run time: "%.%format(Sys.time()-time.start, digits=1))
