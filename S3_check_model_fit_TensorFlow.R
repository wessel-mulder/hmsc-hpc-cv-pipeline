#################################################
# Here we fit the models defined for the system
################################################
cat("Script started\n"); flush.console()

remove(list=ls())
set.seed(369)
library(jsonify,lib='~/Rlibs')
library(colorspace,lib='~/Rlibs')

library(RColorBrewer,lib='~/Rlibs')
library(farver,lib='~/Rlibs')
library(scales,lib='~/Rlibs')
library(Hmsc,lib='~/Rlibs')
library(cli,lib='~/Rlibs')
library(vioplot,lib="~/Rlibs")
cat("Libraries loaded\n"); flush.console()

### Set up directories #### Because I run this on two difference computers this
#setwd(dirname(rstudioapi::getSourceEditorContext()$path))

### Set up directories ####
guilds <- c("Woodpeckers","Plovers")
env_vars <- c('LandusePercs','Climate')
for(guild in guilds){
for(env_var in env_vars){

# assigne properly
model_description = sprintf("2026-01-20_12-40-41_%s_%s_Atlas3",guild,env_var)

localDir = sprintf("./HmscOutputs/%s",model_description)
ModelDir = file.path(localDir, "Models/Fitted")
ResultDir = file.path(localDir, "Results")

showBeta = TRUE
showGamma = TRUE
showOmega = TRUE
showAlpha = TRUE

maxOmega = 100

ma.beta = NULL
lables.beta = NULL
ma.gamma = NULL
lables.gamma = NULL
ma.omega= NULL
lables.omega=NULL

samples_list = c(250)
thin_list = c(10)
transient = 100000
nChains = 4
nParallel = 4
Lst = 1
verbose = 1

#Add back in the loop for multiple model lengths in the future
text.file = file.path(ResultDir,paste0(model_description,"model_fit_details.txt"))
cat(sprintf("%s\nThis file contains human readable information regarding the model fit.\nFitting Date:\t%s\n%1$s\n\n",
            strrep("-",80), strptime(Sys.time(),format = "%Y-%m-%d %H:%M")),
    file=text.file)

while(Lst <= length(samples_list)){
  samples = samples_list[Lst]
  thin = thin_list[Lst]
  transient = transient
  filename = file.path(ModelDir,sprintf("HPC_samples_%.4d_thin_%.2d_chains_%.1d.Rdata", samples, thin, nChains))
  if(file.exists(filename)){
    ptm = proc.time()
    load(file = filename)
    posteriors = fitted_model$posteriors
    fit_times = fitted_model$fit_times
    labels = c("Names","Thin", "Samples", "Chains")
    labels2 = c("Min.   :","1st Qu.:","Median :","Mean   :","3rd Qu.:","Max.   :")
    values = c(model_description,sprintf("%.3d",thin),sprintf("%.3d",samples),sprintf("%.2d",nChains))
    rm(fitted_model)
    cli_h1("Model running")
    #cli_text("Model Info:\n Samples:{samples} Thin:{thin}")
    #cli_text("Current time: {format(Sys.time())}")
    #cli_rule()
    cat(sprintf("%s\n%-35s\t\t\t|\t%-s",strrep("-",80),"Model Description","Fitting Times"),file=text.file,append = TRUE)
    cat(sprintf("\n%-05s\t\t%-30s\t|\tChain %.1d\t\t%.2f s", labels, values, c(1,2,3,4), unlist(fit_times)),file=text.file,append = TRUE)
    cat(sprintf("\n%35s\t\t\t|\t\t\t\t-------\n%1$35s\t\t\t|\tTotal\t\t%.2f s\n","",sum(unlist(fit_times))),file=text.file,append = TRUE)
    #Note that this code can only one one model at a time.
    nm = 1
    mpost = convertToCodaObject(posteriors, spNamesNumbers = c(T,F), covNamesNumbers = c(T,F))
    nr = posteriors$nr
    #plot(mpost$Beta)
    #We use the showBeta, showGamma extra to allow the user to define if they want these values shown
    #I think this loop can be improved
    cat(sprintf("%s\n%50s\n",strrep("-",80),"Parameter information"),file=text.file,append = TRUE)
    #Empty predefined dataframe for Beta and Gamma parameters
    Beta_Gamma = data.frame(Beta_PE = rep(NA,6), Beta_ES = NA, Gamma_PE=NA)
    if(showBeta){
      cli_progress_step("Calculating Beta")
      #Add in effective sampling size to this part...
      eff = effectiveSize(mpost$Beta)
      tmp_eff = summary(eff)
      psrf = gelman.diag(mpost$Beta,multivariate=FALSE)$psrf
      #tmp = summary(psrf)
      Beta_Gamma$Beta_PE = summary(psrf[,1])
      Beta_Gamma$Beta_ES = summary(eff)
      #cat(sprintf("---------------\nBeta\n\nPoint est:\t\t\t\tEffective Size:\n\t%s\t\t%.2f\n\t%s\t\t%.2f\n\t%s\t\t%.2f\n\t%s\t\t%.2f\n\t%s\t\t%.2f\n\t%s\t\t%.2f\n",tmp[1],tmp_eff[[1]],tmp[2],tmp_eff[[2]],tmp[3],tmp_eff[[3]],tmp[4],tmp_eff[[4]],tmp[5],tmp_eff[[5]],tmp[6],tmp_eff[[6]]),file=text.file,append = TRUE)
      if(is.null(ma.beta)){
        ma.beta = psrf[,1]
        lables.beta = paste0(as.character(thin),",",as.character(samples))
      } else {
        ma.beta = cbind(ma.beta,psrf[,1])
        lables.beta = c(lables.beta,paste0(as.character(thin),",",as.character(samples)))
      }
      cli_progress_done()
    }
    if(showGamma){
      cli_progress_step("Calculating Gamma")
      psrf = gelman.diag(mpost$Gamma,multivariate=FALSE)$psrf
      #tmp = summary(psrf)
      Beta_Gamma$Gamma_PE = summary(psrf[,1])
      #cat(sprintf("---------------\nGamma\n\nPoint est:\n\t%s\n\t%s\n\t%s\n\t%s\n\t%s\n\t%s\n",tmp[1],tmp[2],tmp[3],tmp[4],tmp[5],tmp[6]),file=text.file,append = TRUE)
      if(is.null(ma.gamma)){
        ma.gamma = psrf[,1]
        lables.gamma = paste0(as.character(thin),",",as.character(samples))
      } else {
        ma.gamma = cbind(ma.gamma,psrf[,1])
        lables.gamma = c(lables.gamma,paste0(as.character(thin),",",as.character(samples)))
      }
      cli_progress_done()
    }
    if(showBeta|showGamma){
      cat(sprintf("%s\n%-35s\t\t\t|\t%-s\n%35s\t\t\t|\n%s\t\t\t%-s\t\t\t|\t%5$s",strrep("-",80),"Beta","Gamma","","Point est:","Effective Size:"),file=text.file,append = TRUE)
      cat(sprintf("\n %s\t%.4f\t\t\t\t%.2f\t\t|\t\t%.4f",labels2,Beta_Gamma$Beta_PE,Beta_Gamma$Beta_ES,Beta_Gamma$Gamma_PE),file=text.file,append = TRUE)
    }
    if(showOmega & nr>0){
      cli_progress_step("Calculating Omega")
      #Omega = matrix(NA,nrow=6,ncol=nr)
      labels2 = c(" Min.   :"," 1st Qu.:"," Median :"," Mean   :"," 3rd Qu.:"," Max.   :")
      line = c("","","","","","","","","")
      for(k in 1:nr){
        tmp = mpost$Omega[[k]]
        z = dim(tmp[[1]])[2]
        #Randomly sample maxOmega species pairs, since all pairs would take too long
        if(z > maxOmega){
          sel = sample(1:z, size = maxOmega)
          for(i in 1:length(tmp)){
            tmp[[i]] = tmp[[i]][,sel]
          }
        }
        psrf = gelman.diag(tmp, multivariate = FALSE)$psrf
        #tmp = summary(psrf)
        tmp = summary(psrf[,1])
        #This line based text output is because there can be any number of random variables and each on is reported next to the last, so the lines grow in length each iteration
        line[1] = paste0(line[1],sprintf("Factor %i\t\t\t|%3s",k,""))
        line[2] = paste0(line[2],sprintf("%3s\t\t\t\t\t|",""))
        line[3] = paste0(line[3],sprintf("Point est:\t\t\t|%3s",""))
        for(x in 1:6){
          line[x+3] = paste0(line[x+3], sprintf("%s%3s%.2f\t|%2$3s",labels2[x],"",tmp[x]))
        }
        if(is.null(ma.omega)){
          ma.omega = psrf[,1]
          lables.omega = paste0(as.character(thin),",",as.character(samples))
        } else {
          ma.omega = cbind(ma.omega,psrf[,1])
          lables.omega = c(lables.omega,paste0(as.character(thin),",",as.character(samples)))
        }
      }
      cli_progress_done()
      cat(sprintf("\n%s\nOmega\n",strrep("-",80)),file=text.file,append = TRUE)
      for(x in 1:9){
        cat(line[x],"\n",file=text.file,append = TRUE)
      }
    }
    if(!is.null(mpost$Rho)){
      cli_progress_step("Calculating Rho")
      psrf = gelman.diag(mpost$Rho,multivariate=FALSE)$psrf
      cat(sprintf("%s\nRho\n\nPoint est:\nMean:\t%.2f\n",strrep("-",80),psrf[1]),file=text.file,append=TRUE)
    }
    
    if(showAlpha & nr>0){
      cli_progress_step("Calculating Alpha")
      line = c("","","","","","","","","")
      for(k in 1:nr){
        if(posteriors$ranLevels[[k]]$sDim>0){
          psrf = gelman.diag(mpost$Alpha[[k]],multivariate = FALSE)$psrf
          #cat(sprintf("\nRandom factor:\t%s\n\nPoint est:",names(posteriors$ranLevels)[k]),file=text.file,append=TRUE)
          #cat(sprintf("\n%s:\t%.2f",names(psrf[,1]),psrf[,1]),file=text.file,append=TRUE)
          #cat("\n",file=text.file,append=TRUE)
          line[1] = paste0(line[1],sprintf("Random factor: %10s\t\t\t|%1s",names(posteriors$ranLevels)[k],""))
          line[2] = paste0(line[2],sprintf("%25s\t\t\t|",""))
          line[3] = paste0(line[3],sprintf("%-25s\t\t\t|","Point est:"))
          for(x in 1:length(psrf[,1])){
            line[x+3] = paste0(line[x+3], sprintf("%s%s%1$3s%.2f%1$3s\t\t\t|"," ",names(psrf[,1])[x],psrf[x,1]))
          }
        }
      }
      cat(sprintf("%s\nalpha\n%1$s\n",strrep("-",80)),file=text.file,append=TRUE)
      for(x in 1:length(psrf[,1])+3){
        cat(line[x],"\n",file=text.file,append = TRUE)
      }
      cli_progress_done()
    }
    computational.time = proc.time() - ptm
    cat("Time taken:", computational.time[3],"s \nCurrent time:",format(Sys.time(), "%H:%M:%S"),"\n\n")
  } else {
    cat(sprintf("File not found for:\nModel Description\nName:\t\t%s\nThin:\t\t%.3d\nSamples:\t%.4d\nChains:\t\t%.2d\n\n",model_description,thin,samples,nChains),file=text.file,append = TRUE)
  }
  
  cat(sprintf("%s\n%50s\n%1$s\n",strrep("=",80),"END OF CURRENT MODEL"),file=text.file,append=TRUE)
  Lst = Lst + 1
}
pdf(file= file.path(ResultDir,paste0("/",model_description,"MCMC_convergence.pdf")))
par(mfrow=c(3,2))
vioplot(ma.beta,col=rainbow_hcl(nm),names=lables.beta,ylim=c(0,max(ma.beta)),main="psrf(beta)")
#legend("topright",legend = names(models), fill=rainbow_hcl(nm))
vioplot(ma.beta,col=rainbow_hcl(nm),names=lables.beta,ylim=c(0.9,1.1),main="psrf(beta)")

vioplot(ma.gamma,col=rainbow_hcl(nm),names=lables.gamma,ylim=c(0,max(ma.gamma)),main="psrf(gamma)")
#legend("topright",legend = names(models), fill=rainbow_hcl(nm))
vioplot(ma.gamma,col=rainbow_hcl(nm),names=lables.gamma,ylim=c(0.9,1.1),main="psrf(gamma)")

vioplot(ma.omega,col=rainbow_hcl(nm),names=lables.omega,ylim=c(0,max(ma.omega)),main="psrf(omega)")
#legend("topright",legend = names(models), fill=rainbow_hcl(nm))
vioplot(ma.omega,col=rainbow_hcl(nm),names=lables.omega,ylim=c(0.9,1.1),main="psrf(omega)")

}
}