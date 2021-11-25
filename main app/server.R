#R.home()
#/Library/Frameworks/R.framework/Resources/Rscript /Users/xid32/desktop/webserver/DEApp-master/app.R
library(limma)
library(gtools)
library(vsn)
library(naturalsort)
library("shinyWidgets")
library(data.table)
library(plyr)
library("Hmisc")
library(gplots)
library(shinyHeatmaply)
library(matrixStats)
library(shinyShortcut)
library(heatmaply)
library(ber) 
library(shiny)
library(plotly)
library(viridis)
library(jsonlite)
library(RColorBrewer)
library(readxl)
library(DT)
library(xtable)
library(htmltools)
library(htmlwidgets)
library(dplyr)
library(tidyverse)
library(hrbrthemes)
library(viridis)
library(DelayedMatrixStats)
library("xlsx")

Neighbour=function(slidedata,proteinType){ 
  attach(slidedata) 
  No_Block=max(slidedata[,1])
  No_col=max(slidedata[,2]) 
  No_row=max(slidedata[,3]) 
  if(proteinType==532) {
    data=data.frame(Block,Column,Row,B532.Median) 
  } else {
    data=data.frame(Block,Column,Row,B635.Median) 
  }
  Bneighbour=c() 
  for(i in 1:No_Block) {
    Bmatrix=matrix(data[data[,1]==i,4],nrow=No_row,ncol=No_col,byrow=TRUE) 
    Nmatrix=ma3x3.matrix(Bmatrix,Fun=median,na.rm=TRUE) 
    block_vector=as.vector( t(Nmatrix)) 
    Bneighbour=c(Bneighbour,block_vector) 
  }
  return (Bneighbour) 
}


BottomPositive=function( aggr ) {
  library(limma)
  library(naturalsort) 
  elist=list() 
  index=seq(1:nrow(aggr)) 
  aggrx=cbind(aggr,index) 
  pcs=aggrx[grep("Gliadin",rownames(aggrx)),] 
  pcs=as.matrix(pcs) 
  if(nrow(pcs)>0 ) {
    pcsNo=as.numeric(pcs[nrow(pcs),]) 
    aggrx=aggrx[-pcsNo,]  
  }
  pc1=aggrx[grep("ontrol",rownames(aggrx)),] 
  pc2_1=aggrx[grep("Anti-Ig",rownames(aggrx)),]
  pc2_2=aggrx[grep("anti-Ig",rownames(aggrx)),]
  pc2=rbind(pc2_1,pc2_2) 
  
  pc3=aggrx[grep("IgG",rownames(aggrx)),]
  pc4=aggrx[grep("IgM",rownames(aggrx)),]
  pc5=aggrx[grep("IgA",rownames(aggrx)),]
  pc6=aggrx[grep("IgE",rownames(aggrx)),]
  pc=do.call("rbind",list(pc1,pc2,pc3,pc4,pc5,pc6)) 
  positiveNum=as.numeric(pc[,ncol(pc)]) 
  
  No_positive=length(positiveNum) 
  
  if(No_positive > 1 ) { 
    nopositive=aggr[- positiveNum,] 
    nopositive=avereps(nopositive) 
    nopositive=nopositive[naturalorder(rownames(nopositive)),]
    positive=aggr[positiveNum,] 
    positive=avereps(positive) 
    positive=positive[naturalorder(rownames(positive)),]
    elist$No_positive=nrow(positive) 
    elist$newdata=rbind(nopositive,positive) 
  } else {
   elist$newdata=avereps(aggr)  
   elist$No_positive=No_positive 
    
  }
  
  pc2=as.matrix(pc2) 
     
 if(nrow(pc2) ==0) { 
     elist$No_anti_Ig <- 0 } else { 
         if(ncol(pc2)==1) { elist$No_anti_Ig <- 1 } else {
         elist$No_anti_Ig <- nrow(pc2)
       }
   }
  elist$No_Ig_control= elist$No_positive -elist$No_anti_Ig 
   
  return (elist) 
}


## average PBS 
avePBS=function(pbs) {
  pbs=as.matrix(pbs) 
  pbs_ave=c() 
  for(i in 1:nrow(pbs)) { 
    antirow=sort(na.omit(pbs[i,]))
    anti_len=length(antirow) 
    anti_cut=ceiling( length(antirow)*0.1) 
    L=anti_len-anti_cut 
    if( anti_len < 3) {
      pbs_ave[i]=mean(antirow) } else {
        pbs_ave[i]=mean(antirow[1:L]) 
      }
    if(is.na(pbs_ave[i])) pbs_ave[i]=1e-3 }
  return (pbs_ave) 
}

# No negative value in the matrix 
No_minus=function(A) {
  A[A<0]<- 0.001 
  return (A) 
}

# Remove batch effect with user's instruction
Batch_correction_Method=function(Mtx, target, method){ 
  if ('Batch' %in% names(target)){
    Batch = as.factor(target[,4])
    class = NULL
    if ('Group' %in% names(target)){
      class = as.factor(target[,5])
      class<-data.frame(class)
    }
      
    if(method==1){     #ANOVA
      Mtx = t(ber_bg(t(Mtx), batch, class))
    }
    else if(method==2){ #combat_p
      Mtx = t(combat_p(t(Mtx), batch, class))
    }
    
    else if(method==3){ #combat_np
      Mtx = t(combat_np(t(Mtx), batch, class))
    } 
  }
  return(Mtx)
} 

 
# normalization process with user's intruction 

Normalization_Method=function(Intensity,background,target,method) {
         
       # calculate background correction factor 
       for(i in 1:ncol(background)) {
	 		       Lower=quantile(background[,i],probs=0.05) 
	 		       Upper=quantile(background[,i],probs=0.95) 
	 		        
	 		       background[,i]=ifelse(background[,i] > Upper, Upper,background[,i]) 
	 		       background[,i]=ifelse(background[,i] < Lower, Lower,background[,i]) 
	 }
	
	BGcolAve=colMeans(background) 
	BGtotalmean=mean(BGcolAve) 
	BCF=BGtotalmean/BGcolAve 
	
        ## delete PBS and bad samples 
        Num=seq(1:nrow(target)) 
	target=cbind(target,Num) 
	Bad=target[target[,3]=="Bad",] 
	PBS=target[target[,3]=="PBS",] 
	PBS_Bad=rbind(PBS,Bad) 
	delNo=PBS_Bad[,ncol(target)] 
        target=target[-delNo,] 
        ##***********************************
       
        target[,1:2]=apply(target[,1:2],2,as.numeric) 
       
       
       
         
        if(method==1){                         #RLM 
             n.slides <-   max(target[,1])      
             n.blocks <-   max(target[,2])
             n.sample <-   ncol(Intensity)  
      
             pc1=Intensity[grep("ontrol",rownames(Intensity)),] 
             pc2=Intensity[grep("nti-Ig",rownames(Intensity)),] 
             positive=rbind(pc1,pc2) 
             n.controls=nrow(positive) 
             
             
              
             
             if( length(positive) ==  0 )  
             { 
                 Ntx=Intensity
             } else {
                         
               positive=avereps(positive) 
               positive=positive[naturalorder(rownames(positive)),]
               n.antigen=nrow(Intensity) 
               if( n.controls > 2) {
	                    positive =positive[- grep("nti-Ig",rownames(positive)),] 
	                    n.controls=nrow(positive) 
               }
               
               control=positive
      
               if( is.data.frame(control)) {
                 x=as.matrix(log2(control+1),ncol=1)  
               } else  { 
                 x=as.vector(log2(control+1)) 
               }
      
               # Create matrix of explanatory variables for rlm
               # rows: correspond to the rows of x
               dummy <- matrix(0,ncol=n.slides+n.blocks+n.controls,nrow=length(x) ) 
	            
	       for(arrN in 1:n.slides){ 
	         for(controlN in 1:n.controls){
	            for(i in 1:nrow(target)) { 
	                 k=(controlN-1)*n.sample +i 
	                 if(target[i,1]== arrN)  dummy[k,arrN] <- 1 
	             }
	          }
	        }
	             
	            
	        for(blockN in 1:n.blocks){
	           for(controlN in 1:n.controls){
	    	     for(i in 1:nrow(target) ) { 
	    	        k=(controlN-1)*n.sample +i 
	    	        if(target[i,2]== blockN)  dummy[k,n.slides+blockN] <- 1 
	    	        dummy[k,n.slides+n.blocks+controlN] <- 1
	    	               
	             }
	           }
                }
         
                       
               # add 3 rows under the dummies matrix 
               # (likely to denote which columns pertain to arrays, which ones to blocks, 
               # and which ones to controls)
               # with  -1's under first n.slides columns, then under n.blocks.contr columns,
               # then under n.controls columns
             
               dummy <- rbind(dummy, 
                       c(rep(-1,n.slides), rep(0, n.blocks), rep(0,n.controls)),
                       c(rep(0, n.slides), rep(-1,n.blocks), rep(0,n.controls)), 
                       c(rep(0, n.slides), rep(0, n.blocks), rep(-1,n.controls)))
               # add intercept
               dummy <- cbind(rep(1, {length(x)+3}), dummy)
                     
               rlm.result <- rlm(y=c(x,0,0,0), x=dummy, method="M", psi=psi.bisquare)
               signal=log2(Intensity+1)  
             
               # merge all coefficients pertaining to the effect of slides into one vector
     	       coef1vect <- c(sapply(2:(1+n.slides),function(arrNum){ rep(rlm.result$coefficients[arrNum],n.blocks) }))  
               coef1vect <- coef1vect-mean(coef1vect)   
               # merge all coefficients pertaining to the effect of blocks into one vector
               coef2vect <- rep(rlm.result$coefficients[(2+n.slides):(1+n.slides+n.blocks)],n.slides )
               coef2vect <- coef2vect-mean(coef2vect) 
           
               normalized_RLM=sapply(1:ncol(signal),function(x) signal[,x]-coef1vect[x]-coef2vect[x] )   
             
               Ntx=2^(normalized_RLM)  
               
                 }
           
       
        }
         else if(method==2){  ## BCF: Background Correction Factor
	        
	        
	                Ntx=Intensity
	                for(j in 1:ncol(Intensity)) {
	                Ntx[,j]=Intensity[,j]*BCF[j] 
	                }
	                          
	        
	        }
	else if (method==3) { ## scaling
	              pc1=Intensity[grep("ontrol",rownames(Intensity)),] 
                      pc2=Intensity[grep("nti-Ig",rownames(Intensity)),] 
                      positive=rbind(pc1,pc2)   
		      positive=avereps(positive) 
		      positive=positive[naturalorder(rownames(positive)),]
		      n.controls=nrow(positive)
		      
		      if( n.controls > 2) {
		                    control =positive[- grep("nti-Ig",rownames(positive)),] 
		                    n.controls= nrow(control) 
                      } else {
                      
                         control=t( as.matrix(positive[2,])) 
                         n.controls=nrow(control) 
                      
                      }
		    
		    
		    
		    
		      for(i in 1:n.controls) {
		       Lower=quantile(control[i,],probs=0.05) 
		       Upper=quantile(control[i,],probs=0.95) 
		        
		       control[i,]=ifelse(control[i,] > Upper,Upper,control[i,]) 
		       control[i,]=ifelse(control[i,] < Lower, Lower,control[i,]) 
		      }
		      		      
		      colAve=colMeans(control) 
		      overall_mean=mean(colAve)
		      correctFactor=overall_mean/colAve 
		               
		      Ntx=Intensity
		      for(j in 1:ncol(Intensity)) {
		          Ntx[,j]=Intensity[,j]*correctFactor[j] 
	              }
	        
               }
        
        else if(method==4){ #quantile 
          Ntx =normalize.quantiles(Intensity)  
        }
        
        else if(method==5){ #Lowess
          Ntx=normalizeCyclicLoess(Intensity,method="fast") 
        } 
       
       else if(method==6){ # vsn
         vsnnor=normalizeVSN(Intensity) 
         Ntx=2^(vsnnor) 
       } 

       rownames(Ntx)=rownames(Intensity) 
       colnames(Ntx)=colnames(Intensity) 
       
        
    
    
    return (Ntx) 
 
}




Pre_process=function(Elist,measure_method,koc,batch_method,norm_method) {

   elist=list()
   flagdata=Elist$rawdata 
   Num=nrow(flagdata)
   flagdata=cbind(flagdata,Num) 
   
   flagdata[,7]=Elist$rawdata[,6] 
   flagdata[,8]=Elist$rawdata[,7] 
   
   if(measure_method==1) {
      flagdata[,1]=Elist$rawdata[,2] 
      flagdata[,2]=Elist$rawdata[,3] 
      flagdata[,3]=flagdata[,2]-flagdata[,1] 
      flagdata[,4]=flagdata[,2]/flagdata[,1] 
      flagdata[,5]=flagdata[,2]-Elist$rawdata[,1] 
      flagdata[,6]=flagdata[,2]/Elist$rawdata[,1] 
   } else if(measure_method ==2) {
      flagdata[,1]=Elist$rawdata[,4] 
      flagdata[,2]=Elist$rawdata[,5] 
      flagdata[,3]=flagdata[,2]-flagdata[,1] 
      flagdata[,4]=flagdata[,2]/flagdata[,1] 
      flagdata[,5]=flagdata[,2]-Elist$rawdata[,1] 
      flagdata[,6]=flagdata[,2]/Elist$rawdata[,1] 
   } else { 
      flagdata[,1]=(Elist$rawdata[,2]+Elist$rawdata[,4])*0.5
      flagdata[,2]=(Elist$rawdata[,3]+Elist$rawdata[,5])*0.5
      flagdata[,3]=flagdata[,2]-flagdata[,1] 
      flagdata[,4]=flagdata[,2]/flagdata[,1] 
      flagdata[,5]=flagdata[,2]-Elist$rawdata[,1] 
      flagdata[,6]=flagdata[,2]/Elist$rawdata[,1] 
   }
   
   
   
   
   flagdata=cbind(flagdata,Num ) 
   flagNum=flagdata[flagdata[,8]== -100,9]
  
  if(koc==1){
     flagdata[flagNum,1:8]=NA
   }
   
   
  background = matrix(flagdata[,1],nrow=Elist$No_antigen,ncol=Elist$No_sample) 
  foreground = matrix(flagdata[,2],nrow=Elist$No_antigen,ncol=Elist$No_sample) 
  intensity=matrix(flagdata[,3],nrow=Elist$No_antigen,ncol=Elist$No_sample) 
  ratio=matrix(flagdata[,4],nrow=Elist$No_antigen,ncol=Elist$No_sample) 
  intensity_neighbour=matrix(flagdata[,5],nrow=Elist$No_antigen,ncol=Elist$No_sample) 
  ratio_neighbour=matrix(flagdata[,6],nrow=Elist$No_antigen,ncol=Elist$No_sample) 
  snr=matrix(flagdata[,7],nrow=Elist$No_antigen,ncol=Elist$No_sample) 
  
  colnames(background)=colnames(foreground)=colnames(intensity_neighbour)=colnames(ratio_neighbour)=colnames(intensity)=colnames(ratio)=colnames(snr)=Elist$target[,3] 
  rownames(background)=rownames(foreground)=rownames(intensity_neighbour)=rownames(ratio_neighbour)=rownames(intensity)=rownames(ratio)=rownames(snr)=Elist$annotation[1:Elist$No_antigen,4]
  
  intensity=No_minus(intensity)
  intensity_neighbour=No_minus(intensity_neighbour) 
  
  ##aggregate duplicate spots   
  meanIntensity=avereps(intensity) 
  meanRatio=avereps(ratio) 
  meanIntensity_neighbour=avereps(intensity_neighbour) 
  meanRatio_neighbour=avereps(ratio_neighbour)   
  meanSNR=avereps(snr) 
  meanBackground=avereps(background)
  meanForeground=avereps(foreground)
  
  if(koc==1){
      
     # solve NA effect
      meanIntensity[is.na(meanIntensity)] <- 0.001 
      meanRatio[is.na(meanRatio)] <- 0.001 
      meanIntensity_neighbour[is.na(meanIntensity_neighbour)] <- 0.001 
      meanRatio_neighbour[is.na(meanRatio_neighbour)] <- 0.001 
      meanSNR[is.na(meanSNR)] <- 0.001 
      meanBackground[is.na(meanBackground)] <- 0.001 
      meanForeground[is.na(meanForeground)] <- 0.001 
    ## remove bad samples. 
    if( length( which(colnames(meanIntensity)=="Bad")) > 0 ){ 
      meanIntensity=meanIntensity[,-which(colnames(meanIntensity)=="Bad")] 
      meanRatio=meanRatio[,-which(colnames(meanRatio)=="Bad")] 
      meanIntensity_neighbour=meanIntensity_neighbour[,-which(colnames(meanIntensity_neighbour)=="Bad")] 
      meanRatio_neighbour=meanRatio_neighbour[,-which(colnames(meanRatio_neighbour)=="Bad")] 
      meanSNR=meanSNR[,-which(colnames(meanSNR)=="Bad")] 
      meanBackground=meanBackground[,-which(colnames(meanBackground)=="Bad")]
      meanForeground=meanForeground[,-which(colnames(meanForeground)=="Bad")]
    }
    
    ## sort names by row/col and put all positive Abs in the bottom 
      Intensity=BottomPositive(meanIntensity)$newdata 
      Ratio=BottomPositive(meanRatio)$newdata  
      Intensity_neighbour=BottomPositive(meanIntensity_neighbour)$newdata  
      Ratio_neighbour=BottomPositive(meanRatio_neighbour)$newdata  
      SNR=BottomPositive(meanSNR)$newdata  
      Background=BottomPositive(meanBackground)$newdata  
      Foreground=BottomPositive(meanForeground)$newdata   
      
      ## abstract subset of PBS 
      
     if( length( which(colnames(Intensity)=="PBS")) > 0 ) { 
     
      PBS_intensity=Intensity[,which(colnames(Intensity)=="PBS")] 
      PBS_ratio=Ratio[,which(colnames(Ratio)=="PBS")]
      PBS_ratio=as.matrix(PBS_ratio) 
      rownames(PBS_ratio)=rownames(Intensity) 
      colnames(PBS_ratio)=rep("PBS",ncol(PBS_ratio)) 
      PBS_intensity_neighbour=Intensity_neighbour[,which(colnames(Intensity_neighbour)=="PBS")] 
      PBS_ratio_neighbour=Ratio_neighbour[,which(colnames(Ratio_neighbour)=="PBS")] 
      PBS_snr=SNR[,which(colnames(SNR)=="PBS")] 
      PBS_background=Background[,which(colnames(Background)=="PBS")] 
      PBS_foreground=Foreground[,which(colnames(Foreground)=="PBS")] 
      
      nonPBS_intensity=Intensity[,-which(colnames(Intensity)=="PBS")] 
      nonPBS_ratio=Ratio[,-which(colnames(Ratio)=="PBS")] 
      nonPBS_intensity_neighbour=Intensity_neighbour[,-which(colnames(Intensity_neighbour)=="PBS")] 
      nonPBS_ratio_neighbour=Ratio_neighbour[,-which(colnames(Ratio_neighbour)=="PBS")] 
      nonPBS_snr=SNR[,-which(colnames(SNR)=="PBS")] 
      nonPBS_background=Background[,-which(colnames(Background)=="PBS")] 
      nonPBS_foreground=Foreground[,-which(colnames(Foreground)=="PBS")] 
     
      pbs_intensity_ave=avePBS(PBS_intensity) 
      pbs_ratio_ave=avePBS(PBS_ratio) 
      pbs_intensity_neighbour_ave=avePBS(PBS_intensity_neighbour) 
      pbs_ratio_neighbour_ave=avePBS(PBS_ratio_neighbour) 
      pbs_snr_ave=avePBS(PBS_snr) 
      pbs_background_ave=avePBS(PBS_background) 
      pbs_foreground_ave=avePBS(PBS_foreground) 
           
      nonPBS_minimum = backgroundCorrect.matrix(nonPBS_foreground, nonPBS_background, method="minimum", offset=0, printer=NULL,
                                                normexp.method="saddle", verbose=TRUE)
      nonPBS_normexp = backgroundCorrect.matrix(nonPBS_foreground, nonPBS_background, method="normexp", offset=0, printer=NULL,
                                                normexp.method="saddle", verbose=TRUE)
      pbs_minimum = backgroundCorrect.matrix(pbs_foreground_ave, pbs_background_ave, method="minimum", offset=0, printer=NULL,
                                                normexp.method="saddle", verbose=TRUE)
      pbs_normexp = backgroundCorrect.matrix(pbs_foreground_ave, pbs_background_ave, method="normexp", offset=0, printer=NULL,
                                                normexp.method="saddle", verbose=TRUE)
                                             
                                             
      TT=nrow(nonPBS_ratio)-BottomPositive(meanIntensity)$No_positive
      
      for( i in 1:TT) { 
        nonPBS_intensity[i,]=nonPBS_intensity[i,]-pbs_intensity_ave[i] 
        nonPBS_ratio[i,]=nonPBS_ratio[i,]/pbs_ratio_ave[i] 
        nonPBS_intensity_neighbour[i,]=nonPBS_intensity_neighbour[i,]-pbs_intensity_neighbour_ave[i] 
        nonPBS_ratio_neighbour[i,]=nonPBS_ratio_neighbour[i,]/pbs_ratio_neighbour_ave[i]  
        nonPBS_minimum[i,]=nonPBS_minimum[i,]-pbs_minimum[i]
        nonPBS_normexp[i,]=nonPBS_normexp[i,]-pbs_normexp[i]
      } 
      
      nonPBS_intensity=No_minus(nonPBS_intensity) 
      nonPBS_ratio=No_minus(nonPBS_ratio) 
      nonPBS_intensity_neighbour=No_minus(nonPBS_intensity_neighbour) 
      nonPBS_ratio_neighbour=No_minus(nonPBS_ratio_neighbour) 
      
      elist$pbs_signal=PBS_intensity 
      elist$pbs_ratio=PBS_ratio 
      elist$pbs_intensity_neighbour=PBS_intensity_neighbour
      elist$pbs_ratio_neighbour= PBS_ratio_neighbour
      
     } else { 
      
       nonPBS_intensity=Intensity
       nonPBS_ratio =Ratio 
       nonPBS_intensity_neighbour=Intensity_neighbour
       nonPBS_ratio_neighbour=Ratio_neighbour
       nonPBS_snr=SNR
       nonPBS_minimum=backgroundCorrect.matrix(Foreground, Background, method="minimum", offset=0, printer=NULL,
                                                normexp.method="saddle", verbose=TRUE)
       nonPBS_normexp=backgroundCorrect.matrix(Foreground, Background, method="normexp", offset=0, printer=NULL,
                                                normexp.method="saddle", verbose=TRUE)
      
     }
    
          
  } 
  
  # End of quality control     
   
 
  layout= Elist$target 
  attach(layout) 
  # do batch correction with users' instruction
  
        layout_1 = layout[ - which(layout$sampleID=="Bad") , ]  
        layout_2 = layout_1[ - which(layout_1$sampleID=="PBS"), ] 
            
     if(batch_method != 0) {
         meanIntensity=Batch_correction_Method(meanIntensity, layout,batch_method)
         meanRatio=Batch_correction_Method(meanRatio, layout, batch_method)
         meanIntensity_neighbour=Batch_correction_Method(meanIntensity_neighbour,layout,batch_method)
         meanRatio_neighbour=meanIntensity_neighbour=Batch_correction_Method(meanRatio_neighbour, layout, batch_method)
      } else {
      
       meanIntensity=nonPBS_intensity
       meanRatio=nonPBS_ratio 
       meanIntensity_neighbour=nonPBS_intensity_neighbour 
       meanRatio_neighbour=nonPBS_ratio_neighbour
           
      }
    
    # do normalization with users' instruction
       
    
      if(norm_method != 0) { 
       Norm_intensity=Normalization_Method(nonPBS_intensity,nonPBS_background,  layout,norm_method) 
       Norm_ratio=Normalization_Method( nonPBS_ratio,nonPBS_background, layout,norm_method) 
       Norm_signal_neighbour=Normalization_Method(nonPBS_intensity_neighbour,nonPBS_background,  layout,norm_method)
       Norm_ratio_neighbour=Normalization_Method(nonPBS_ratio_neighbour,nonPBS_background, layout,norm_method)
      } else {
      
       Norm_intensity=meanIntensity
       Norm_ratio=meanRatio 
       Norm_signal_neighbour=meanIntensity_neighbour 
       Norm_ratio_neighbour=meanRatio_neighbour
     }
   
   
   
  elist$signal=No_minus(Norm_intensity) 
  elist$ratio=No_minus(Norm_ratio) 
  elist$signal_neighbour=No_minus(Norm_signal_neighbour) 
  elist$ratio_neighbour=No_minus(Norm_ratio_neighbour) 
  
  elist$snr=No_minus(nonPBS_snr) 
  elist$minimum=No_minus(nonPBS_minimum) 
  elist$normexp=No_minus(nonPBS_normexp) 
  elist$No_sample=ncol(nonPBS_ratio)
     
  elist$No_Ab=nrow(nonPBS_ratio) 
  elist$score=log2(elist$signal * elist$snr +1) 
  #elist$score=apply(target[,1:2],2,as.numeric) 
  elist$score_neighbour=log2(elist$signal_neighbour * elist$snr+1) 
  
 
  
  
  return (elist) 
}   


Make_RepList=function(raw_data, NAME, num_anti, target, aggre_way){
  foreground = makeMatrix(raw_data[,3],NAME,num_anti,target)
  background = makeMatrix(raw_data[,2],NAME,num_anti,target)
  minimum = backgroundCorrect.matrix(foreground, background, method="minimum", offset=0, printer=NULL,
                                     normexp.method="saddle", verbose=TRUE)
  normexp = backgroundCorrect.matrix(foreground, background, method="normexp", offset=0, printer=NULL,
                                     normexp.method="saddle", verbose=TRUE)
  Final=list(F_B_sub=makeMatrix(raw_data$F_B,NAME,num_anti,target), F_B_div = makeMatrix(raw_data$FB_ratio,NAME,num_anti,target), SNR=makeMatrix(raw_data$SNR,NAME,num_anti,target), min= minimum, normexp= normexp)
  Final$F_B_sub = order_rearrange(Final$F_B_sub)
  Final$F_B_div = order_rearrange(Final$F_B_div)
  Final$SNR = order_rearrange(Final$SNR)
  Final$min = order_rearrange(Final$min)
  Final$normexp = order_rearrange(Final$normexp)
  return(Final)
}

#Caculate Table Summary
Table_Sum=function(mtx, type){ 
  idName = "SampleID"
  if(type==1){
    mtx = t(mtx)
  }
  D2 = as.data.frame(mtx)
  row_names = rownames(mtx)
  col_names= colnames(mtx)
  mean = colMeans(mtx, na.rm = FALSE, dims = 1)
  min = colQuantiles(mtx, probs = 0)
  first_quantile = colQuantiles(mtx, probs = 0.25)
  median = colQuantiles(mtx, probs = 0.5)
  third_quantile = colQuantiles(mtx, probs = 0.75)
  max = colQuantiles(mtx, probs = 1)
  table_sum <- data.frame("SampleID" = col_names, "Min" = min, "FirstQuantile" = first_quantile, "Median"=median, "Mean"= mean,  "ThirdQuantile"=third_quantile, "Max"=max)
  if(type==1){
    names(table_sum)[names(table_sum) == "SampleID"] <- "AntigenID"
  }
  return (table_sum) 
}



server <- shinyServer(function(input, output,session) {
  
  
  # import 
  d=data(package='datasets')$results[,'Item']
  d=d[!grepl('[\\()]',d)]
  d=d[!d%in%c('UScitiesD','eurodist','sleep','warpbreaks')]
  d=d[unlist(lapply(d,function(d.in) eval(parse(text=paste0('ncol(as.data.frame(datasets::',d.in,'))')))))>1]
  d=d[-which(d=='mtcars')]
  d=c('mtcars',d)
  TEMPLIST<-new.env()
  TEMPLIST$d <- d
  #############################################   Data Upload   #####################################################
  
  #return the the data frame of layout file
  getTarget <- reactive({
    if (is.null(input$file2)){
      return(NULL)
    }
    layout = read.table(input$file2[[1, 'datapath']], header=T, sep="\t")
    layout_code=seq(from=1,to=nrow(layout),by=1)
    target = layout[layout_code,]
    # attach(target)
    return(target)
  })
  
  
  #read in gpr files from user import and return rawdata needed (loadGPR)
  getGPR <- reactive({
    proteinType = input$pro_type
    inFile <- input$file1
    if (is.null(inFile)){
      return(NULL)
    }else {
      inFile = input$file1[naturalorder(input$file1$name),]
      numfiles = nrow(inFile)
      sample_layout=getTarget()
      No_slide=length(unique(sample_layout[,1]))    
      elist=list()
      elist$target=sample_layout 
      elist$No_slide=length(unique(elist$target[,1])) 
      alldata=c() 
      Bneighbour_Median =c() 
      for( s in 1: numfiles) {
        row2nd=read.table(inFile[[s, 'datapath']], skip=1, nrow=1,header=F, sep="\t")
        skipRow=as.numeric(row2nd[1])+2 
        slidedata=read.table(inFile[[s, 'datapath']], header = TRUE,sep="\t",skip=skipRow) 
        slide_neighbour_median=Neighbour(slidedata,proteinType) 
        Bneighbour_Median=c(Bneighbour_Median,slide_neighbour_median) 
        alldata=rbind(alldata,slidedata)  
      }
      elist$annotation=alldata[,1:5] 
      elist$Block_per_slide=length(unique(elist$annotation[,1])) 
      column_per_block=length(unique(elist$annotation[,2])) 
      row_per_block=length(unique(elist$annotation[,3])) 
      elist$No_antigen=row_per_block*column_per_block
      elist$antigenID=elist$annotation[1:elist$No_antigen,4] 
      elist$No_sample=No_slide * elist$Block_per_slide 
      attach(alldata) 
      if(proteinType==532) { 
              elist$rawdata=data.frame(Bneighbour_Median,B532_Median=B532.Median, F532_Median=F532.Median,B532_Mean=B532.Mean,F532_Mean=F532.Mean, SNR=SNR.532, Flag=Flags)
            }else { elist$rawdata=data.frame(Bneighbour_Median,B635_Median=B635.Median, F635_Median=F635.Median,B635_Mean=B635.Mean, F635_Mean=F635.Mean,SNR=SNR.635,Flag=Flags)
      }
      return(elist)
    }
  })
  
  
 #construct summmary table of general infomation(displayed in Data Upload)
   basicInfo <- reactive({
     m = matrix(nrow=9,ncol=1)
     rownames(m) <- c("Number of Slides", "Number of bad samples", "Number of PBS controls", "Number of good samples", "Number of total positive controls","Number of Ig control","Number of anti-Ig controls","Number of negative control","Number of antibodies")
     colnames(m) <- c("Quantity")
     m[,1]=c(0,0,0,0,0,0,0,0,0)
     if(!is.null(input$file1)&(!is.null(input$file2))){
       Mtx_list = get_Mtx_list()
       sig = Mtx_list$signal
       target=getTarget()
       num_anti = nrow(sig)
       m[5,1]=BottomPositive(sig)$No_positive
       m[6,1]=BottomPositive(sig)$No_Ig_control 
       m[7,1]=BottomPositive(sig)$No_anti_Ig
       m[8,1]= 0 
       m[9,1]= num_anti-m[5,1]-m[8,1]
     }
     if(!is.null(input$file2)){
       target=getTarget()
       m[1,1]=max(target$Slide)
       pbs <- target[ which(target$sampleID=='PBS'), ]
       bad <- target[ which(target$sampleID=='Bad'), ]
       m[2,1]=as.numeric(dim(bad)[1])
       m[3,1]=as.numeric(dim(pbs)[1])
       m[4,1]=nrow(target)-m[3,1]-m[2,1] 
     }
     return(m)
   })
  
  
  
  #display basic info Table 
  output$basicInfo <- renderTable({
    basicInfo()
  }, rownames = TRUE, digits=0)
  
  
  
    #return the list of all matrix needed
      get_Mtx_list <-reactive({
             target = getTarget()
             Mtx  =   getGPR()
             Mtx  = Pre_process(Mtx,input$measure_type, input$qc_type, input$batch_type,input$norm_type)
      })
    
  
  
  #############################################   QC report   #####################################################                                  
    #obtain the control tables
    getControls <- reactive({
       
      sig = get_Mtx_list()$signal
      sortedsignal=BottomPositive(sig)
      
      mtx=sortedsignal$newdata 
      n.controls=sortedsignal$No_positive 
      
      n.Igcontrol=sortedsignal$No_Ig_control 
      n.AntiIg=sortedsignal$No_Ig_anti_Ig
      
      TX=nrow(mtx)-n.controls+1  
      control=mtx[TX:nrow(mtx),] 
      
      if(n.controls <3) {
        pc_control=control[2,]
        pc_control=t(as.matrix(pc_control)) 
        rownames(pc_control)="Ig control"
        pc_anti=control[1,] 
        pc_anti=t(as.matrix(pc_anti)) 
        rownames(pc_anti)="Anti_Ig control"
        
      }else {
        pc_anti=control[grepl("nti-Ig", rownames(control)), ]
        pc_control=control[-grep("nti-Ig",rownames(control)),]
      
      }
       
       pbs_ratio = get_Mtx_list()$pbs_ratio
           
           
           removeControl=control
           removeControl=avereps(removeControl) 
           removeName =row.names(removeControl)
           
           n.controls=nrow(removeControl) 
           
           n.pbs=ncol(pbs_ratio) 
           
           if(n.pbs==1) {
             TT=nrow(pbs_ratio)-n.controls 
             pbs_ratio=t(pbs_ratio[1:TT,]) 
           } else {
             pbs_ratio =t( pbs_ratio[!row.names(pbs_ratio)%in%removeName,]) 
      }
      
      
      pc_control_to_graph <- switch(input$graph_type, '1'=pc_control, '2'=pc_anti, '3'=pbs_ratio)
      return(pc_control_to_graph)
    })
    
    #generate control table
    output$ControlsTable <- renderTable({
      getControls()
    }, rownames = TRUE, digits=0)
    
    #generate line chart graph
    output$line_chart_plot <- renderPlot({
      input$qc_graph
      mtx=getControls()
      if(input$graph_type<3){
        matplot(t(mtx), type = "l", lty = 6, lwd = 3, lend = par("lend"),
                pch = NULL,
                col = 1:100, cex = NULL, bg = NA,
                xlab = "Sample Index", ylab = "Values" , xlim = NULL, ylim = NULL,
                log = "", add = FALSE, verbose = getOption("verbose"))
        LegendPosition <- switch(input$legend_position, '0'= "topleft", '1'=0)
        if(LegendPosition!=0){
           if(nrow(mtx)<2) {
             legend(LegendPosition, legend=rownames(mtx), title="Names of the Controls",pch=15, horiz=TRUE, col=1:as.numeric(ncol(mtx)[1]))
           
           } else {
            legend(LegendPosition, inset=.0, legend=rownames(mtx), title="Names of the Controls",pch=15, horiz=TRUE, col=1:as.numeric(ncol(mtx)[1]))
            }
        }
        target=getTarget()
        target[,1:2]=apply(target[,1:2],2,as.numeric) 
        num_slides=max(target[,1])
          
        num_samples=length(target$Slide)
        samples_cumulative=0
        for (i in 1:num_slides){
          #number of samples in a slide - number of PBS in a slide -number of Bad samples in a slide
          total_valid_sample_inSlide = length(which(target$Slide==i)) - length(which((target$Slide==i)&target$sampleID=='PBS')) - length(which((target$Slide==i)&target$sampleID=='Bad'))
          samples_cumulative= samples_cumulative + total_valid_sample_inSlide
          slide_number ='Slide'
          slide_number = paste(slide_number,i, sep="")
          slide_number = paste(slide_number,"end", sep=" ")
          abline(h = 10, v = samples_cumulative, col = "gray60", lwd=3, lty=2)
          #default make slide indicator at 25% of the max value
          show_slide_indicator= switch(input$legend_position, '0'=25, '1'=NA)
          text(samples_cumulative,(show_slide_indicator/100)*getControls_max_value(),slide_number)
        }
      }
      else{
        #outliers not shown
        boxplot(getControls(), use.cols = TRUE, outline=TRUE, main="PBS Ratio Boxplot", col=(c("pink")))
        abline(h=1, col="red")
         
     }
    })
    
    
    #generate graph of controls
    Plot_control<- reactive({
        mtx=getControls()
        if(input$graph_type<3){
          matplot(t(mtx), type = "l", lty = 6, lwd = 3, lend = par("lend"),
                  pch = NULL,
                  col = 1:100, cex = NULL, bg = NA,
                  xlab = "Sample Index", ylab = "Values" , xlim = NULL, ylim = NULL,
                  log = "", add = FALSE, verbose = getOption("verbose"))
          LegendPosition <- switch(input$legend_position, '0'= "topleft", '1'=0)
          if(LegendPosition!=0){
              if(nrow(mtx)<2 ) {
	                legend(LegendPosition, legend=rownames(mtx), title="Names of the Controls",pch=15, horiz=TRUE, col=1:as.numeric(ncol(mtx)[1]))
	               
	               } else {
	                legend(LegendPosition, inset=.0, legend=rownames(mtx), title="Names of the Controls",pch=15, horiz=TRUE, col=1:as.numeric(ncol(mtx)[1]))
            }
          }
          target=getTarget()
          target[,1:2]=apply(target[,1:2],2,as.numeric) 
          num_slides=max(target[,1])
          num_samples=length(target$Slide)
          samples_cumulative=0
          for (i in 1:num_slides){
            #number of samples in a slide - number of PBS in a slide -number of Bad samples in a slide
            total_valid_sample_inSlide = length(which(target$Slide==i)) - length(which((target$Slide==i)&target$sampleID=='PBS')) - length(which((target$Slide==i)&target$sampleID=='Bad'))
            samples_cumulative= samples_cumulative + total_valid_sample_inSlide
            slide_number ='Slide'
            slide_number = paste(slide_number,i, sep="")
            slide_number = paste(slide_number,"end", sep=" ")
            abline(h = 10, v = samples_cumulative, col = "gray60", lwd=3, lty=2)
            #default make slide indicator at 25% of the max value
            show_slide_indicator= switch(input$legend_position, '0'=25, '1'=NA)
            text(samples_cumulative,(show_slide_indicator/100)*getControls_max_value(),slide_number)
          }
        }
        else{
          boxplot(getControls(), use.cols = TRUE, outline=TRUE, main="PBS Ratio Boxplot", col=(c("pink")))
          abline(h=1, col="red")
           
        }
         
    })
    
    
   #get max value from the control matrix
    getControls_max_value <- reactive({
      max_value=10000
      if(!is.null(input$file1)&(!is.null(input$file2))){
        mtx = getControls()
        max = apply(mtx, c(1,2), max)
        max_value=max(max)
      }
      return(max_value)
    })
  
  #download graph
    output$download_control_graph <- downloadHandler(
      filename = function() { paste('control_analysis', '.png', sep='') },
      content = function(file) {
        device <- function(..., width, height) grDevices::png(..., width = 2*width, height = 1.4*height, res = 300, units = "in")
        ggsave(file, plot = Plot_control(), device = device)
      }
    )
  #obtain antigen's table summary
    get_Table_Summary <- reactive({
      input$generate_table_sum
      List = get_Mtx_list()
      ##use substraction only
      mtx = List$signal
      Table_Sum(mtx, input$summary_table_type)
    })
    
    
    output$Table_Summary <- renderTable(
      get_Table_Summary(), align = 'c', width=920)
    
    #table summary download type
    output$antibodies_table_summary_data <- downloadHandler(
      filename = function() {
        save_type = switch(input$antibodies_table_summary_download_type, '0'= '.csv', '1'='.xlsx', '2'='.txt')
        paste('data-', Sys.Date(), save_type, sep='')
      },
      
      content = function(file) {
        if(input$antibodies_table_summary_download_type==2){
          write.table(get_Table_Summary(), file, sep = "\t", row.names = FALSE)
        }
        if(input$antibodies_table_summary_download_type==0){
          write.csv(get_Table_Summary(), file, row.names = FALSE)
        }
        if(input$antibodies_table_summary_download_type==1){
          write.xlsx(get_Table_Summary(),file,sheetName = "Table_Summary",col.names = TRUE,row.names = FALSE)
        }
      }
  )
  
  
  
#############################################   Data Preprocessing   #####################################################
 
    
  #provide users the matrix they wanted according to the choice from radio button(input$protein_type, input$data_type)
  #trigger by click show action button
  
  getMatrix <- reactive({
    input$show
    List = get_Mtx_list()
    Signal_Intensity <- switch(input$protein_type, '0'=List$signal, '1'=List$ratio, '2'=List$minimum, '3'=List$normexp, '4'=List$signal_neighbour)
    data_matrix <- switch(input$data_type, "Signal"=Signal_Intensity, "SNR"=List$snr, "Ab_score"= List$score)
    threshold = input$decimal_percent#for threshld
    
    vac = c()#an empty for hold delete
    SNR3 = c() #hold for SNR>3
    # loop all and find all > 3
    for(row in 1:nrow(List$snr)){
      total_col = 0
      total_greather3 = 0
      for(col in 1:ncol(List$snr)){
        total_col = total_col + 1
        if(List$snr[row,col]>=3){
          total_greather3 = total_greather3 + 1
        }
      }
      percent_3 = total_greather3/total_col * 100
      SNR3 <- c(SNR3, percent_3)
    }
    data_matrix <- cbind(data_matrix, c(SNR3)) #add a new col for SNR>3 rate
    for(row in 1:nrow(data_matrix)){
      if(data_matrix[row,ncol(data_matrix)]<threshold){
        vac <- c(vac,row) #append
      }
    }
    if(threshold != 0){
      data_matrix <- data_matrix[-c(vac), ]
    }
    colnames(data_matrix)[ncol(data_matrix)] <- "SNR%>=threshld"
    #delete that row
    return(data_matrix)
  })
   
    
  #output the matrix wanted, with rowname and 3 significant figures
  output$matrix <- renderTable({
    getMatrix()
  },rownames = TRUE, digits=3)
  
   
  #make a progressbar for generating matrix
  observeEvent(input$show, {
    for (i in 1:100) {
      updateProgressBar(
        session = session,
        id = "pb2",
        value = i, total = 100,
        title = paste("Your matrix is generating")
      )
      #sys sleep time
      # Sys.sleep(0.01)
    }
  })
  
  
  #check if the files are uploaded to display error message
  output$fileUploaded <- reactive({
    return((!is.null(getGPR()))&(!is.null(getTarget())))
  })
  outputOptions(output, 'fileUploaded', suspendWhenHidden=FALSE)
   
  #Download button to download files with corresponding matrix and name
  output$preprocessing_data_download <- downloadHandler(
    filename = function() {
      save_type = switch(input$preprocessing_download_type, '0'= '.csv', '1'='.xlsx', '2'='.txt')
      paste(input$data_type, save_type, sep='')
    },
    
    content = function(file) {
      if(input$preprocessing_download_type==2){
        write.table(getMatrix(), file, sep = "\t", quote=FALSE,row.names = TRUE)
      }
      if(input$preprocessing_download_type==0){
        write.csv(getMatrix(), file, quote=FALSE,row.names = TRUE)
      }
      if(input$preprocessing_download_type==1){
        write.xlsx(getMatrix(),file,sheetName = input$data_type, row.names = TRUE)
      }
    }
  )
  
  #get the matrix chosen by the user
  datasetInput <- reactive({
    List = get_Mtx_list()
    Signal_Intensity <- switch(input$protein_type, '0'=List$signal, '1'=List$ratio, '2'=List$minimum, '3'=List$normexp, '4'=List$signal_neighbour)
    mtx = switch(input$matrix_given,
                 '0' = Signal_Intensity,
                 '1'=  List$ratio,
                 '2' = List$snr,
                 '3' = List$score)
    return(mtx)
  })
  
  
  # Use DT table, to display matrix chosen in "Data Preprocessing"
  output$tables2=DT::renderDataTable(round(getMatrix(),2),server = T,filter='top',
                                     extensions = c('Scroller','FixedHeader','FixedColumns','Buttons','ColReorder'),
                                     options = list(
                                       dom = 't',
                                       buttons = c('copy', 'csv', 'excel', 'pdf', 'print','colvis'),
                                       colReorder = TRUE,
                                       scrollX = TRUE,
                                       fixedColumns = TRUE,
                                       fixedHeader = TRUE,
                                       deferRender = TRUE,
                                       scrollY = 725,
                                       scroller = TRUE
                                     ))
                                     
  
 
  
#############################################   Heatmap   #####################################################  
#Annotation Variable UI ----
observeEvent(data.sel(),{
  output$annoVars<-renderUI({
    data.in=data.sel()
    NM=NULL
    
    if(any(sapply(data.in,class)=='factor')){
      NM=names(data.in)[which(sapply(data.in,class)=='factor')]  
    } 
    column(width=4,
           selectizeInput('annoVar','Annotation',choices = names(data.in),selected=NM,multiple=T,options = list(placeholder = 'select columns',plugins = list("remove_button")))
    )
  })

#Sampling UI ----  
  output$sample<-renderUI({
    list(
      column(4,textInput(inputId = 'setSeed',label = 'Seed',value = sample(1:10000,1))),
      column(4,numericInput(inputId = 'selRows',label = 'Number of Rows',min=1,max=pmin(500,nrow(data.sel())),value = pmin(500,nrow(data.sel())))),
      column(4,selectizeInput('selCols','Columns Subset',choices = names(data.sel()),multiple=T))
    )
  })
})


#Data Selection UI ----
output$data=renderUI({
  if(!is.null(input$mydata)) TEMPLIST$d=c(input$mydata$name,TEMPLIST$d)
  selData=head(TEMPLIST$d,1)
  selectInput("data","Select Data",TEMPLIST$d,selected = selData)
})


#Color Pallete UI ----
output$colUI<-renderUI({
  
  colSel='Vidiris'
  if(input$transform_fun=='cor') colSel='RdBu'
  if(input$transform_fun=='is.na10') colSel='grey.colors'
  
  selectizeInput(inputId ="pal", label ="Select Color Palette",
                 choices = c('Vidiris (Sequential)'="viridis",
                             'Magma (Sequential)'="magma",
                             'Plasma (Sequential)'="plasma",
                             'Inferno (Sequential)'="inferno",
                             'Magma (Sequential)'="magma",
                             'Magma (Sequential)'="magma",
                             
                             'RdBu (Diverging)'="RdBu",
                             'RdYlBu (Diverging)'="RdYlBu",
                             'RdYlGn (Diverging)'="RdYlGn",
                             'BrBG (Diverging)'="BrBG",
                             'Spectral (Diverging)'="Spectral",
                             
                             'BuGn (Sequential)'='BuGn',
                             'PuBuGn (Sequential)'='PuBuGn',
                             'YlOrRd (Sequential)'='YlOrRd',
                             'Heat (Sequential)'='heat.colors',
                             'Grey (Sequential)'='grey.colors'),
                 selected=colSel)
})

#Manual Color Range UI ----
output$colRng=renderUI({
  if(!is.null(data.sel())) {
    rng=range(data.sel(),na.rm = TRUE)
  }else{
    rng=range(mtcars) # TODO: this should probably be changed
  }
  sliderInput("colorRng", "Set Color Range", min = round(rng[1],1), max = round(rng[2],1), step = .1, value = rng)
  n_data = nrow(data.sel())
  
  min_min_range = ifelse(input$transform_fun=='cor',-1,-Inf)
  min_max_range = ifelse(input$transform_fun=='cor',1,rng[1])
  min_value = ifelse(input$transform_fun=='cor',-1,rng[1])
  
  max_min_range = ifelse(input$transform_fun=='cor',-1,rng[2])
  max_max_range = ifelse(input$transform_fun=='cor',1,Inf)
  max_value = ifelse(input$transform_fun=='cor',1,rng[2])
  
  a_good_step = 0.1 # (max_range-min_range) / n_data
  
  list(
    numericInput("colorRng_min", "Set Color Range (min)", value = min_value, min = min_min_range, max = min_max_range, step = a_good_step),
    numericInput("colorRng_max", "Set Color Range (max)", value = max_value, min = max_min_range, max = max_max_range, step = a_good_step)
  )
})

 
  button_values <- reactiveValues(example = 0, previous = 0)
  observeEvent(input$show2, {
    button_values$previous <- 1
    button_values$example <- 0
  })

  observeEvent(input$show3, {
    button_values$previous <- 0
    button_values$example <- 1
  })


  #button for 3 given matrix generate
  # data.sel=eventReactive(input$data,{
  #   if(input$data%in%d){
  #     List = get_Mtx_list()
  #     bg_correction <- switch(input$protein_type, '0'=List$signal, '1'=List$ratio, '2'=List$minimum, '3'=List$normexp, '4'=List$signal_neighbour)
  #     data.in = switch(input$matrix_given,
  #                      '0' = bg_correction,
  #                      '1' = List$snr,
  #                      '2' = List$score)
  #   }else{
  #     data.in = read.csv(input$mydata[[1, 'datapath']], sep=",", row.names=1)
  #     data.in= as.matrix(data.in)
  #   }
  #   return(data.in)
  # })
  
  
  #read data from user or dataset
  data.sel=eventReactive(input$data,{
    if(input$data%in%d){
      eval(parse(text=paste0('data.in=as.data.frame(datasets::',input$data,')')))
    }else{
      data.in=importSwitch(input$mydata[input$mydata$name%in%input$data,])
      data.in = read.csv(input$mydata[[1, 'datapath']], sep=",", row.names=1)
      data.in= as.matrix(data.in)
    }
    data.in=as.data.frame(data.in)
    # data.in=data.in[,sapply(data.in,function(x) class(x))%in%c('numeric','integer')] # no need for this
    return(data.in)
  })

  
  
  
#Building heatmaply ----
  
interactiveHeatmap<- reactive({
  data.in=data.sel()
  if(input$showSample){
    if(!is.null(input$selRows)){
        set.seed(input$setSeed)
      if((input$selRows >= 2) & (input$selRows < nrow(data.in))){
        # if input$selRows == nrow(data.in) then we should not do anything (this save refreshing when clicking the subset button)
        if(length(input$selCols)<=1) data.in=data.in[sample(1:nrow(data.in),pmin(500,input$selRows)),]
        if(length(input$selCols)>1) data.in=data.in[sample(1:nrow(data.in),pmin(500,input$selRows)),input$selCols]
      }
    }
  }
  # ss_num = sapply(data.in,function(x) class(x)) %in% c('numeric','integer') # in order to only transform the numeric values
  
  if(length(input$annoVar)>0){
    if(all(input$annoVar%in%names(data.in))) 
      data.in <- data.in%>%mutate_at(funs(factor),.vars=vars(input$annoVar))
  } 
  ss_num =  sapply(data.in, is.numeric) # in order to only transform the numeric values
  
  if(input$transpose) data.in=t(data.in)
  if(input$transform_fun!='.'){
    if(input$transform_fun=='is.na10'){
      updateCheckboxInput(session = session,inputId = 'showColor',value = T)
      data.in[, ss_num]=is.na10(data.in[, ss_num])
    } 
    if(input$transform_fun=='cor'){
      updateCheckboxInput(session = session,inputId = 'showColor',value = T)
      updateCheckboxInput(session = session,inputId = 'colRngAuto',value = F)
      data.in=cor(data.in[, ss_num],use = "pairwise.complete.obs")
    }
    
    if(input$row_or_col==1){
      if(input$transform_fun=='mean')
      {
        data.in = sweep(data.in,1,apply(data.in,1,mean), '-')
      }
      if(input$transform_fun=='sd') {
        data.in = sweep(data.in,1,apply(data.in,1,sd),'/')
      }
      if(input$transform_fun=='normalize') {
        
        datamean = apply(data.in,1,mean)
        datasd = apply(data.in,1,sd)
        data.in = sweep(data.in,1,datamean,'-')
        data.in = sweep(data.in,1,datasd,'/')
      }
      if(input$transform_fun=='median') {
        data.in = sweep(data.in,1,apply(data.in,1,median), '-')
      }
    }
    if(input$row_or_col==2){
      if(input$transform_fun=='mean')
      {
        data.in = sweep(data.in,2,apply(data.in,2,mean), '-')
      }
      if(input$transform_fun=='sd') {
        data.in = sweep(data.in,2,apply(data.in,2,sd),'/')
      }
      if(input$transform_fun=='normalize') {
        
        datamean = apply(data.in,2,mean)
        datasd = apply(data.in,2,sd)
        data.in = sweep(data.in,2,datamean,'-')
        data.in = sweep(data.in,2,datasd,'/')
      }
      if(input$transform_fun=='median') {
        data.in = sweep(data.in,2,apply(data.in,2,median), '-')
      }
    }
    if(input$transform_fun=='percentize') data.in=heatmaply::percentize(data.in)
  } 
  
        
  if(!is.null(input$tables_true_search_columns)) 
    data.in=data.in[activeRows(input$tables_true_search_columns,data.in),]
  if(input$colRngAuto){
    ColLimits=NULL 
  }else{
    ColLimits=c(input$colorRng_min, input$colorRng_max)
  }
  
  distfun_row = function(x) dist(x, method = input$distFun_row)
  distfun_col =  function(x) dist(x, method = input$distFun_col)
  
  hclustfun_row = function(x) hclust(x, method = input$hclustFun_row)
  hclustfun_col = function(x) hclust(x, method = input$hclustFun_col)
  p <- heatmaply(data.in,
            main = input$main,xlab = input$xlab,ylab = input$ylab,
            row_text_angle = input$row_text_angle,
            column_text_angle = input$column_text_angle,
            dendrogram = input$dendrogram,
            branches_lwd = input$branches_lwd,
            seriate = input$seriation,
            colors=eval(parse(text=paste0(input$pal,'(',input$ncol,')'))),
            distfun_row =  distfun_row,
            hclustfun_row = hclustfun_row,
            distfun_col = distfun_col,
            hclustfun_col = hclustfun_col,
            k_col = input$c, 
            k_row = input$r,
            limits = ColLimits) %>% 
    layout(margin = list(l = input$l, b = input$b, r='0px'))
    
  p$elementId <- NULL
  
  p
})

#Render Plot ----
observeEvent(input$data,{
output$heatout <- renderPlotly({
  if(!is.null(input$data))
    interactiveHeatmap()
})
})



#Clone Heatmap ----
observeEvent({interactiveHeatmap()},{
  h<-interactiveHeatmap()
  
  l<-list(main = input$main,xlab = input$xlab,ylab = input$ylab,
          row_text_angle = input$row_text_angle,
          column_text_angle = input$column_text_angle,
          dendrogram = input$dendrogram,
          branches_lwd = input$branches_lwd,
          seriate = input$seriation,
          colors=paste0(input$pal,'(',input$ncol,')'),
          distfun_row =  input$distFun_row,
          hclustfun_row = input$hclustFun_row,
          distfun_col = input$distFun_col,
          hclustfun_col = input$hclustFun_col,
          k_col = input$c, 
          k_row = input$r,
          limits = paste(c(input$colorRng_min, input$colorRng_max),collapse=',')
  )
  
  #l=l[!l=='']
  l=data.frame(Parameter=names(l),Value=do.call('rbind',l),row.names = 1,stringsAsFactors = F)
  l[which(l$Value==''),2]='NULL'
  paramTbl=print(xtable::xtable(l),type = 'html',include.rownames=FALSE,print.results = F,html.table.attributes = c('border=0'))
  
  
  h$width='100%'
  h$height='800px'
  s<-tags$div(style="position: relative; bottom: 5px;",
              HTML(paramTbl),
                tags$em('This heatmap visualization was created using',
                  tags$a(href="https://github.com/yonicd/shinyHeatmaply/",target="_blank",'shinyHeatmaply'),
                  Sys.time()
                        )
              )
  
  output$downloadData <- downloadHandler(
    filename = function() {
      paste("heatmaply-", gsub(' ','_',Sys.time()), ".html", sep="")
    },
    content = function(file) {
      libdir <- paste(tools::file_path_sans_ext(basename(file)),"_files", sep = "")

      htmltools::save_html(htmltools::browsable(htmltools::tagList(h,s)),file=file,libdir = libdir)
      if (!htmlwidgets:::pandoc_available()) {
          stop("Saving a widget with selfcontained = TRUE requires pandoc. For details see:\n", 
          "https://github.com/rstudio/rmarkdown/blob/master/PANDOC.md")
      }
      
      htmlwidgets:::pandoc_self_contained_html(file, file)
      
      unlink(libdir, recursive = TRUE)
    }
  )
  
})


 #################################   FAQ   #####################################

  output$evaluation <- renderText({
    if(input$good |input$bad |input$fine)
      text = "Thanks for your evaluation of our application. Your reply is important to us."
  })
  
  output$ggAperturaTiendas=renderUI({
    HTML('<iframe src="http://127.0.0.1:5353/" style="width:100%;height:900px">')
    
  })

#End of Code ----
})



