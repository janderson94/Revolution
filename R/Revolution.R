#Packages that are neccesary for life
#Input_df is Ngenes by X Variables
#perm_df is nGenes x Ypermutations by 1 variable
#Pvals_col_name is the name of the column in Input_df with pvalues for empirical FDR correction
#Name is output file name FDR_"name"

perm.fdr=function(input_df,perm_df,Pvals_col_name,name){
  pvals_index=which(colnames(input_df)==Pvals_col_name)
  ro<-input_df[order(input_df[,pvals_index]),]
  p_obs <- data.frame(pvalue=ro[,pvals_index])
  p_vector<-matrix(as.matrix(perm_df),ncol=1)
  p_vector=data.frame(p_vector[order(p_vector)])

  F<-p_obs[,1]
  F_o<-p_obs[,1]
  pi_hat<-p_obs[,1]

  j=1
  observed<-length(p_obs[,1])
  randoms<-length(p_vector[,1])

  for(i in 1:observed)
  {
    repeat
    {
      if((p_vector[j,1]<p_obs[i,1])&j<randoms){j<-j+1}else{break}
    }
    F[i]=i/observed
    F_o[i]=(j-1)/randoms
    if(F_o[i]<1){pi_hat[i]=(1-F[i])/(1-F_o[i])}else{pi_hat[i]=1}
  }
  tabla <-data.frame(pi_hat,pval=p_obs[,1])

  tabla[1,]=c(1,0)
  last_percentile_average=mean(tabla$pi_hat[as.integer(min((length(tabla[,1])*0.99),(nrow(tabla)-1)):length(tabla[,1]))])
  tabla[nrow(tabla),]=c(last_percentile_average,1)
  constraint_matrix=as.matrix(data.frame(c(0,2),c(0,1),c(1,0)))
  f_hat<-suppressWarnings(cobs::cobs(tabla$pval,tabla$pi_hat,constraint="convex",pointwise=constraint_matrix,maxiter=1000,print.warn=FALSE,print.mesg=FALSE))

  f_hat_serie=f_hat$fitted
  pi_o=f_hat_serie[length(f_hat_serie)]
  pi_o=min(pi_o,1)
  pi_o=max(pi_o,0)

  Fdr_ST_perm=pi_o*F_o/F

  for(i in 1:length(p_obs[,1]))
  {
    Fdr_ST_perm[i]=pi_o*F_o[i]/F[i]
    if(i>1)
    {
      for(j in 1:(i-1))
      {
        if(Fdr_ST_perm[i-j]>Fdr_ST_perm[i]){Fdr_ST_perm[i-j]=Fdr_ST_perm[i]}else{break}
      }
    }
    if(Fdr_ST_perm[i]>1)  Fdr_ST_perm[i]=1
  }

  fdrs_df <-data.frame(ro,q_ST_perm=Fdr_ST_perm)
  rownames(fdrs_df)=rownames(ro)
  colnames(fdrs_df)[ncol(fdrs_df)]=paste0("fdr_",name)

  return(fdrs_df)
}






#Themes for ggplot2
library(ggplot2)
new_theme<-theme(axis.text=element_text(size=10)  , axis.title=element_text(size=12,face="bold") , legend.text=element_text(size=16),plot.title=element_text(size=18,hjust=.5),panel.grid.major = element_blank(),axis.line = element_line(colour = "black"),
                 panel.grid.minor = element_blank(),
                 panel.border = element_blank(),
                 panel.background = element_blank(),
                 legend.key = element_blank(),
                 legend.title = element_blank())




#Mode function
mode<-function(x){
  tab<-table(x)
  return(names(tab[tab==max(tab)]))
}







#Parallel age predictions
library(glmnet)
library(parallel)

#PARALLEL CODE FOR RUNNING IT ON YOUR COMPUTER
#Better to just run it on HARDAC if you are using many samples, or are building an epigenetic model, but:

#First define Z-score function if using GE data
get_z<-function(x){
  mean_x<-mean(x)
  sd_x<-sd(x)
  z<-(x-mean_x)/sd_x
  return(z)
}

parallel_predict<-function(age_vector,counts_file,sample_type,no_cores=NULL,alphas){

  # Calculate the number of cores if not provided
  if(!is.null(no_cores)){
    no_cores<-no_cores
  } else {no_cores <- detectCores() - 1}

N<-length(age_vector)
# Initiate cluster
cl <- makeCluster(no_cores)
clusterExport(cl=cl, varlist=ls(),envir = environment())
clusterEvalQ(cl,library(glmnet))
predicted<-matrix(NA,nrow=N,ncol=length(alphas))
colnames(predicted)<-alphas

par_out<-parLapply(cl,X = alphas,
                   function(z){
                     for (i in 1:N){
                       #Quantile normalize all samples and remove test subject
                       #If GE, log2 CPM +1 normalize then QQ normalize
                       if(sample_type=="GE"){
                       norm_counts<-as.matrix(apply(counts_file,2,function(y){return(qqnorm(log2(y/sum(y)+1),plot=F)$x)}))
                       }
                       #If GE, QQ normalize
                        else if (sample_type=="epi"){
                          norm_counts<-as.matrix(apply(counts_file,2,function(y){return(qqnorm(y,plot=F)$x)}))
                        }
                       else {
                         return("invalid sample_type")
                       }

                       norm_train<-norm_counts[,-i]
                       norm_test<-norm_counts[,i]

                       #Quantile or Z-score normalize by rows (by gene or by site)
                       if(sample_type=="epi"){
                       trainreads_norm<-as.matrix(apply(norm_train,1,function(y){return(qqnorm(y,plot=F)$x)}))}

                       else if(sample_type=="GE"){
                       trainreads_norm<-as.matrix(apply(norm_train,1,function(x){return(get_z(x))}))}

                       else {
                         return("invalid sample_type")
                       }

                       trainage<-age_vector[-i]
                       testage<-age_vector[i]

                       #For each gene in the test set, define the quantile relative to the training set
                       #Following this, use this quantile to give it a value from a normal distribution

                       if(sample_type=="epi"){
                       testreads_normalized<-norm_test
                       for (d in 1:dim(norm_counts)[1]){
                         a<-ecdf(norm_train[d,])
                         probs<-a(norm_test[d])
                         probs[probs==1]<-.99
                         probs[probs==0]<-.01
                         testreads_normalized[d]<-qnorm(probs)
                       }
                       }

                       else if (sample_type=="GE"){
                       #Z-score normalize each row
                       testreads_normalized<-norm_test
                       for (d in 1:dim(norm_counts)[1]){
                         mean_d<-mean(norm_train[d,])
                         sd_d<-sd(norm_train[d,])
                         testreads_normalized[d]<-(norm_test[d]-mean_d)/sd_d
                       }
                       }

                       else {
                         return("invalid sample_type")
                       }

                       #Always using N-fold internal CV -- Not asking GLMNET to standardize
                       model<-cv.glmnet(trainreads_norm,trainage,nfolds=(N-1),alpha=z,standardize=F)
                       predicted[i,colnames(predicted)==z]<-predict(model,newx=t(testreads_normalized),s="lambda.min")
                     }
                     return(predicted)
                     })
stopCluster(cl)
return(par_out)
}


















