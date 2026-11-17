# Revision analysis, 23 October 2026. Base R only. No patient data are simulated.
# Run from the project folder: Rscript analysis/reproduce_analysis.R
# Optional verified aggregate inputs go in analysis/input (see README).
options(stringsAsFactors=FALSE, warn=1)
dir.create('analysis/results',recursive=TRUE,showWarnings=FALSE)
dir.create('analysis/input',recursive=TRUE,showWarnings=FALSE)
sink('analysis/results/run_log.txt',split=TRUE)
cat('Rerun of available published aggregate counts\n')
cat('Sources are included under analysis/source.\n')
status <- data.frame(analysis=character(),status=character(),detail=character())
record <- function(a,s,d) {status[nrow(status)+1,] <<- c(a,s,d);cat(a,':',s,'-',d,'\n')}
savecsv <- function(x,n) write.csv(x,file.path('analysis/results',n),row.names=FALSE)
# Values transcribed from manuscript Table 1, Table 2 and Supplementary Table S5.
endpoints <- data.frame(year=c(1999,2013,2024),lung=c(152063,156178,128665),ckd=c(508,1110,1826),
 provenance=c('Main Table 1 and Table 2 Panel B','Recovered Desktop C34 overall export; not inferred from percentage','Main Table 1 and Table 2 Panel A'))
endpoints$crude_percent <- 100*endpoints$ckd/endpoints$lung
savecsv(endpoints,'national_endpoint_checks.csv')
age <- data.frame(age=c('45-54','55-64','65-74','75-84','85+'),lung=c(4049,21201,43871,41602,17086),ckd=c(20,153,462,768,417),reported_percent=c(.494,.722,1.053,1.846,2.441))
age$recomputed_percent <- 100*age$ckd/age$lung
age$matches_3dp <- round(age$recomputed_percent,3)==age$reported_percent
stopifnot(all(age$matches_3dp),all(age$ckd<=age$lung))
savecsv(age,'age_2024_checks.csv')
record('Published endpoint arithmetic','RUN','Five age-specific proportions agree at three decimals; national 1999 and 2024 proportions agree.')
issues <- data.frame(id=c('D01','D02','D03','D05','D06'),issue=c('National <55 counts recovered from archive; standalone originals not located','Rounded endpoints subtract to 0.702 pp; full contrast is 0.702912775 pp','Printed age weights sum to 1.00001','Residual r and P differ from manuscript','Native APC upper CI is 5.9370 percent rather than 5.9570 percent'),action=c('Resolved in round two: all 120 counts match new official CDC query responses','Resolved as rounding','Resolved: exact weights computed from counts','Verified correction: R and Python agree at r 0.589912 and P 0.056093','Verified native transcription correction: upper limit 5.9370'),data_changed=FALSE)
savecsv(issues,'data_review_register.csv')
savecsv(data.frame(year=2024,all_lung=128665,age45plus_lung=sum(age$lung),remainder_under45_or_unknown=128665-sum(age$lung),ckd=1826,age45plus_ckd=sum(age$ckd),remainder_under45_or_unknown_ckd=1826-sum(age$ckd)),'age_count_reconciliation.csv')
record('Data changes','NONE','No source counts changed. Proposed result corrections D05/D06 and fresh sensitivities are labeled separately.')
grad <- function(f,newdata,weights) {
 X <- model.matrix(delete.response(terms(f)),newdata,contrasts.arg=f$contrasts,xlev=f$xlevels)
 p <- plogis(drop(X%*%coef(f)))
 list(p=sum(weights*p),g=colSums(X*as.vector(weights*p*(1-p))))
}
hac <- function(f,year,lag=1) {
 X <- model.matrix(f); y <- f$y; n <- f$prior.weights; p <- fitted(f)
 # glm with two-column response stores y as success fraction and trial totals in prior.weights.
 scores <- X*as.vector(n*(y-p))
 S <- rowsum(scores,year,reorder=TRUE); B <- solve(crossprod(X,X*as.vector(n*p*(1-p))))
 M <- crossprod(S)
 for(k in seq_len(lag)) {G<-crossprod(S[(k+1):nrow(S),,drop=FALSE],S[1:(nrow(S)-k),,drop=FALSE]);M<-M+(1-k/(lag+1))*(G+t(G))}
 B%*%M%*%B # No dispersion multiplier: score meat already reflects observed variability.
}
or_ci <- function(b,v) exp(c(OR=b,lower=b-1.96*sqrt(v),upper=b+1.96*sqrt(v)))
infile <- 'analysis/input/national_age_counts.csv'
if(file.exists(infile)) {
 d<-read.csv(infile);stopifnot(all(c('year','age','lung','ckd')%in%names(d)))
 stopifnot(all(complete.cases(d)),!anyDuplicated(d[c('year','age')]),all(d$lung>0),all(d$ckd>=0),all(d$ckd<=d$lung),all(d$lung==floor(d$lung)),all(d$ckd==floor(d$ckd)))
 d$age<-factor(d$age,levels=c('<55','55-64','65-74','75-84','85+'));stopifnot(!anyNA(d$age),setequal(d$year,2013:2024),nrow(d)==60)
 d<-d[order(d$year,d$age),]; d$year_c<-d$year-2013
 f<-glm(cbind(ckd,lung-ckd)~year_c+age,data=d,family=quasibinomial())
 stopifnot(nobs(f)==60,f$rank==6,df.residual(f)==54)
 phi<-summary(f)$dispersion; V<-vcov(f); j<-which(names(coef(f))=='year_c')
 savecsv(data.frame(method=c('Quasi binomial','HAC lag 1','HAC lag 2'),rbind(or_ci(coef(f)[j],V[j,j]),or_ci(coef(f)[j],hac(f,d$year,1)[j,j]),or_ci(coef(f)[j],hac(f,d$year,2)[j,j]))),'national_annual_OR.csv')
 w<-tapply(d$lung,d$age,sum);w<-w/sum(w)
 savecsv(data.frame(age=names(w),weight=as.vector(w)),'national_exact_weights.csv')
 standards<-list(pooled=w,year2013=with(subset(d,year==2013),setNames(lung/sum(lung),age)),year2024=with(subset(d,year==2024),setNames(lung/sum(lung),age)))
 pred<-data.frame();contrast<-data.frame()
 for(s in names(standards)) {
  ww<-standards[[s]][levels(d$age)]; pp<-list()
  for(y in 2013:2024) {nd<-data.frame(year_c=y-2013,age=factor(levels(d$age),levels=levels(d$age)));q<-grad(f,nd,ww);se<-sqrt(drop(t(q$g)%*%V%*%q$g));pp[[as.character(y)]]<-q;pred<-rbind(pred,data.frame(standard=s,year=y,percent=100*q$p,lower=100*(q$p-1.96*se),upper=100*(q$p+1.96*se)))}
  dd<-pp[['2024']]$p-pp[['2013']]$p;gg<-pp[['2024']]$g-pp[['2013']]$g;se<-sqrt(drop(t(gg)%*%V%*%gg));contrast<-rbind(contrast,data.frame(standard=s,change_pp=100*dd,lower=100*(dd-1.96*se),upper=100*(dd+1.96*se)))
 }
 savecsv(pred,'national_standardized_predictions.csv');savecsv(contrast,'national_standardized_contrasts.csv')
 fi<-glm(cbind(ckd,lung-ckd)~year_c*age,data=d,family=quasibinomial());ix<-grep(':',names(coef(fi)));bb<-coef(fi)[ix];W<-drop(t(bb)%*%solve(vcov(fi)[ix,ix],bb))
 savecsv(data.frame(Wald=W,df=length(ix),P=pchisq(W,length(ix),lower.tail=FALSE)),'national_age_interaction.csv')
 fc<-glm(cbind(ckd,lung-ckd)~factor(year)+age,data=d,family=quasibinomial())
 cp<-do.call(rbind,lapply(2013:2024,function(y){q<-grad(fc,data.frame(year=y,age=factor(levels(d$age),levels=levels(d$age))),w);data.frame(year=y,percent=100*q$p)}));savecsv(cp,'categorical_year_predictions.csv')
 pear<-rowsum(residuals(f,type='pearson'),d$year);ac<-cor.test(pear[-nrow(pear),1],pear[-1,1])
 capture.output(summary(f),file='analysis/results/national_model_summary.txt');savecsv(data.frame(n=nobs(f),parameters=f$rank,residual_df=df.residual(f),dispersion=phi,residual_lag1=unname(ac$estimate),P=ac$p.value),'national_diagnostics.csv')
 record('National grouped models','RUN','Primary quasi binomial, age interaction, alternative weights, categorical year and HAC lags 1 and 2.')
} else record('National grouped models','NOT RUN','Missing national_age_counts.csv; printed endpoints cannot identify the 60 annual age cells.')
rf<-'analysis/input/rate_series.csv'
if(file.exists(rf)) {
 d<-read.csv(rf);stopifnot(all(c('series','year','rate','se','reliable')%in%names(d)))
 stopifnot(!anyDuplicated(d[c('series','year')]))
 result<-data.frame()
 for(s in unique(d$series)) {
  z<-subset(d,series==s & year>=2013 & year<=2024 & reliable & is.finite(rate) & rate>0 & is.finite(se) & se>0)
  if(nrow(z)<6) {record(paste('Rate sensitivity',s),'NOT RUN','Fewer than six usable annual rates');next}
  z$hinge<-pmax(0,z$year-2020);f<-lm(log(rate)~year,data=z);g<-lm(log(rate)~year+hinge,data=z,weights=(rate/se)^2)
  ci<-confint(f)['year',];result<-rbind(result,data.frame(series=s,n=nrow(z),OLS_APC=100*expm1(coef(f)['year']),lower=100*expm1(ci[1]),upper=100*expm1(ci[2]),fixed2020_pre_APC=100*expm1(coef(g)['year']),fixed2020_post_APC=100*expm1(coef(g)['year']+coef(g)['hinge'])))
 }
 savecsv(result,'rate_sensitivity.csv')
 record('Rate sensitivity','RUN','Exploratory OLS and fixed 2020 continuous hinge WLS; neither reproduces NCI permutation selection or empirical quantile intervals.')
} else record('Rate sensitivity','NOT RUN','Missing annual rate-and-SE input; no interpolation of figure trajectories.')
record('NCI native execution outside R','SEPARATE EXECUTED WORKFLOW','Three native runs were executed in Joinpoint 5.4.0.0: overall all-lung full-period AAPC and two cleaned regional inputs. See native exports and continued_native_validation.json; R does not execute NCI selection or empirical-quantile algorithms.')
record('Comparator and renal-code revision analyses','SEPARATE EXECUTED WORKFLOW','See round2_sensitivities.R and revision_sensitivity_models.csv. New official national data were retrieved; the original primary inputs were not overwritten.')
uf<-'analysis/input/urbanicity_age_counts.csv'
if(file.exists(uf)) {
 u<-read.csv(uf);excluded<-subset(u,urban=='Not Available');savecsv(excluded,'urbanicity_excluded_not_available.csv')
 u<-subset(u,urban!='Not Available');u$lung<-as.numeric(u$lung);u$ckd<-as.numeric(u$ckd)
 u$urban<-factor(u$urban,levels=c('Large Central Metro','Large Fringe Metro','Medium Metro','Small Metro','Micropolitan (Nonmetro)','NonCore (Nonmetro)'))
 u$age<-factor(u$age,levels=c('55-64','65-74','75-84','85+'));u$year_c<-u$year-2018
 stopifnot(nrow(u)==168,!anyNA(u),!anyDuplicated(u[c('year','urban','age')]),all(u$ckd<=u$lung))
 f<-glm(cbind(ckd,lung-ckd)~year_c*urban+age,data=u,family=quasibinomial());V<-vcov(f);ix<-grep(':',names(coef(f)));b<-coef(f)[ix];W<-drop(t(b)%*%solve(V[ix,ix],b))
 w<-tapply(u$lung,u$age,sum);w<-w/sum(w);up<-data.frame();ors<-data.frame()
 for(uu in levels(u$urban)) {
  for(y in 2018:2024){nd<-data.frame(year_c=y-2018,urban=factor(uu,levels=levels(u$urban)),age=factor(levels(u$age),levels=levels(u$age)));q<-grad(f,nd,w);se<-sqrt(drop(t(q$g)%*%V%*%q$g));up<-rbind(up,data.frame(urban=uu,year=y,percent=100*q$p,lower=100*(q$p-1.96*se),upper=100*(q$p+1.96*se)))}
  nd<-data.frame(year_c=c(0,1),urban=factor(uu,levels=levels(u$urban)),age=factor('55-64',levels=levels(u$age)))
  X<-model.matrix(delete.response(terms(f)),nd,contrasts.arg=f$contrasts);v<-X[2,]-X[1,];bb<-sum(v*coef(f));vv<-drop(t(v)%*%V%*%v)
  ors<-rbind(ors,data.frame(urban=uu,t(or_ci(bb,vv)),P=2*pnorm(-abs(bb/sqrt(vv)))))
 }
 savecsv(up,'urbanicity_predictions.csv');savecsv(ors,'urbanicity_OR.csv');savecsv(data.frame(n=nobs(f),parameters=f$rank,residual_df=df.residual(f),dispersion=summary(f)$dispersion,Wald=W,df=length(ix),P=pchisq(W,length(ix),lower.tail=FALSE),lung=sum(u$lung),ckd=sum(u$ckd)),'urbanicity_diagnostics.csv')
 savecsv(data.frame(age=names(w),weight=as.vector(w)),'urbanicity_weights.csv')
 png('analysis/results/urbanicity_trajectories.png',width=1800,height=1100,res=180);par(mar=c(4.2,4.5,3,1));cols<-c('#264653','#2a9d8f','#b08b00','#e76f51','#6a4c93','#457b9d');plot(range(up$year),range(up$percent),type='n',xlab='Year',ylab='Age-standardized co-recording (%)',main='Conditional CKD co-recording among lung-cancer deaths aged 55 years or older')
 for(j in seq_along(levels(u$urban))){z<-subset(up,urban==levels(u$urban)[j]);lines(z$year,z$percent,col=cols[j],lwd=2)}
 legend('topleft',legend=levels(u$urban),col=cols,lwd=2,bty='n',cex=.8);mtext(sprintf('Global year by urbanicity interaction P = %.3f; curves do not establish differential trends',pchisq(W,length(ix),lower.tail=FALSE)),side=3,line=.2,cex=.65);dev.off()
 record('Urbanicity model','RUN','168 original WONDER cells; quasi binomial model, pooled standardization and global interaction reproduced in R.')
} else record('Urbanicity model','NOT RUN','Missing source grouped cells.')
record('Race and state standardized co-recording','NOT RUN','Year by age by race/state cells absent; interpretation restricted to crude exploratory descriptions.')
savecsv(status,'analysis_status.csv')
capture.output(sessionInfo(),file='analysis/results/sessionInfo.txt')
sink()
