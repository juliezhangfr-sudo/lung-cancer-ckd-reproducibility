# Run from project directory after prepare_round2_inputs.py.
# Public aggregate counts only; no suppression reconstruction or simulated data.
options(stringsAsFactors=FALSE)
args<-commandArgs(trailingOnly=TRUE);root<-if(length(args)) args[1] else if(dir.exists('analysis/input')) 'analysis' else 'round2_revision/analysis'
dir.create(file.path(root,'results'),recursive=TRUE,showWarnings=FALSE)
sink(file.path(root,'results/round2_R_log.txt'),split=TRUE)
d<-read.csv(file.path(root,'input/revision_sensitivity_counts.csv'))
savecsv<-function(x,n)write.csv(x,file.path(root,'results',n),row.names=FALSE)
levels_age<-c('55-64','65-74','75-84','85+')
den<-subset(d,outcome=='denominator'&age%in%levels_age)
names(den)[names(den)=='count']<-'lung'
result<-data.frame(); predictions<-data.frame();quality<-data.frame()
for(out in c('N18','diabetes','broad_renal','stage5')) {
 q<-subset(d,outcome==out&age%in%levels_age)
 quality<-rbind(quality,data.frame(outcome=out,returned=nrow(q),non_numeric=sum(!is.finite(q$count)),min_deaths=if(any(is.finite(q$count)))min(q$count,na.rm=TRUE) else NA,below20=sum(q$count<20,na.rm=TRUE)))
 if(nrow(q)!=48||nrow(den)!=48||any(!is.finite(q$count))||any(!is.finite(den$lung))) {cat(out,': full 48-cell analysis unavailable; no suppressed cells imputed\n');next}
 x<-merge(q,den[c('year','age','lung')],by=c('year','age'));stopifnot(nrow(x)==48,all(x$count>=0),all(x$count<=x$lung))
 x$age<-factor(x$age,levels=levels_age);x$year_c<-x$year-2013
 f<-glm(cbind(count,lung-count)~year_c+age,data=x,family=binomial(),control=glm.control(epsilon=1e-10,maxit=200))
 stopifnot(f$converged,df.residual(f)==43)
 phi<-sum(residuals(f,type='pearson')^2)/df.residual(f);V<-vcov(f)*phi;b<-coef(f)['year_c'];se<-sqrt(V['year_c','year_c'])
 # Independent family implementation verifies coefficients and dispersion-scaled covariance.
 fq<-glm(cbind(count,lung-count)~year_c+age,data=x,family=quasibinomial(),control=glm.control(epsilon=1e-10,maxit=200))
 stopifnot(max(abs(coef(f)-coef(fq)))<1e-8,max(abs(V-vcov(fq)))<1e-8)
 # Year-aggregated score HAC with lag 1, no finite sample correction.
 X<-model.matrix(f);p<-fitted(f);U<-rowsum(X*as.vector(x$count-x$lung*p),x$year)
 meat<-crossprod(U)+.5*(crossprod(U[-1,,drop=FALSE],U[-nrow(U),,drop=FALSE])+crossprod(U[-nrow(U),,drop=FALSE],U[-1,,drop=FALSE]))
 bread<-solve(crossprod(X,X*as.vector(x$lung*p*(1-p))));H<-bread%*%meat%*%bread
 w<-as.numeric(tapply(x$lung,x$age,sum));w<-w/sum(w);grads<-list();ep<-list()
 for(y in 2013:2024) {
  nd<-data.frame(year_c=y-2013,age=factor(levels_age,levels=levels_age));A<-model.matrix(~year_c+age,nd);pr<-as.numeric(plogis(A%*%coef(f)));g<-colSums(A*as.vector(w*pr*(1-pr)));s<-sqrt(as.numeric(t(g)%*%V%*%g));v<-sum(w*pr)
  ep[[as.character(y)]]<-v;grads[[as.character(y)]]<-g
  predictions<-rbind(predictions,data.frame(outcome=out,year=y,percent=100*v,lower=100*(v-1.96*s),upper=100*(v+1.96*s)))
 }
 delta<-ep[['2024']]-ep[['2013']];g<-grads[['2024']]-grads[['2013']];ds<-sqrt(as.numeric(t(g)%*%V%*%g))
 # Finite-difference gradient independently checks the marginal contrast calculation.
 pred<-function(beta,y) {nd<-data.frame(year_c=y-2013,age=factor(levels_age,levels=levels_age));sum(w*plogis(model.matrix(~year_c+age,nd)%*%beta))}
 fd<-sapply(seq_along(coef(f)),function(j){u<-v<-coef(f);u[j]<-u[j]+1e-6;v[j]<-v[j]-1e-6;((pred(u,2024)-pred(u,2013))-(pred(v,2024)-pred(v,2013)))/2e-6})
 stopifnot(max(abs(fd-g))<1e-7)
 result<-rbind(result,data.frame(outcome=out,n=48,parameters=5,residual_df=43,phi=phi,events=sum(x$count),denominator=sum(x$lung),OR=exp(b),lower=exp(b-1.96*se),upper=exp(b+1.96*se),P=2*pnorm(-abs(b/se)),HAC_lower=exp(b-1.96*sqrt(H['year_c','year_c'])),HAC_upper=exp(b+1.96*sqrt(H['year_c','year_c'])),start_percent=100*ep[['2013']],end_percent=100*ep[['2024']],change_pp=100*delta,change_lower=100*(delta-1.96*ds),change_upper=100*(delta+1.96*ds)))
 cat(out,': OR',exp(b),'dispersion',phi,'events',sum(x$count),'\n')
}
savecsv(result,'revision_sensitivity_models.csv');savecsv(predictions,'revision_sensitivity_predictions.csv');savecsv(quality,'revision_suppression_audit.csv')
# Validate all newly retrieved national primary cells without overwriting primary inputs.
checks<-read.csv(file.path(root,'input/primary_count_comparison.csv'))
checks$agreement<-tolower(as.character(checks$agreement))=='true'
savecsv(data.frame(checked=nrow(checks),agreements=sum(checks$agreement),differences=sum(!checks$agreement)),'primary_count_validation.csv')
if(nrow(checks)==120)stopifnot(all(checks$agreement))
# Release overlap checked by outcome and age; never silently substitute vintages.
all<-read.csv(file.path(root,'input/new_cdc_counts_all_releases.csv'))
a<-subset(all,database=='D77'&year%in%2018:2020);b<-subset(all,database=='D157'&year%in%2018:2020)
ov<-merge(a,b,by=c('outcome','year','age'),suffixes=c('_bridged','_single'))
ov$difference<-ov$count_single-ov$count_bridged;savecsv(ov,'new_release_overlap.csv')
capture.output(sessionInfo(),file=file.path(root,'results/round2_sessionInfo.txt'))
cat('Primary source counts were not edited. New sensitivity inputs are separate files.\n');print(result);sink()
