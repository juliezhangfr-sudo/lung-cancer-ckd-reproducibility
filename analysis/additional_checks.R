# Additional revision analyses using located records and archived native exports.
options(stringsAsFactors=FALSE)
source('analysis/reproduce_analysis.R')
sink('analysis/results/additional_log.txt',split=TRUE)
ar<-read.csv('analysis/input/archived_records.csv');ar$deaths<-suppressWarnings(as.numeric(ar$deaths));ar$rate<-suppressWarnings(as.numeric(ar$rate));ar$se<-suppressWarnings(as.numeric(ar$se));ar$population<-suppressWarnings(as.numeric(ar$population))
national<-subset(ar,source%in%c('Lung_all_1999_2020.csv','Lung_CKD_1999_2020.csv','Lung_all_2018_2024.csv','Lung_CKD_2018_2024.csv'))
national<-subset(national,!(grepl('2018_2024',source)&year<=2020));national$kind<-ifelse(grepl('Lung_CKD',national$source),'ckd','lung')
n<-reshape(national[c('year','kind','deaths','rate','se')],idvar='year',timevar='kind',direction='wide');n<-n[order(n$year),];n$crude_percent<-100*n$deaths.ckd/n$deaths.lung
savecsv(n,'national_annual_verified_records.csv')
d<-read.csv('analysis/input/national_age_counts.csv');z<-aggregate(cbind(lung,ckd)~year,d,sum);z<-merge(z,n[c('year','deaths.lung','deaths.ckd')]);z$unknown_age_lung<-z$deaths.lung-z$lung;z$unknown_age_ckd<-z$deaths.ckd-z$ckd;savecsv(z,'unknown_age_reconciliation.csv')
cat('Unknown-age exclusions',sum(z$unknown_age_lung),sum(z$unknown_age_ckd),'\n')
# Verify age-model point estimates against native AAPC exports, not NCI interval reruns.
ai<-read.csv('analysis/source/age_audit/CKD_2013_2024/Lung_CKD_Age_2013_2024_Joinpoint_FINAL.csv')
ei<-read.delim('analysis/source/age_audit/CKD_2013_2024/exports/Lung_CKD_Age_2013_2024_Joinpoint_EXPORT_FINAL.Export.AAPC.txt',check.names=FALSE)
aa<-unique(ai$Age_Group);ao<-data.frame()
for(i in seq_along(aa)) {q<-subset(ai,Age_Group==aa[i]);f<-lm(log(Rate)~Year,data=q,weights=(Rate/Standard_Error)^2);b<-coef(f)['Year'];orig<-ei$AAPC[i];ao<-rbind(ao,data.frame(age=aa[i],n=nrow(q),R_WLS_APC=100*expm1(b),native_AAPC=orig,agrees_4dp=abs(100*expm1(b)-orig)<.000051,year2015_present=2015%in%q$Year))}
stopifnot(all(ao$agrees_4dp));savecsv(ao,'age_native_point_checks.csv')
q<-subset(ai,Age_Group=='45-54 years');q<-rbind(q,data.frame(Age_Group='45-54 years',Year=2015,Rate=17/43188161*100000,Standard_Error=.009546842,Deaths=17,Source='Archived 17-death export; unreliable rate derived from disclosed count/population for sensitivity only'))
f<-lm(log(Rate)~Year,data=q,weights=(Rate/Standard_Error)^2)
savecsv(data.frame(analysis=c('Exclude flagged 2015 point','Include count-derived unreliable 2015 point for sensitivity only'),n=c(11,12),APC=c(ao$R_WLS_APC[1],100*expm1(coef(f)['Year']))),'age45_unreliable_sensitivity.csv')
# Original within-era OLS sensitivity results.
race<-data.frame()
for(g in c('NHWhite','NHBlack','Hispanic','Other'))for(era in c('1999_2020','2018_2024')){
 q<-subset(ar,source==paste0('Lung_CKD_',g,'_',era,'.csv')&is.finite(rate)&rate>0);f<-lm(log(rate)~year,data=q);ci<-confint(f)['year',];race<-rbind(race,data.frame(group=g,era=era,n=nrow(q),APC=100*expm1(coef(f)['year']),lower=100*expm1(ci[1]),upper=100*expm1(ci[2]),P=summary(f)$coefficients['year','Pr(>|t|)']))
}
savecsv(race,'race_era_OLS.csv')
# Overall pre-transition OLS sensitivity.
q<-subset(n,year<=2010);f<-lm(log(rate.ckd)~year,data=q);ci<-confint(f)['year',];savecsv(data.frame(APC=100*expm1(coef(f)['year']),lower=100*expm1(ci[1]),upper=100*expm1(ci[2]),P=summary(f)$coefficients['year','Pr(>|t|)']),'pretransition_OLS.csv')
# Recreate Figure 3 from the reproduced national model, not graph digitization.
p<-read.csv('analysis/results/national_standardized_predictions.csv');p<-subset(p,standard=='pooled');catp<-read.csv('analysis/results/categorical_year_predictions.csv');raw<-subset(n,year>=2013)
png('analysis/results/figure3_reproduced.png',width=2400,height=1050,res=200);par(mfrow=c(1,2),mar=c(4.2,4.5,3.5,1),cex=.9)
plot(p$year,p$percent,type='n',ylim=range(p$lower,p$upper,raw$crude_percent),xlab='Year',ylab='Percentage of lung-cancer deaths',main='A  Observed crude and marginal standardized proportions',cex.main=.9)
polygon(c(p$year,rev(p$year)),c(p$lower,rev(p$upper)),col='#dce9f0',border=NA);lines(p$year,p$percent,col='#23658c',lwd=2);lines(raw$year,raw$crude_percent,col='#b45a24',lwd=2,lty=2);legend('topleft',c('Marginal age standardized','Observed crude'),col=c('#23658c','#b45a24'),lty=c(1,2),lwd=2,bty='n',cex=.8)
plot(p$year,p$percent,type='l',ylim=range(p$percent,catp$percent),col='#23658c',lwd=2,xlab='Year',ylab='Percentage of lung-cancer deaths',main='B  Marginal standardized model-form sensitivity',cex.main=.9);lines(catp$year,catp$percent,col='#51713a',lwd=2,lty=2);legend('topleft',c('Continuous year','Categorical year'),col=c('#23658c','#51713a'),lty=c(1,2),lwd=2,bty='n',cex=.8);dev.off()
cat('All five CKD age point estimates agree with native exports. New sensitivity results are separate from NCI selection and intervals.\n')
sink()
