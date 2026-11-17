# Verify named regional point estimates after excluding unlabelled rows in a separate copy.
# Knots are fixed at saved native selections. This is NOT a new NCI selection or CI run.
d<-read.delim('analysis/input/region_CKD_named_cohorts_only.txt',check.names=FALSE)
a<-read.delim('analysis/source/joinpoint_native/NCI_region_CKD_2013_2024.Export.APC.txt',check.names=FALSE)
b<-read.delim('analysis/source/joinpoint_native/NCI_region_CKD_2013_2024.Export.AAPC.txt',check.names=FALSE)
groups<-c('Northeast','Midwest','South','West');out<-data.frame()
for(i in seq_along(groups)) {
 x<-d[d$Region==groups[i],];x$t<-x$Year-2013
 seg<-a[a$Region==i-1 & is.finite(a$APC),];knots<-sort(unique(seg[['Segment Start']][seg[['Segment Start']]>2013]))
 rhs<-'t'
 if(length(knots))for(j in seq_along(knots)) {n<-paste0('hinge',j);x[[n]]<-pmax(0,x$Year-knots[j]);rhs<-paste(rhs,n,sep='+')}
 f<-lm(as.formula(paste('log(Rate)~',rhs)),data=x,weights=(Rate/x[['Standard Error']])^2)
 slope<-cumsum(coef(f)[-1]);width<-diff(c(2013,knots,2024));point<-100*expm1(sum(slope*width)/11)
 z<-b[b$Region==i-1 & b[['Start Obs']]==2013 & b[['End Obs']]==2024,]
 out<-rbind(out,data.frame(region=groups[i],n=nrow(x),R_fixed_knots_AAPC=point,native_AAPC=z$AAPC,agrees_4dp=abs(point-z$AAPC)<.000051))
}
stopifnot(nrow(out)==4,all(out$n==12),all(out$agrees_4dp))
write.csv(out,'analysis/results/region_review_copy_point_checks.csv',row.names=FALSE)
print(out)
