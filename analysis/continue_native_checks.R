# Run from the repository root. Base R only.
# Saved native selections are fixed here: this does not rerun NCI selection or CIs.
options(stringsAsFactors=FALSE)
native <- 'analysis/source/joinpoint_native'
dir.create('analysis/results',showWarnings=FALSE)
sink('analysis/results/continued_R_log.txt',split=TRUE)
groups <- list(overall='Overall',sex=c('Female','Male'),
 race=c('Hispanic','NH Black','NH White','Other'),
 region=c('Northeast','Midwest','South','West'),urbanicity=c('Metro','Nonmetro'))
baseline <- read.csv('analysis/input/retained_AAPC_baseline.csv',check.names=FALSE)
point_checks <- apc_checks <- data.frame()
for (jpo in list.files(native,pattern='^NCI_.*[.]jpo$',full.names=TRUE)) {
 if(grepl('_age_',jpo))next
 stem <- sub('[.]jpo$','',jpo)
 kind <- names(groups)[vapply(names(groups),function(k)grepl(paste0('_',k,'_'),basename(jpo)),logical(1))]
 stopifnot(length(kind)==1)
 disease <- if(grepl('_CKD_',jpo))'CKD' else 'all'
 d <- read.delim(paste0(stem,'.txt'),check.names=FALSE)
 a <- read.delim(paste0(stem,'.Export.APC.txt'),check.names=FALSE)
 for (i in seq_along(groups[[kind]])) {
  group <- groups[[kind]][i]
  x <- if(kind=='overall')d else d[d[[1]]==group,]
  seg <- if(kind=='overall')a else a[a[[1]]==i-1,]
  seg <- seg[is.finite(seg$APC) & seg[["Segment End"]]>seg[["Segment Start"]],]
  start <- min(x$Year); end <- max(x$Year)
  knots <- sort(unique(seg[['Segment Start']][seg[['Segment Start']]>start]))
  x$t <- x$Year-start; rhs <- 't'
  if(length(knots))for(k in seq_along(knots)) {
   name <- paste0('hinge',k); x[[name]]<-pmax(0,x$Year-knots[k]);rhs<-paste(rhs,name,sep='+')
  }
  fit <- lm(as.formula(paste('log(Rate)~',rhs)),x,weights=(Rate/x[['Standard Error']])^2)
  slopes <- cumsum(coef(fit)[-1]); widths <- diff(c(start,knots,end))
  apcs <- 100*expm1(slopes); aapc <- 100*expm1(sum(widths*slopes)/(end-start))
  old <- baseline[baseline$Group==group & baseline$Disease==disease & baseline$Period==paste(start,end,sep='-'),]
  if(length(apcs)!=nrow(seg))stop(paste(basename(jpo),group,length(apcs),nrow(seg),paste(seg$APC,collapse=",")));stopifnot(nrow(old)==1)
  point_checks <- rbind(point_checks,data.frame(session=basename(jpo),group,disease,n=nrow(x),start,end,
   R_AAPC=aapc,retained_AAPC=old$AAPC,absolute_difference=abs(aapc-old$AAPC),agrees_4dp=abs(aapc-old$AAPC)<0.000051))
  apc_checks <- rbind(apc_checks,data.frame(session=basename(jpo),group,segment=seq_along(apcs),
   R_APC=unname(apcs),native_APC=seg$APC,absolute_difference=abs(apcs-seg$APC),agrees_4dp=abs(apcs-seg$APC)<0.000051))
 }
}
write.csv(point_checks,'analysis/results/continued_fixed_knot_AAPC_checks.csv',row.names=FALSE)
write.csv(apc_checks,'analysis/results/continued_fixed_knot_APC_checks.csv',row.names=FALSE)
provenance <- data.frame()
for(disease in c('CKD','alllung')) {
 years <- if(disease=='CKD')'2013_2024' else '1999_2024'
 f <- file.path(native,paste0('NCI_region_',disease,'_',years,'.txt'))
 d <- read.delim(f,check.names=FALSE)
 o <- read.delim(file.path(native,paste0('NCI_overall_',disease,'_',years,'.txt')),check.names=FALSE)
 ix <- which(d$Region=='')
 for(i in ix) {
  national <- o[o$Year==d$Year[i],]
  stopifnot(nrow(national)==1)
  rmatch <- d$Rate[i]==national$Rate; sematch <- d[['Standard Error']][i]==national[['Standard Error']]
  provenance<-rbind(provenance,data.frame(disease,original_file=basename(f),original_line=i+1,year=d$Year[i],
   unlabelled_rate=d$Rate[i],national_rate=national$Rate,unlabelled_SE=d[['Standard Error']][i],national_SE=national[['Standard Error']],
   exact_rate_match=rmatch,exact_SE_match=sematch,action='Exclude national total from separate named-region analysis copy; preserve original'))
 }
 clean <- d[d$Region %in% groups$region,]
 stopifnot(nrow(clean)==4*nrow(o),all(table(clean$Region)==nrow(o)))
 write.table(clean,file.path('analysis/input',paste0('region_',disease,'_named_cohorts_only.txt')),sep='\t',row.names=FALSE,quote=FALSE)
}
write.csv(provenance,'analysis/results/regional_total_rows_one_by_one.csv',row.names=FALSE)
print(point_checks);print(apc_checks)
cat('\nRegional total row exact matches:',sum(provenance$exact_rate_match & provenance$exact_SE_match),'of',nrow(provenance),'\n')
cat('AAPC point matches:',sum(point_checks$agrees_4dp),'of',nrow(point_checks),'\n')
cat('APC point matches:',sum(apc_checks$agrees_4dp),'of',nrow(apc_checks),'\n')
stopifnot(nrow(point_checks)==26,all(point_checks$agrees_4dp),all(apc_checks$agrees_4dp),
 nrow(provenance)==38,all(provenance$exact_rate_match),all(provenance$exact_SE_match))
sessionInfo()
sink()
