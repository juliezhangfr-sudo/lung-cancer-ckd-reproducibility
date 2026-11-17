# Round-two regenerated observed-rate figures. Base R only; no counts changed.
# Input annual rates from recovered C34_N18_REPRODUCTION; original export provenance retained in round-one manifest.
args<-commandArgs(trailingOnly=TRUE);root<-if(length(args)) args[1] else if(dir.exists('analysis/input')) 'analysis' else 'round2_revision/analysis'
annual<-read.csv(file.path(root,'input/annual_data.csv'),stringsAsFactors=FALSE)
st<-read.csv(file.path(root,'input/state_pooled_2018_2024.csv'),stringsAsFactors=FALSE)
AGES<-c('45-54 years','55-64 years','65-74 years','75-84 years','85+ years')
# Explicitly preserve the unreliable 2015 age-specific gap; no interpolation.
annual$Rate[annual$Disease=='CKD' & annual$Dimension=='Age' & annual$Group=='45-54 years' & annual$Year==2015]<-NA
save_figure<-function(name,w,h,draw) {
 png(file.path(root,'results',paste0(name,'.png')),width=w,height=h,units='in',res=300,type='cairo');draw();dev.off()
 svg(file.path(root,'results',paste0(name,'.svg')),width=w,height=h);draw();dev.off()
 pdf(file.path(root,'results',paste0(name,'.pdf')),width=w,height=h);draw();dev.off()
}
# Sourced by 02_analysis.R. Scientific plots exported as 300-dpi PNG and editable SVG.
palette<-c('#226494','#C76D36','#459887','#9465A5','#C59B37','#687787')
draw_trends<-function(d,dim,title,ckd=FALSE) {
  q<-d[d$Dimension==dim | (dim=='Sex' & d$Dimension=='Overall'),];if(dim=='Age') q<-q[q$Group%in%AGES,]
  groups<-if(dim=='Age')AGES else sort(unique(q$Group))
  ymax<-max(q$Rate,na.rm=TRUE)*1.18
  plot(NA,xlim=range(q$Year),ylim=c(0,ymax),xlab='Year',ylab=if(dim=='Age')'CMR per 100,000' else 'AAMR per 100,000',main=title,bty='l',las=1,cex.main=1.02)
  if(ckd && dim!='Race')rect(2010.5,par('usr')[3],2012.5,par('usr')[4],col='#E7E7E7',border=NA)
  if(dim=='Race')abline(v=2020.5,col='#999999',lty=3)
  for(i in seq_along(groups)) {
    z<-q[q$Group==groups[i],];z<-z[order(z$Year),]
    # Expand the calendar grid to preserve omitted, suppressed, and unreliable gaps.
    grid<-merge(data.frame(Year=seq(min(q$Year),max(q$Year))),z,by='Year',all.x=TRUE)
    if(dim=='Race') {
      for(era in c('Bridged','Single')) {
        v<-grid;v$Rate[is.na(v$Era)|v$Era!=era]<-NA
        lines(v$Year,v$Rate,col=palette[i],lwd=1.6)
      }
    } else lines(grid$Year,grid$Rate,col=palette[i],lwd=1.6)
    if(ckd && dim=='Age') points(z$Year,z$Rate,pch=16,cex=.46,col=palette[i])
  }
  if(ckd && dim=='Age') { abline(v=2015,lty=3,col='#999999'); mtext('45-54 years: 2015 omitted, unreliable (17 deaths)',side=3,line=.05,cex=.70,col='#555555') }
  legend('topright',legend=gsub('85+ years','\u226585 years',groups,fixed=TRUE),col=palette[seq_along(groups)],lwd=2,bty='n',cex=.80,ncol=if(length(groups)>3)2 else 1)
}
main_panel<-function(disease) {
  d<-annual[annual$Disease==disease,];ckd<-disease=='CKD'
  par(mfrow=c(3,2),mar=c(3.3,4.1,2.3,1),mgp=c(2.5,.6,0),oma=c(1,0,2.5,0),family='sans')
  dims<-c('Sex','Age','Race','Region','Urbanicity')
  titles<-c('A. Overall and sex','B. Age','C. Race/ethnicity','D. Census region','E. Historical urbanicity')
  for(i in seq_along(dims))draw_trends(d,dims[i],titles[i],ckd)
  metric<-if(ckd)'CKD_AAMR' else 'All_lung_AAMR'
  s<-head(st[order(-st[[metric]]),],8);s<-s[rev(seq_len(nrow(s))),]
  par(mar=c(3.3,6.3,2.3,1))
  barplot(s[[metric]],names.arg=s$Group,horiz=TRUE,las=1,col=palette[if(ckd)2 else 1],border=NA,
          xlab='Pooled 2018-2024 AAMR per 100,000',main='F. Highest observed state point estimates',cex.names=.90,cex.main=.92)
  mtext(if(ckd)'Lung cancer mortality with co-recorded CKD' else 'Overall lung cancer mortality',outer=TRUE,font=2,cex=1.25,line=.5)
  mtext('Observed rates. Race series broken at 2021; historical urbanicity ends in 2020.',outer=TRUE,side=1,cex=.7,line=0)
}
save_figure('Figure1_all_lung',9,9,function()main_panel('all'))
save_figure('Figure2_lung_CKD',9,9,function()main_panel('CKD'))

