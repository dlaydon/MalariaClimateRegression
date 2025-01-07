
library(raster)
d=read.csv('/home/sjbhatt/pfpr.csv')
d=d[d$continent=="Africa",]
d=d[d$year_start>=2001,]
d=d[d$year_start<=2022,]
d$month_year = paste0(d$month_start,"-",d$year_start)
un = unique(d$month_year)
d$ts = NA
for(i in 1:length(un)){
	month = d$month_start[d$month_year==un[i]][1]
	year = d$year_start[d$month_year==un[i]][1]
	month = sprintf("%02d", month)
	raster_string = paste0('/home/sjbhatt/ts/TSI-Martens2-Pf.',year,month,'.Data.5km.Data.tif')
	r=raster(raster_string)
	NAvalue(r)=-9999
	
	d$ts[d$month_year==un[i]]=extract(r,cbind(d$longitude[d$month_year==un[i]],d$latitude[d$month_year==un[i]]))
}
d$rain = NA
for(i in 1:12){
	month = i
	month = sprintf("%02d", month)
	set_strings=c()
	for(year in 2000:2017){
		raster_string = paste0('/home/sjbhatt/rainfall/Rainfall/5km/Monthly/chirps-v2-0.',year,'.',month,'.sum.5km.NN.tif')
		if(file.exists(raster_string)){
			set_strings = c(set_strings,raster_string)
		}
	}
	st=stack(set_strings)
	r <- calc(st, fun = mean)
	NAvalue(r)=-9999
	d$rain[d$month_start==i]=extract(r,cbind(d$longitude[d$month_start==i],d$latitude[d$month_start==i]))
}

d$rain = asinh(d$rain)


emplogit<-function(Y,N){
	top=Y*N+0.5
	bottom=N*(1-Y)+0.5
	return(log(top/bottom))
}

d_t = d[complete.cases(d[,c('ts','rain','pf_pr','examined'),]),]
d_t$PfPr_logit<-emplogit(d_t$pf_pr,d_t$examined)

m = lm(PfPr_logit~1,data=d_t)
mean(abs(plogis(m$fitted.values)-plogis(d_t$PfPr_logit)))
cor(plogis(m$fitted.values),plogis(d_t$PfPr_logit))


m2 = lm(PfPr_logit~ts+rain,data=d_t)
mean(abs(plogis(m2$fitted.values)-plogis(d_t$PfPr_logit)))
cor(plogis(m2$fitted.values),plogis(d_t$PfPr_logit))

#################
data=d_t

ll.to.xyz<-function(ll){
	if(is.null(colnames(ll))){
		colnames(ll)<-c('longitude','latitude')	
	}
	if(colnames(ll)[1]=='x' & colnames(ll)[2]=='y'){
		colnames(ll)<-c('longitude','latitude')
	}
	if(colnames(ll)[1]=='lon' & colnames(ll)[2]=='lat'){
		colnames(ll)<-c('longitude','latitude')
	}

	ll[,'longitude']<-ll[,'longitude']*(pi/180)
	ll[,'latitude']<-ll[,'latitude']*(pi/180)
	
	x = cos(ll[,'latitude']) * cos(ll[,'longitude'])
	
	y = cos(ll[,'latitude']) * sin(ll[,'longitude'])
	
	z = sin(ll[,'latitude'])

	return(cbind(x,y,z))
}  

xyz<-ll.to.xyz(data[,c('longitude','latitude')])
data<-cbind(data,xyz)


comb<-rbind(data[,c('longitude','latitude')])
xyz<-as.data.frame(ll.to.xyz(comb))
un<-paste(xyz$x,xyz$y,xyz$z,sep=':')
dup<-!duplicated(un)		

mesh = inla.mesh.2d(loc=cbind(xyz[dup,'x'],xyz[dup,'y'],xyz[dup,'z']),
  cutoff=0.01, #should be 0.003
  min.angle=c(25,25),
  max.edge=c(0.015,1) )


 
spde = inla.spde2.matern(mesh,alpha=2)

start_year=2001
end_year=2022
mesh1d=inla.mesh.1d(seq(start_year,end_year,by=4),interval=c(start_year,end_year),degree=2, boundary=c('free'))

est.cov<-as.list(data.frame(ts=d_t$ts,rain=d_t$rain))

A.est =
  inla.spde.make.A(mesh, loc=cbind(data[,'x'],data[,'y'],data[,'z']),group=data[,'year_start'],group.mesh=mesh1d)
  
	#-- Create index matrix --#
field.indices =
  inla.spde.make.index("field", n.spde=mesh$n,n.group=mesh1d$m)

  stack.est =
inla.stack(data=list(response=data$PfPr_logit),
		   A=list(A.est,1),
		   effects=
			 list(c(field.indices),
					c(est.cov)
				  ),
		   tag="est", remove.unused=FALSE,compress=FALSE)

formula<- as.formula(paste(
	paste("response ~ -1 + ts + rain + "),
	paste("f(field, model=spde,group=field.group, control.group=list(model='ar1'))",sep=""),
	sep="")
)

stack.est<-stack.est

#-- Call INLA and get results --#
mod.pred =   inla(formula,
			 data=inla.stack.data(stack.est),
			 family="gaussian",
			 ##################################################
			# control.mode=list(theta=thetac, restart=TRUE,fixed=FALSE), ###
			 ##################################################
			 control.predictor=list(A=inla.stack.A(stack.est), compute=TRUE,quantiles=NULL),
			 control.compute=list(cpo=TRUE, dic=TRUE,config=TRUE),
			 keep=FALSE, verbose=TRUE,#num.threads=5,
			 control.inla= list(strategy = 'gaussian',
			 int.strategy='eb',
			 verbose=TRUE,fast=TRUE,dz=1,
			 step.factor=1,
			 stupid.search=FALSE)
	 )      

index= inla.stack.index(stack.est,"est")$data
lp=mod.pred$summary.linear.predictor$mean[index]

mean(abs(plogis(lp)-plogis(data$PfPr_logit)))
cor(plogis(lp),plogis(data$PfPr_logit))


match.cols<-function(val){
    n=1000
    col<-data.frame(val=seq(min(val),max(val),length.out=n),col=colfunc(n))
    out<-rep(NA,length(col))
    for(i in 1:length(val)){
        out[i]<-as.character(col[which.min(abs(col$val-val[i])),'col'])
    }	
    return(out)
}
bias=1
colfunc <- colorRampPalette(c('yellow','orange','red','brown'),bias=bias)
#colfunc <- colorRampPalette(c("snow1","snow2","snow3","seagreen","orange","firebrick","darkred"), space = "rgb",bias=bias)#colors
#colfunc <- colorRampPalette(c("blue","cyan","yellow","orange","red"), space = "rgb",bias=bias)


s0=plogis(data$PfPr_logit)
s1=plogis(m$fitted.values)
s2=plogis(m2$fitted.values)
s3=plogis(lp)
lat=data$latitude
lon=data$longitude
s0=c(s0,1,0)
s1=c(s1,1,0)
s2=c(s2,1,0)
s3=c(s3,1,0)
lat=c(lat,10,10)
lon=c(lon,10,10)


c0 = match.cols(s0)
c1 = match.cols(s1)
c2 = match.cols(s2)
c3 = match.cols(s3)
library(scales)
cex=0.5
pch=16
par(mfrow=c(2,2),mar=c(2,2,2,2))
plot(lon,lat,col=alpha(c0,0.5),cex=cex,pch=pch,xlab="",ylab="")
plot(lon,lat,col=alpha(c1,0.5),cex=cex,pch=pch,xlab="",ylab="")
plot(lon,lat,col=alpha(c2,0.5),cex=cex,pch=pch,xlab="",ylab="")
plot(lon,lat,col=alpha(c3,0.5),cex=cex,pch=pch,xlab="",ylab="")

library("maps")
library("raster")

# Set up the layout with reduced margins
par(mfrow = c(2, 2), mai = c(1, 1, 1, 1), mar = c(2, 2, 2, 2))

n <- 5 # number in legend

# Loop through your data and create plots
for (i in 0:3) {
  dat <- data.frame(lon, lat, pr = get(paste0("s", i)))
 # map('world', xlim = c(-25, 55), ylim = c(-35, 40), lwd = 0.5, col = "grey95", fill = TRUE, interior = TRUE)
 # title("")
 # map.axes()
  plot(dat$lon, dat$lat, cex = 0.5, pch = 18, col = get(paste0("c", i)))
  legend("bottomleft", title = "PfPr", legend = round(seq(min(get(paste0("s", i))), max(get(paste0("s", i))), length.out = n), 1), col = colfunc(n), pch = 20)
}

out = cbind(lon,lat,s0,s1,s2,s3)
write.csv(out,'/home/sjbhatt/climate.csv')


 d=read.csv('~/Downloads/climate.csv')
 d=d[,2:ncol(d)]
 d2=d[1:(nrow(d)-2),]
 
 p1 = ggplot(data.frame(x=d[,1],y=d[,2],color_variable=d[,3]), aes(x, y, color = color_variable)) +
  geom_point(size = 0.5,alpha=0.5) +
  scale_color_viridis_c(name = expression(italic("PfPr"))) +
  theme_minimal() +
    theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank()
  )

 p2 = ggplot(data.frame(x=d[,1],y=d[,2],color_variable=d[,5]), aes(x, y, color = color_variable)) +
  geom_point(size = 0.5,alpha=0.5) +
  scale_color_viridis_c(name = expression(italic("PfPr"))) +
  theme_minimal()+
    theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank()
  )

  
   p3 = ggplot(data.frame(x=d2[,1],y=d2[,2],color_variable=d2[,5]), aes(x, y, color = color_variable)) +
  geom_point(size = 0.5,alpha=0.5) +
  scale_color_viridis_c(name = expression(italic("PfPr"))) +
  theme_minimal()+
    theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank()
  )

  
   p4 = ggplot(data.frame(x=d[,1],y=d[,2],color_variable=d[,6]), aes(x, y, color = color_variable)) +
  geom_point(size = 0.5,alpha=0.5) +
  scale_color_viridis_c(name = expression(italic("PfPr"))) +
  theme_minimal()+
    theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank()
  )

  
  grid.arrange(p1,p2,p3,p4,ncol=2)