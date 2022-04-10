
setwd('D:\\spatial_scales\\global_scales')

rm(list=ls())

dist1<-function(x){
  (geoDist(x[1],x[2],x[3],x[4]))/1000
}

library(ncdf4);library(geosphere);library(reshape2);library(gstat);library(SoDA) 

lon<-seq(0.125,359.875,0.25)
lon[which(lon>180)]<-lon[which(lon>180)]-360
lat<-seq(59.375,-59.375,-0.25)


land_mask<-read.csv('land_mask.csv',) #CRU land mask
land_mask<-land_mask[,-1]
s=seq(0,1000,25)
label=seq(1,40,1)

mon_ping<-c(0,31,59,90,120,151,181,212,243,273,304,334,365)
mon_run<-c(0,31,60,91,121,152,182,213,244,274,305,335,366)

x=1
x2=8*x-7

for (x in 1:476){ ##by the lats
  l<-which(land_mask[x,]==1)

  ll<-lon[l] #lon locations in this lat.
  per_l=rbind(t(l),t(ll))
  write.csv(per_l,paste0('lat_',lat[x],'.csv'))
  
  

  land_mask_1<-land_mask
  
  if( lat[x]>=54.5 ){
    land_mask_1[which(lat<= lat[x]-6 ),]<-NA  
    
    lat_range<-c(1:(x+20))
    
  }else if( lat[x]<54.5 & lat[x]> -54.5 ){
    land_mask_1[which( lat>=lat[x]+6  | lat<= lat[x]-6 ),]<-NA  
    
    lat_range<-c( (x-20) :(x+20))
    
  }else{
    land_mask_1[which(lat>= lat[x]+6 ),]<-NA  
    lat_range<-c( (x-20) : 476 )
  }
  
  la<-floor((abs(lat[x])/10)+5) 
  
  
  # by 12 months, 
  for (m in 1:12) {
    mean_ex<-read.csv(paste0('mean_extreme_monthly/mean_extreme_',m,'.csv')) # mean extreme precipition over the 40 years
    mean_ex<-mean_ex[,-1] 
    
    SUM<-{}
    DATE<-{}
    
    sum_0<-sum_0.5<-matrix(0,nrow = 40,ncol = length(l))
    
    for (y in 1983:2018) {
      if(y%%4==0){#闰年
        ma=mon_run[m]+1
        mb=mon_run[m+1]
      }else{
        ma=mon_ping[m]+1
        mb=mon_ping[m+1]
      }
      
      
      ## Convert PERSIANN precipitation to CSV files because of memory reasons, there are 365/366 files for each [name]file
      name<-list.files( paste0('D:\\spatial_scales\\github\\persiann-cdr-csv\\',y )  ) 
      
      for (date in ma:mb) {
        
        

        p=read.csv(paste0('D:\\spatial_scales\\github\\persiann-cdr-csv\\',y,'\\',name[date] ))
        p=p[,-1]
        p<-p*land_mask_1 #land precipitation in that day

        b<- p>=mean_ex #land precipitation in that day 

        b[which(b==T)]=1 
        b=b*p
        
        
        ex=which(b>0, arr.ind=TRUE) #grid point with extreme precipitation
        ex=subset(ex,ex[,1]==x ) #(lat,lon=476*1440) grid point in lat x
        
        
        if( nrow(ex)!=0  ){
          nn=1:nrow(ex)
          for (k in nn) {
            
            
            if( lon[ex[k,2]] <=la & 0 < lon[ex[k,2]]  ){
              d=b[lat_range,  c(1:(ex[k,2]+(la/0.25))  ) ] 
              
              lon1=lon[ c(1:(ex[k,2]+(la/0.25))  )]
              lat1=lat[ lat_range  ]
              
            }else if(-la < lon[ex[k,2]]&lon[ex[k,2]] < 0 ){
              
              d=b[lat_range,  c((ex[k,2]-(la/0.25)):1440 ) ] 
              lon1=lon[ c((ex[k,2]-(la/0.25)):1440 )]
              lat1=lat[lat_range]
            }else{
              
              d=b[lat_range,  c((ex[k,2]-(la/0.25)):(ex[k,2]+(la/0.25)) ) ] 
              lon1=lon[ c((ex[k,2]-(la/0.25)):(ex[k,2]+(la/0.25)) )]
              lat1=lat[lat_range]
            }
            
            
            
            
            e=d
            DIST<-{} #distance between grids 
            for (i in 1:length(lon1)) {
              for (j in 1:length(lat1)) {
                
                dist=(distm( c(lon[ex[k,2]],lat[ex[k,1]]), c(lon1[i],lat1[j]) ))
                DIST=rbind(DIST,dist)
                if(dist>500000){
                  e[j,i]=NA #lat，lon
                }
              }
            }
            #e : all grid points <=500km
            if(sum(e>=0,na.rm = T)>=200){ #if Meet the conditions
              if( sum(e>0,na.rm = T)/sum(e>=0,na.rm = T) >=0.1 ){
                colnames(e)<-lon1
                rownames(e)<-lat1
                
                
                h<-melt(as.matrix(e),na.rm = T)
                #The first column is the lata, the center point is removed and all points in the neighborhood around it are included
                h<-h[-which(h[,1]==lat[x] & h[,2]==lon[ex[k,2]]),] 
                h[which(h[,3]>0),3]<-1
                
                h1lat=rep(as.matrix(h[,1]),each=nrow(h))
                h1lon=rep(as.matrix(h[,2]),each=nrow(h))
                h2lat=rep(as.matrix(h[,1]),times=nrow(h))
                h2lon=rep(as.matrix(h[,2]),times=nrow(h))
                hh=cbind(h1lat,h1lon,h2lat,h2lon)
                
                h31=rep(as.matrix(h[,3]),each=nrow(h))
                h32=rep(as.matrix(h[,3]),times=nrow(h))
                h3=abs(h31-h32)
 
                a=nrow(h)
                aa<-{}
                aa<-rbind(aa,(a*1+1))
                for (n in 2:(a-1)) {
                  aa<-append(aa, c( (a*n+1):(a*n+n)  ))
                }
                n=c(aa , which(h31==0 & h32==0), which(hh[,1]==hh[,3] & hh[,2]==hh[,4]))
                hh<-hh[-n,]
                h3=h3[-n]
                
                g=apply(hh, 1, dist1)
                
                
                g_0<-g[which(h3==0)]
                g_0.5<-g[which(h3==1)]
                
                a<-cut(g_0 ,breaks = s,labels = label)

                s_0<-as.data.frame(table(as.numeric(as.character(a))))

                
                a<-cut(g_0.5,breaks = s,labels = label)

                s_0.5<-as.data.frame(table(as.numeric(as.character(a))))
                
                lon_situ<-which(ll== lon[ex[k,2]]) 
                
                if(nrow(s_0)!=0){
                  s_0[,1]<-as.numeric(as.character(s_0[,1]))
                  s_0[,2]<-as.numeric(as.character(s_0[,2]))
                  sum_0[s_0[,1], lon_situ ]<-sum_0[s_0[,1], lon_situ ]+s_0[,2]
                }
                if(nrow(s_0.5)!=0){
                  s_0.5[,1]<-as.numeric(as.character(s_0.5[,1]))
                  s_0.5[,2]<-as.numeric(as.character(s_0.5[,2]))
                  sum_0.5[s_0.5[,1], lon_situ ]<- sum_0.5[s_0.5[,1], lon_situ ]+s_0.5[,2]
                }
                
                
                
              }
            }

            
          }
          
          
        }
        
      }
      
      
    }
    #Here, all Category 1–1 and Category 0–1 for the experimental semi-variograms within this latitude x over 40 years are obtained
    write.csv(sum_0,paste0('semi-variograms/sum_0_',m,'_',lat[x],'.csv'))
    write.csv(sum_0.5,paste0('semi-variograms/sum_0.5_',m,'_',lat[x],'.csv'))
    
  }
  
  
} 




