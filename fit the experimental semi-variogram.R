

library(ggplot2);library(geoR);library(automap);library(gstat);library(reshape2);library(fields);library(VIM);library(mice)
#via touma
season=c('DJF','MAM','JJA','SON')

all_lon<-seq(0.125,359.875,0.25)
all_lon[which(all_lon>180)]<-all_lon[which(all_lon>180)]-360
all_lat<-seq(59.375,-59.375,-0.25)
all_lat8=all_lat[seq(1,460,8)] #lat列表，取1/8后的
new_land_mask<-read.csv('land_mask.csv')

colnames(new_land_mask)<-all_lon
rownames(new_land_mask)<-all_lat



lat=list.files('D:\\spatial_scales\\data\\7semi_again\\lat') # lat

lat8=as.numeric(substr(lat,5,nchar(lat)-4))#The order of lat
ll=strsplit(l,'_')

s=do.call(rbind,ll)[,2] #find 0/0.5
m=as.numeric( do.call(rbind,ll)[,3])#month
la=as.numeric(substr(do.call(rbind,ll)[,4],1,nchar(do.call(rbind,ll)[,4])-4))

data(meuse)
coordinates(meuse)<-c('x','y')

r=new_land_mask


m1=matrix(c(12,1:11),nrow = 4,byrow = T) ##season


for (k in 1:12) { 
  #range=r[seq(1,460,8),]
  lat1=lat2=lat3=lon1=lon2=lon3={}
  
  for (j in 1:nrow(new_land_mask)) {
    lon=read.csv(paste0('D:\\spatial_scales\\data\\7semi_again\\lat\\',lat[j])) #la这个纬度下的陆地面积
    # 第一个纬度的1月份
    
    s_0.5=read.csv(paste0('D:\\spatial_scales\\data\\7semi_again\\semi_again_1\\',l[which(la==lat8[j] & m==k & s==0.5)])) 
    s_0=read.csv(paste0('D:\\spatial_scales\\data\\7semi_again\\semi_again_1\\',l[which(la==lat8[j] & m==k & s==0)]))


    
    
    sum=which(colSums(s_0)!=0) #找到有数据的列
    sum=sum[-1]
    for (i in sum) {
      name<-c('np','dist','gamma','dir.hor','dir.ver','id')
      
      gamma=0.5*s_0.5[,i]/(s_0.5[,i]+s_0[,i])
      dist=c(1:40)*25-12.5
      np=s_0.5[,i]+s_0[,i]
      vario=matrix(0,nrow = 40,ncol=6)
      vario[,1]=np
      vario[,2]=dist
      vario[,3]=gamma
      vario[,4]=vario[,5]=0
      vario[,6]='var1'
      vario<-as.data.frame(vario)
      colnames(vario)<-name
      class(vario)<-c("gstatVariogram","data.frame")
      vario[,1]<-as.numeric(as.character(vario[,1]))
      vario[,2]<-as.numeric(as.character(vario[,2]))
      vario[,3]<-as.numeric(as.character(vario[,3]))
      
      if(is.na(sum(vario[,3]))==T ){
        vario=vario[-which(is.na(vario[,3]==T)),]
      }
      
      
      vfit=fit.variogram(vario,vgm(model = 'Exp',nugget = 0.5))
      
      # if((3*vfit[2,3])>800){
      #   
      #   png(file=paste0('D:\\spatial_scales\\global_scales\\7range\\',season[mm],'\\more_than_800\\',season[mm],'_',j,'_',lon[1,i],'_',round(3*vfit[2,3],0),'km.png'),width = 240*3, height = 240) 
      #   print(plot(vario,vfit))
      #   dev.off()
      #   lat1=rbind(lat1,j)
      #   lon1=rbind(lon1,lon[1,i])
      # }else if((3*vfit[2,3])<800 & (3*vfit[2,3])>600 ){
      #   png(file=paste0('D:\\spatial_scales\\global_scales\\7range\\',season[mm],'\\600_800\\',season[mm],'_',j,'_',lon[1,i],'_',round(3*vfit[2,3],0),'km.png'),width = 240*3, height = 240) 
      #   print(plot(vario,vfit))
      #   dev.off()
      #   lat2=rbind(lat2,j)
      #   lon2=rbind(lon2,lon[1,i])
      # }else{
      #   png(file=paste0('D:\\spatial_scales\\global_scales\\7range\\',season[mm],'\\less_than_600\\',season[mm],'_',j,'_',lon[1,i],'_',round(3*vfit[2,3],0),'km.png'),width = 240*3, height = 240) 
      #   print(plot(vario,vfit))
      #   dev.off()
      #   
      #   lat3=rbind(lat3,j)
      #   lon3=rbind(lon3,lon[1,i])
      #   
      # }
      
      range[which(all_lat8==lat8[j]),lon[1,i] ]= 3*vfit[2,3]
      
    }
    
    print(j)
    
  }
  # write.csv(cbind(lat1,lon1),paste0('7range\\more_than_800_',season[mm],'.csv'))
  # write.csv(cbind(lat2,lon2),paste0('7range\\600_800_',season[mm],'.csv'))
  # write.csv(cbind(lat3,lon3),paste0('7range\\less_than_600_',season[mm],'.csv'))
  write.csv(range,paste0('7range\\range_',season[mm],'.csv'))
  
}