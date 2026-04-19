library(ggplot2)
library(lubridate)
library(dplyr)
library(pracma)

# **Step0 API for Step5**
findSOSEOS <- function(df_day_x,peak_thresh,peak_nums){
  
  # con
  df_day_SOS_index = NULL
  df_day_EOS_index = NULL
  df_day_POS_index = NULL
  
  # peak thresh range
  df_day_x$sg_peak_thresh <- df_day_x$sg
  df_day_x$sg_peak_thresh[df_day_x$ratio<peak_thresh]=0
 
  df_day_sg <- data.frame(sgNG = df_day_x$sg[2:length(df_day_x$sg)]-df_day_x$sg[1:(length(df_day_x$sg)-1)])
  df_day_sg$sgNG[df_day_sg$sgNG>=0]=1
  df_day_sg$sgNG[df_day_sg$sgNG<0]=0
  df_day_sg <- rbind(data.frame(sgNG=1),df_day_sg)
  df_day_sg$sgNG[df_day_x$ratio<peak_thresh]=NA
  df_day_x$sgNG <- df_day_sg$sgNG
  
  
  df_day_SOSEOS <- df_day_x[which(df_day_x$sg_peak_thresh!=0),]
  for (peak in 1:peak_nums) {
    print(peak)
    df_day_SOSEOS <- df_day_SOSEOS[order(df_day_SOSEOS$TIMESTAMP),]
    df_day_SOS <- df_day_SOSEOS[which(df_day_SOSEOS$sgNG==1),][1,]
    df_day_SOS$index <- peak
    df_day_SOS_index <- rbind(df_day_SOS_index,df_day_SOS)
    
    df_day_POS <- df_day_SOSEOS[which(df_day_SOSEOS$TIMESTAMP>df_day_SOS$TIMESTAMP & df_day_SOSEOS$sgNG==0),][1,]
    
    df_day_SOSEOS <- df_day_SOSEOS[which(df_day_SOSEOS$TIMESTAMP>=df_day_POS$TIMESTAMP),]
  }
  
 
  df_day_SOSEOS <- df_day_x[which(df_day_x$sg_peak_thresh!=0),]
  df_day_SOSEOS$sgNG <- df_day_SOSEOS$sgNG-1
  for (peak in peak_nums:1) {
    print(peak)
    df_day_SOSEOS <- df_day_SOSEOS[order(df_day_SOSEOS$TIMESTAMP,decreasing = T),]
    df_day_EOS <- df_day_SOSEOS[which(df_day_SOSEOS$sgNG== -1),][1,]
    df_day_EOS$index <- peak
    df_day_EOS_index <- rbind(df_day_EOS_index,df_day_EOS)
    
    df_day_POS <- df_day_SOSEOS[which(df_day_SOSEOS$TIMESTAMP<df_day_EOS$TIMESTAMP & df_day_SOSEOS$sgNG==0),][1,]
    df_day_SOSEOS <- df_day_SOSEOS[which(df_day_SOSEOS$TIMESTAMP<=df_day_POS$TIMESTAMP),]
  }
  
  df_day_SOSEOS_index <- rbind(df_day_SOS_index,df_day_EOS_index)
  df_day_SOSEOS_index <- df_day_SOSEOS_index[order(df_day_SOSEOS_index$TIMESTAMP),]
  
  df_day_SOSPOSEOS_index <- rbind(df_day_SOSEOS_index,df_day_POS_index)
  df_day_SOSPOSEOS_index$sgNG <- as.character(df_day_SOSPOSEOS_index$sgNG)
  
  # plot
  p_x <- ggplot(df_day_x,aes(x=TIMESTAMP)) + 
    scale_color_manual(values = c("grey60", "black", "blue"),
                       breaks=c('GPP','SG',"sg_peak_thresh"),
                       labels=c('GPP','SG',"sg_peak_thresh"))+
    geom_point(aes(y=GPP_DT_VUT_REF,color='GPP'),shape=16, alpha=.7)+
    geom_line(aes(y=sg_peak_thresh,color='sg_peak_thresh'),size=1)+
    geom_line(aes(y=sg,color='SG'))+
    geom_vline(data = df_day_SOSPOSEOS_index,aes(xintercept = TIMESTAMP,linetype=sgNG))
  print(p_x)
  
  df_day_SOSPOSEOS_index[order(df_day_SOSPOSEOS_index$TIMESTAMP),]
} 

site_name <- 'US-Tw3' 
file_name <- list.files(path = r'(.\fluxnet\Crop_site\)',recursive=T,pattern = sprintf('.*%s.*_DD_.*.csv',site_name),full.names=T)
df <- read.csv(file_name)
df$TIMESTAMP <- ymd(as.character(df$TIMESTAMP))
df <- df[,c('TIMESTAMP','GPP_DT_VUT_REF')]
df$year <- year(df$TIMESTAMP)
print(levels(factor(df$year)))

df_day_year <-df %>% 
  group_by(year) %>% 
  summarise(counts = length(year))
df_day_year <- df_day_year[which(df_day_year$counts>0.8*365),]
# filter gte 2003
df_day_year <- df_day_year[which(df_day_year$year>1999),]
if(length(df_day_year$year)==0){print("!!!This site is invalid")}
print(df_day_year)

df <- dplyr::left_join(df,df_day_year,by='year')
df <- df[which(df$counts!=''),]


ggplot(df,aes(x=as.POSIXct(TIMESTAMP),y=GPP_DT_VUT_REF))+
  geom_point()+
  geom_vline(data = df_day_year,aes(xintercept=as.POSIXct(paste(year,'-01-01',sep=''))),
             linetype='dashed',size=1,alpha=.8)+
  scale_x_datetime(date_breaks='1 year')

if(aggregate_flag){
  
  
  dn=8
  df$year_start <- '0101'
  df$year_start <- ymd(paste(as.character(df$year),df$year_start,sep=""))
  df$dn <- as.integer(floor(difftime(df$TIMESTAMP,df$year_start,units = 'days')/dn)+1)
  df$Yeardn <- paste(as.character(df$year),as.character(df$dn),sep="")
  
 
  # aggregate by dn
  df <-df %>% 
    group_by(Yeardn) %>% 
    summarise(TIMESTAMP = min(TIMESTAMP),GPP_DT_VUT_REF=mean(GPP_DT_VUT_REF),dn=first(dn),Yeardn = first(Yeardn),year = first(year))
  df <- df[order(df$TIMESTAMP),]
}
# *************************************
# ***     Start for yearly          ***
# *************************************

debug_year = '2014' #par for year

# filter yearly
{
start_date <- as.Date(paste(debug_year,'-01-01',sep = ''))
end_date <- as.Date(paste(debug_year,'-12-31',sep = ''))
df_day <- df[which(df$TIMESTAMP>=start_date & df$TIMESTAMP<=end_date),]

# sg filter
df_day$sg <- savgol(df_day$GPP_DT_VUT_REF, 13,3)
df_day$sg[df_day$sg<0]=0

# ratio
GPP_max <- max(df_day$sg)
GPP_min <- min(df_day$sg)
df_day$ratio <- (df_day$sg-GPP_min)/(GPP_max-GPP_min)

# *step5-1 identify ture peak numbers*
peak_nums = 1 # par for peaks numers
peak_thresh = 0.8 # par for peaks numers
df_day_SOSPOSEOS_nums_peak <- findSOSEOS(df_day,peak_thresh,peak_nums)

# *step5-2 identify potential sos eos pos numbers*
SOSEOS_thresh <- 0.01 #par for  SOS and EOS
SOSEOS_thresh_peaks = 1 #par for  SOS and EOS
df_day_SOSPOSEOS_nums_soseos <- findSOSEOS(df_day,SOSEOS_thresh,SOSEOS_thresh_peaks)

# step5-3 update df_day_SOSPOSEOS_nums_peak based on the df_day_SOSPOSEOS_nums_soseos*
df_day_SOSPOSEOS_nums_soseos_update <- NULL
if(SOSEOS_thresh_peaks != 1){
  for (i in 1:length(df_day_SOSPOSEOS_nums_peak$Yeardn)) {
    print(i)
    df_day_SOSPOSEOS_nums_peak_x <- df_day_SOSPOSEOS_nums_peak[i,]
    df_day_SOSPOSEOS_nums_soseos_x <- df_day_SOSPOSEOS_nums_soseos
    df_day_SOSPOSEOS_nums_soseos_x$dif <- df_day_SOSPOSEOS_nums_peak_x$TIMESTAMP
    df_day_SOSPOSEOS_nums_soseos_x$dif <- abs(as.numeric(difftime(df_day_SOSPOSEOS_nums_soseos_x$dif,
                                                                  df_day_SOSPOSEOS_nums_soseos_x$TIMESTAMP)))
    df_day_SOSPOSEOS_nums_soseos_x <- df_day_SOSPOSEOS_nums_soseos_x[which(df_day_SOSPOSEOS_nums_soseos_x$sgNG==df_day_SOSPOSEOS_nums_peak_x$sgNG),]
    
    df_day_SOSPOSEOS_nums_soseos_x <- df_day_SOSPOSEOS_nums_soseos_x[order(df_day_SOSPOSEOS_nums_soseos_x$dif),][1,]
    df_day_SOSPOSEOS_nums_soseos_update <- rbind(df_day_SOSPOSEOS_nums_soseos_update,df_day_SOSPOSEOS_nums_soseos_x)
  }
}else{
  df_day_SOSPOSEOS_nums_soseos_update <- df_day_SOSPOSEOS_nums_soseos
}
if(length(df_day_SOSPOSEOS_nums_soseos_update$index) == 3){
  df_day_SOSPOSEOS_nums_soseos_update$index <- 1
}

ggplot(df_day,aes(x=TIMESTAMP)) + 
  scale_color_manual(values = c("grey60","red"),
                     breaks=c('GPP','SG'),
                     labels=c('GPP','SG'))+
  geom_point(aes(y=GPP_DT_VUT_REF,color='GPP'),shape=16, alpha=.7)+
  geom_line(aes(y=sg,color='SG'),size=1)+
  geom_vline(data = df_day_SOSPOSEOS_nums_soseos_update,aes(xintercept = TIMESTAMP,linetype=sgNG))


# p <- ggplot(df_day,aes(x=TIMESTAMP)) + 
#   scale_color_manual(values = c("grey60","red", "black", "blue",'green'),
#                     breaks=c('GPP','ratio','SG',"sg_peak_thresh","sgNG"),
#                     labels=c('GPP','ratio','SG',"sg_peak_thresh","sgNG"))+
#   geom_point(aes(y=GPP_DT_VUT_REF,color='GPP'),shape=16, alpha=.7)+
#   geom_point(aes(y=ratio,color='ratio'),shape=16)+
#   geom_line(aes(y=sg,color='SG'))+
#   geom_line(aes(y=sg_peak_thresh,color='sg_peak_thresh'),size=1)+
#   geom_line(aes(y=sgNG,color='sgNG'),size=1)+
#   geom_vline(xintercept = df_day_SOSPOSEOS_nums$TIMESTAMP)+
#   labs( y = "GPP ("~m^-2~")",x = NULL)+
#   theme_bw()+
#   theme(panel.grid.minor=element_blank(),
#         panel.grid.major = element_line(colour = "white",size=1),
#         panel.background = element_rect(fill = "#E1E1E1"),
#         axis.text = element_text(size = 17, face = "bold",color = 'black'),
#         axis.title = element_text(size = 18, face = "bold",color = 'black'))+
#   theme(
#     legend.text = element_text(size = 16, color = 'black',face = 'bold'),
#     legend.title=element_blank(),
#     legend.position =c(0.5,0.9), #The other is to give numerical position coordinates：
#     legend.justification = c(0.5, 0),
#     legend.background = element_blank(),
#     legend.key.size = unit(20, "pt"),
#     legend.key = element_blank()
#   )+
#   guides(color = guide_legend(nrow = 1))
# # plot(df_day_sg$sgNG, type="l", col="navy")
# p
