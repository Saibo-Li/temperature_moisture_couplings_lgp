# #######################################################################################
# Determine the required verification points according to the base class and time range #
# Acquire EOS，SOS，LSP，crop identity for each site                                    #
#########################################################################################

library(ggplot2)
library(lubridate)
library(phenocamr)
library("phenor")
library("maps")
library("raster")
library(lubridate)
library(dplyr)

# parameters
################################### 
# site_index;1....                #
# site_name;'aafcottawacfiaf14e'  #
# roi_id;e.g.,1000||2000          #
# period_time;1||3                #
# gcc_x;gcc_50;                   #
# transition_x;transition_50      #
################################### 

sites <- list_sites()
# filter crop sites
sites_AG <- sites[which(sites$primary_veg_type=='AG'),]
sites_AG$date_start = as.Date(sites_AG$date_start)
sites_AG$date_end = as.Date(sites_AG$date_end)
# filter all date range > one year
sites_AG$days = difftime(sites_AG$date_end,sites_AG$date_start,units = 'days')
sites_AG <- sites_AG[which(sites_AG$days > 365),]
# get rois
rois <- list_rois()

site_index <- 1 #par 1
sites_AG_0 <- sites_AG[site_index,]
site_name <- sites_AG_0$site #'shangqiu'
rois_0 <- rois[which(rois$site==site_name),]

if(length(rois_0$roi_id_number) == 1){
  roi_id <- rois_0$roi_id_number #par 2
}else{
  print("Need to manually set roi_id!!!")
  print(rois_0$roi_id_number)
}

# download the three day time series for deciduous broadleaf data at the 
# Bartlett site and will estimate the phenophases (spring + autumn). 
period_time = 3 #par 3
phenocamr::download_phenocam(
  frequency = period_time,
  veg_type = "AG",
  roi_id = roi_id,
  site = site_name,
  phenophase = TRUE,
  out_dir = "."
)

# filter valid date range > 90% of one year from **1day.csv or **3day.csv 
gcc_x <- 'gcc_50' # par 4
day_csv = sprintf("./phenocam/%s_AG_%s_%sday.csv",site_name,roi_id,period_time)
df_day <- read.table(day_csv, header = TRUE, sep = ",")
df_day <- df_day[which(df_day[,gcc_x]!=''),]
# aggregate by year
df_day_year <-df_day %>% 
  group_by(year) %>% 
  summarise(counts = length(doy))
df_day_year <- df_day_year[which(df_day_year$counts<0.8*365/period_time),]
# filter gte 2003
df_day_year <- df_day_year[which(df_day_year$year>2002),]
if(length(df_day_year$year)==0){print("!!!This site is invalid")}
print(df_day_year)

df_day <- dplyr::left_join(df_day,df_day_year,by='year')
df_day <- df_day[which(df_day$counts!=''),]

# read in the transition date file
transition_dates_csv <- sprintf("./phenocam/%s_AG_%s_%sday_transition_dates.csv",site_name,roi_id,period_time)
# "archboldavir_AG_1000_3day_transition_dates.csv"
td <- read.table(transition_dates_csv,
                 header = TRUE,
                 sep = ",")

# select the rising (spring dates) and falling for 50% threshold of Gcc 50
td <- td[td$gcc_value == gcc_x,]
transition_x <- 'transition_10' #par 5
td$year = year(td[,transition_x])
td <- dplyr::left_join(td,df_day_year,by='year')
td <- td[which(td$counts!=''),]

# gcc_90
# falling
# transition_100
if (T){
  # read in the transition date file
  phenophases_x <- function (thresh,gcc_thresh_x){
    tdPOS <- phenophases(day_csv,
                         internal = TRUE,
                         upper_thresh = thresh
    )
    tdPOS_falling <- tdPOS$falling[tdPOS$falling$gcc_value == gcc_thresh_x,]
    if(length(tdPOS_falling$gcc_value) != 0){
      transition_thresh_x <- paste('transition_',thresh*100,sep="")
      tdPOS_falling$POS <- as.Date(tdPOS_falling[,transition_thresh_x], origin = "1970-01-01")
      tdPOS_falling$transition <- thresh
      tdPOS_falling$year <- year(tdPOS_falling$POS)
      tdPOS_falling <- tdPOS_falling[c('year','POS','transition')]
    }else{
      tdPOS_falling <- data.frame(year= numeric(), POS= numeric(), transition=numeric())
    }
    tdPOS_falling
  }
  tdPOS_falling100 <- phenophases_x(1,"gcc_90")
  tdPOS_falling95 <- phenophases_x(0.95,"gcc_90")
  tdPOS_falling90 <- phenophases_x(0.9,"gcc_90")
  
  # merge
  tdPOS_falling_merge <- rbind(rbind(tdPOS_falling100,tdPOS_falling95),tdPOS_falling90)
  tdPOS_falling_merge_min <-tdPOS_falling_merge %>% 
    group_by(year) %>% 
    summarise(POS_min = min(POS))
  tdPOS_falling_merge_min <- dplyr::left_join(tdPOS_falling_merge_min,df_day_year,by='year')
  tdPOS_falling_merge_min <- tdPOS_falling_merge_min[which(tdPOS_falling_merge_min$counts!=''),]
}

roistats_csv = sprintf("%s_AG_%s_roistats.csv",site_name,roi_id)
csvurl = sprintf("%s/%s/ROI/%s",baseurl,site_name,roistats_csv)
df = read.csv(url(csvurl),comment.char="#",header=TRUE)
df$date = as.Date(df$date)
df$year = year(df$date)

df$datetime = with(df, as.POSIXct(paste(date,local_std_time), format="%Y-%m-%d %H:%M:%S"))
# filter daytime
df = df[which(df$solar_elev>0),]
df <- dplyr::left_join(df,df_day_year,by='year')
df <- df[which(df$counts!=''),]

# **Step5 Plot**
df_day_yearly <- df_day
df_yearly <- df

P <- ggplot(df_yearly,aes(x=datetime,y=gcc)) + 
  geom_point(size=0.01, shape=16, alpha=.7,color='grey60')+
  geom_point(data=df_day_yearly,aes(x=as.POSIXct(date),y=gcc_50),
             alpha=.7,size=1.5, shape=21, color = "black",fill='yellow')+
  geom_line(data=df_day_yearly,aes(x=as.POSIXct(date),y=smooth_gcc_50),
  alpha=.7, color='black')+
  geom_vline(data = td,aes(xintercept=as.POSIXct(transition_10),color=direction),
                                     linetype='dashed',size=1,alpha=.8)+
  scale_color_manual(values=c("#005a32","#fc8d59"))+
  geom_vline(data = tdPOS_falling_merge_min,aes(xintercept=as.POSIXct(POS_min)),color='black',
                                     linetype='dashed',size=1,alpha=.8)

exportAll_check_valid = NULL
export_years <- c(2021)
td <- td[order(td$transition_10),]
tdPOS_falling_merge_min <- tdPOS_falling_merge_min[order(tdPOS_falling_merge_min$POS_min),]
species <- c("Pasture")#,"Miscanthus","Glycine max","corn"
species_info <- data.frame(year=export_years,
                               species=species)

exportAll_check_valid <- export_data(export_years,species_info,
                                     tdPOS_falling_merge_min,
                                     td,exportAll_check_valid,F,F,F)
exportAll <- rbind(exportAll,exportAll_check_valid)

# **Step 6-2 manual process for invalid years and export data**
# invalid years mean tha this year miss one of SOS, EOS, and POS
Debug = T
Debug_years <- c(2016) # set invalid years
# set invalid period for this invalid year
Debug_years_info <- data.frame(year=Debug_years,
                               fill_invalid=c('rising'), #falling,rising,both,test
                               fill_invalid_POS=c(F))
species <- c('wheat')
species_info <- data.frame(year=Debug_years,
                           species=species)

if(Debug){
  df_day_yearly_echart <- df_day_yearly[which(df_day_yearly$year %in% Debug_years),]
  tdPOS_falling_merge_min_echart <- tdPOS_falling_merge_min
  td_echart <- td
  
  for(Debug_year in Debug_years){
    print(Debug_year)
    df_Debug_year <- df_day_yearly_echart[which(df_day_yearly_echart$year == Debug_year),]
    td_debug <- td_echart[1,]
    Debug_year_info <- Debug_years_info[which(Debug_years_info$year == Debug_year),]
    
    smooth_max <- max(df_Debug_year$smooth_gcc_50)
    smooth_min <- min(df_Debug_year$smooth_gcc_50)
    smooth_10 <- (smooth_max-smooth_min)*0.1+smooth_min
    
    smooth_max_doy <- df_Debug_year[which(df_Debug_year$smooth_gcc_50 == smooth_max),]$doy
    
    # fill SOS or EOS 
    # Ratio 10%
    fill_invalid =Debug_year_info$fill_invalid
    if(fill_invalid == 'falling'){
      df_Debug_year_EOS <- df_Debug_year[which(df_Debug_year$doy>smooth_max_doy & df_Debug_year$smooth_gcc_50<=smooth_10),]
      df_Debug_year_EOS <- df_Debug_year_EOS[order(df_Debug_year_EOS$doy),][1,]$date
      
      # update td
      td_debug$direction <- 'falling'
      td_debug$transition_10 <- df_Debug_year_EOS
      td_debug$year <- Debug_year
     
      
    }else if(fill_invalid == 'rising'){
      df_Debug_year_SOS <- df_Debug_year[which(df_Debug_year$doy<smooth_max_doy & df_Debug_year$smooth_gcc_50<=smooth_10),]
      if(length(df_Debug_year_SOS$date) == 0){
        df_Debug_year_SOS <- paste(Debug_year,'-01-01',sep = '')
      }else{
        df_Debug_year_SOS <- df_Debug_year_SOS[order(df_Debug_year_SOS$doy,decreasing = T),][1,]$date
      }
      
      # update td
      td_debug$direction <- 'rising'
      td_debug$transition_10 <- df_Debug_year_SOS
      td_debug$year <- Debug_year
      
    }else if(fill_invalid == 'both'){
      df_Debug_year_SOS <- df_Debug_year[which(df_Debug_year$doy<smooth_max_doy & df_Debug_year$smooth_gcc_50<=smooth_10),]
      if(length(df_Debug_year_SOS$date) == 0){
        df_Debug_year_SOS <- paste(Debug_year,'-01-01',sep = '')
      }else{
        df_Debug_year_SOS <- df_Debug_year_SOS[order(df_Debug_year_SOS$doy,decreasing = T),][1,]$date
      }
      # update td SOS
      td_debug_SOS <- td_debug
      td_debug_SOS$direction <- 'rising'
      td_debug_SOS$transition_10 <- df_Debug_year_SOS
      td_debug_SOS$year <- Debug_year
      
      df_Debug_year_EOS <- df_Debug_year[which(df_Debug_year$doy>smooth_max_doy & df_Debug_year$smooth_gcc_50<=smooth_10),]
      df_Debug_year_EOS <- df_Debug_year_EOS[order(df_Debug_year_EOS$doy),][1,]$date
      
      # update td EOS
      td_debug_EOS <- td_debug
      td_debug_EOS$direction <- 'falling'
      td_debug_EOS$transition_10 <- df_Debug_year_EOS
      td_debug_EOS$year <- Debug_year
      
      td_debug <- rbind(td_debug_SOS,td_debug_EOS)
      
    }else{
      td_debug <- NULL
    }
    
    td_echart <- rbind(td_echart,td_debug)
    
    # fill POS
    fill_invalid_POS = Debug_year_info$fill_invalid_POS
    if(fill_invalid_POS == T){
      tdPOS_falling_merge_min_echart <- rbind(tdPOS_falling_merge_min_echart,data.frame(
                                         year=Debug_year,
                                         POS_min=df_Debug_year[which(df_Debug_year$smooth_gcc_50 == smooth_max),]$date,
                                         counts=122
                                       ))
    }
  }
  
  # *plot R echart*
  p_echart <- echartr(df_day_yearly_echart, x=date, y=smooth_gcc_50, type='curve') %>% setToolbox(pos=3)
  p_echart <- p_echart %>% setXAxis(name='y/m/d', axisLabel=list(
    formatter='%y%m/%d', rotate=90)) 
  
  # POS
  data <- tdPOS_falling_merge_min_echart
  data$series <- 'POS'
  data$name1 <- 'POS'
  data$type <- NA
  data$POS_min <- as.character(data$POS_min)
  data$echart_POS_min <- data$POS_min
  # search nearest date for POS
  for(i in data$POS_min){
    df_day_yearly_echart$echart_ML_POS <- abs(difftime(df_day_yearly_echart$date,i,units = 'days'))
    data[which(data$POS_min == i),]$echart_POS_min <- df_day_yearly_echart[order(df_day_yearly_echart$echart_ML_POS),][1,]$date
  }
  
  data$xAxis1 <- data$echart_POS_min
  data$yAxis1 <- 0.5
  data$xAxis2 <- data$echart_POS_min
  data$yAxis2 <- 0.0
  
  # SOS EOS
  dataSOSEOS <- td_echart
  dataSOSEOS$series <- dataSOSEOS$direction
  dataSOSEOS$name1 <- dataSOSEOS$direction
  dataSOSEOS$type <- NA
  dataSOSEOS$echart_transition_10 <- dataSOSEOS$transition_10
  # search nearest date for SOS and EOS
  for(i in dataSOSEOS$transition_10[1:2]){
    df_day_yearly_echart$echart_ML <- abs(difftime(df_day_yearly_echart$date,i,units = 'days'))
    dataSOSEOS[which(dataSOSEOS$transition_10 == i),]$echart_transition_10 <- df_day_yearly_echart[order(df_day_yearly_echart$echart_ML),][1,]$date
  }
  
  dataSOSEOS$xAxis1 <- dataSOSEOS$echart_transition_10
  dataSOSEOS$yAxis1 <- 0.5
  dataSOSEOS$xAxis2 <- dataSOSEOS$echart_transition_10
  dataSOSEOS$yAxis2 <- 0.0
  
  p_echart <- p_echart %>% addML(series='POS', data=data) %>% 
    addML(series=c('rising','falling'), data=dataSOSEOS) %>% 
    setXAxis(splitLine=list(show=FALSE)) %>%
    setYAxis(splitLine=list(show=FALSE))
  p_echart
  
  # *export the left data*
  exportAll_check_invalid = NULL
  reverse=F # for SOS and EOS
  
  # for POS
  alter_POS = F
  pickup_after = F
  
  exportAll_check_invalid = export_data(Debug_years,species_info,
                                        tdPOS_falling_merge_min_echart,
                                        td_echart,exportAll_check_invalid,reverse,alter_POS,pickup_after)
  
}

manually = F

# **optional Step 8 K curve**
# if transition_10 is invalid, you can use this code
df_day_K <- df_day[,c("date","doy","smooth_gcc_50")]
colnames(df_day_K) <- c("date","doy","smooth")
df_day_K_iPlus <- rbind(df_day_K[2:length(df_day_K$doy),],data.frame(date=df_day_K$date[length(df_day_K$date)],
                                                                     doy=df_day_K$doy[length(df_day_K$doy)],smooth=0))
# (i+1)-i
df_day_K$dir <- (df_day_K_iPlus$smooth-df_day_K$smooth)
df_day_K <- df_day_K[1:length(df_day_K$doy)-1,]
df_day_K$dirFlag <-  as.numeric(as.character(cut(df_day_K$dir,breaks=c(-Inf,0,Inf),labels=c(-1,1))))

# movemean for decrease or increase flag
mav <- function(x,n){stats::filter(x,rep(1/n,n), sides=1)}
df_day_K$dirFlag_movemean <- mav(df_day_K$dirFlag,2)

# search knee
doy_flags <- df_day_K[which(df_day_K$dirFlag_movemean==0),]$doy
knee_all <- NULL
knee_flag <- 0
df_day_K_iter <- df_day_K[which(df_day_K$dirFlag_movemean !=''),]
for(doy_flag in doy_flags){
  print(doy_flag)
  df_day_K_flag <- df_day_K_iter[which(df_day_K_iter$doy<doy_flag),]
  if(length(df_day_K_flag$doy) < 10) next
  
  smooth_max <- max(df_day_K_flag$smooth)
  smooth_min <- min(df_day_K_flag$smooth)
  df_day_K_flag$ratio <- (df_day_K_flag$smooth-smooth_min)/(smooth_max-smooth_min)
  
  # next iter
  knee_flag <- knee_flag+1
  df_day_K_iter <- df_day_K_iter[which(df_day_K_iter$doy>doy_flag),]
  
  if(mean(df_day_K_flag$dirFlag) == 1){
    knee <- data.frame(t(kneedle(df_day_K_flag$doy,df_day_K_flag$smooth)))
    if(length(knee) == 0) next
    colnames(knee) <- c('doy','smooth')
    
    knee_10 <- df_day_K_flag[which(df_day_K_flag$ratio>=0.1),]
    knee_10 <- knee_10[order(knee_10$ratio),][1,c('doy','smooth')]
    
    knee <- rbind(knee,knee_10)
    knee <- knee[order(knee$doy),][1,]
  }else{
    knee <- data.frame(t(kneedle(df_day_K_flag$doy,df_day_K_flag$smooth, concave = TRUE,decreasing = TRUE)))
    if(length(knee) == 0) next
    colnames(knee) <- c('doy','smooth')
    
    knee_10 <- df_day_K_flag[which(df_day_K_flag$ratio>=0.1),]
    knee_10 <- knee_10[order(knee_10$ratio),][1,c('doy','smooth')]
    
    knee <- rbind(knee,knee_10)
    knee <- knee[order(knee$doy,decreasing = T),][1,]
  }
  
  knee$knee_flag <- knee_flag
  knee_all <- rbind(knee_all,knee)
  
}

P+
  geom_point(data = df_day_K, aes(x=as.POSIXct(date),y=smooth,fill=as.character(df_day_K$dirFlag)),shape=21)+
  geom_point(data = df_day_K, aes(x=as.POSIXct(date),y=sg),shape=1)+
  scale_fill_manual(values = c('black','blue'))+
  theme(legend.position = 'None')
  geom_point(data = df_day_K, aes(x=as.POSIXct(date),y=smooth),shape=1)+
  geom_point(data = knee_all, aes(x=as.POSIXct(date),y=smooth),color='blue',shape=3)