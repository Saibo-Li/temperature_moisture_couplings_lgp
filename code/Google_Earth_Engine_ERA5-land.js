/**** Start of imports. If edited, may not auto-convert in the playground. ****/
var ESA_crop = ee.Image("users/seibertli602/05Global/1992_2020_cropland_final"),
    geometry = /* color: #98ff00 */ee.Geometry.Polygon(
        [[[-180, 70],
          [-180, -60],
          [180, -60],
          [180, 70]]], null, false),
    geometry_large = /* color: #98ff00 */ee.Geometry.Polygon(
        [[[-180, 85],
          [-180, -60],
          [180, -60],
          [180, 85]]], null, false);
/***** End of imports. If edited, may not auto-convert in the playground. *****/
// accumulated:
// mean:

var ERA_hourly = ee.ImageCollection("ECMWF/ERA5_LAND/HOURLY")
                .select(['temperature_2m','total_precipitation','surface_runoff','snowmelt',
                'volumetric_soil_water_layer_1','volumetric_soil_water_layer_2','volumetric_soil_water_layer_3',
                'evaporation_from_vegetation_transpiration',
                'evaporation_from_bare_soil','evaporation_from_the_top_of_canopy',
                'potential_evaporation'])
                .filter(ee.Filter.date('1990-01-01', '2021-01-02'));

var accumulated_var = ['total_precipitation','surface_runoff','snowmelt','evaporation_from_vegetation_transpiration',
                'evaporation_from_bare_soil','evaporation_from_the_top_of_canopy',
                'potential_evaporation']
var mean_var = ['temperature_2m','volumetric_soil_water_layer_1','volumetric_soil_water_layer_2','volumetric_soil_water_layer_3']
var start = 1991,end = 2020;
var years = ee.List.sequence(start,end,1).map(function(i){
  //**accumulated**
  var startDate_acc = ee.Date.fromYMD(i,1,2)
  var endDate_acc = startDate_acc.advance(1,'year')
  var era5_acc_yearly = ERA_hourly.filterDate(startDate_acc,endDate_acc)
  .select(accumulated_var)
  .filter(ee.Filter.eq('hour',0))
  .sum()
  
  var Esc = era5_acc_yearly.select(['evaporation_from_bare_soil','evaporation_from_the_top_of_canopy'])
  .reduce(ee.Reducer.sum())
  .rename('Esc')
  .multiply(-1)
  var Ei = era5_acc_yearly.select('evaporation_from_vegetation_transpiration').rename('Ei').multiply(-1)
  var ET0 = era5_acc_yearly.select('potential_evaporation').rename('ET0').multiply(-1)
  
  //**mean**
  var startDate_mean = ee.Date.fromYMD(i,1,1)
  var endDate_mean = startDate_mean.advance(1,'year')
  var days =  endDate_mean.difference(startDate_mean,'day')
  
  var era5_mean_yearly = ERA_hourly.filterDate(startDate_mean,endDate_mean)
  .select(mean_var)
  
  days = ee.ImageCollection(ee.List.sequence(0,days.subtract(1)).map(function(i){
    var current_day = startDate_mean.advance(i,'day')
    return era5_mean_yearly.filterDate(current_day,current_day.advance(1,'day'))
    .mean()
    // .set('system:time_start',current_day.millis())
    // .set('system:index',current_day.format('YYYY_MM_DD'))
  }))

  var mean_days = days.map(function(img){
    var SM = img.expression('layer1*0.07+layer2*0.21+layer3*0.72',{
                layer1:img.select('volumetric_soil_water_layer_1'),
                layer2:img.select('volumetric_soil_water_layer_2'),
                layer3:img.select('volumetric_soil_water_layer_3'),
              }).rename('SM')
              
    var TP = img.select('temperature_2m').subtract(273.15)
    
    var flag_0 = TP.gte(0)
    var flag_5 = TP.gte(5)
    var flag_10 = TP.gte(10)
    var flag_MAB = TP.gt(0).and(TP.lt(30)).rename('flag_MAB')
    
    var TS_0 = flag_0.multiply(TP).rename('TS_0')
    var TS_5 = flag_5.multiply(TP).rename('TS_5')
    var TS_10 = flag_10.multiply(TP).rename('TS_10')
    var MAB = flag_MAB.multiply(TP).rename('MAB')
    
    return TP.addBands([SM,TS_0,TS_5,TS_10,MAB,flag_MAB])
    // .copyProperties(img)
  })
  
  var sum_mean_days = mean_days.select(['TS_0','TS_5','TS_10','MAB','flag_MAB']).sum()
  var MAB_yearly = sum_mean_days.select('MAB').divide(sum_mean_days.select('flag_MAB')).rename('MAB')
  
  return era5_acc_yearly.select(['total_precipitation','surface_runoff','snowmelt'])
    .addBands([Esc,Ei,ET0])
    .addBands(sum_mean_days.select(['TS_0','TS_5','TS_10']))
    .addBands(MAB_yearly)
    .addBands(mean_days.select(['SM','temperature_2m']).mean())
    .clip(geometry_large)
    .set('system:time_start',startDate_mean.millis(),'system:index',startDate_mean.format('YYYY_MM_DD'))
    .toFloat()
})

var ERA_var = ee.ImageCollection(years)
print(ERA_var)

// var vis2mt = {
//   min: 0,
//   max: 1,
//   palette: [
//     '#000080', '#0000D9', '#4000FF', '#8000FF', '#0080FF', '#00FFFF', '#00FF80',
//     '#80FF00', '#DAFF00', '#FFFF00', '#FFF500', '#FFDA00', '#FFB000', '#FFA400',
//     '#FF4F00', '#FF2500', '#FF0A00', '#FF00FF'
//   ]
// };

// Map.addLayer(ERA_var.select('Ei').first().add(ERA_var.select('Esc').first()),vis2mt)

/*test*/
// var dataset1 = ee.ImageCollection("ECMWF/ERA5_LAND/HOURLY")
//                 .filter(ee.Filter.date('1990-01-01', '1991-01-01'))
//                 .select('potential_evaporation')
//                 .filter(ee.Filter.eq('hour',0))
//                 .sum()
//                 // .subtract(273.15)
// Map.addLayer(dataset1,vis2mt)

var ERA_var_list = ERA_var.toList(40)
for(var i=0;i<ERA_var_list.size().getInfo();i++){
  print(i)
  var img = ee.Image(ERA_var_list.get(i))
  var task = ee.Date(img.get('system:time_start')).format('YYYY-MM-DD').getInfo()
  print(task)
  Export.image.toAsset({
    image: img.toFloat(), 
    description: task,
    assetId:"users/seibertli602/05Global/ERA_Other/"+task,
    region: geometry_large,
    scale:27830,//11132(1deg),
    crs: "EPSG:4326",
    maxPixels: 1e13
  }); 
}

