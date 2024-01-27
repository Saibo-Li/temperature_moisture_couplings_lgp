var GLDAS = ee.ImageCollection("NASA/GLDAS/V021/NOAH/G025/T3H"),
    geometry_large = /* color: #98ff00 */ee.Geometry.Polygon(
        [[[-180, 85],
          [-180, -60],
          [180, -60],
          [180, 85]]], null, false);
// var GLDAS = ee.ImageCollection("NASA/GLDAS/V20/NOAH/G025/T3H") 

// GLDAS ET kg/m2/s = mm/s
// ESoil_tavg W/m^2 Direct evaporation from bare soil
// ECanop_tavg W/m^2 Canopy water evaporation
// Evap_tavg kg/m^2/s Evapotranspiration
 
// SoilMoi0_10cm_inst kg/m^2 Soil moisture
// SoilMoi10_40cm_inst kg/m^2 Soil moisture
// SoilMoi40_100cm_inst kg/m^2 Soil moisture

function gldas_inputs_daily(date_begin, date_end, dailyImg_iters){
    date_end = ee.Date(date_end).advance(1, 'day'); // for filterDate
    var gldas_raw = ee.ImageCollection(GLDAS).select(['Tair_f_inst','Evap_tavg', 'PotEvap_tavg'])
        .filterDate(date_begin, date_end)
        .map(pkg_date.add_dn(true, 1));
    // print(gldas_raw);
     
    var gldas_day = pkg_agg.aggregate_prop(
        gldas_raw.select(['Tair_f_inst','Evap_tavg', 'PotEvap_tavg']), 
        'dn', 'mean')  
        .map(function(img){
            var Tavg = img.select('Tair_f_inst').subtract(273.15).rename('Tavg');
            var Aet = img.select('Evap_tavg').multiply(86400).rename('Aet'); // ['Evap_tavg'], kg/m2/s = mm/s to mm/day
            var Pet = img.select('PotEvap_tavg').rename('Pet')
            var Pet_mm = Tavg.addBands(Pet).expression('Pet*(1/((2501 - 2.361 * tmean)*1000))', {tmean:Tavg,Pet:Pet}).rename('Pet_mm').multiply(86400) //W/m2 to mm/day
            return Tavg.addBands([Aet, Pet, Pet_mm]).clip(geometry_large).copyProperties(img);
        });
    return gldas_day;
}

// var date_begin = ee.Date.fromYMD(2020,1,1)
// var date_end = date_begin.advance(1,'day')
// var imgcol_gldas = ee.ImageCollection(gldas_inputs_daily(date_begin, date_end, null)).aside(Map.addLayer);

// by month
var year = 2020
for(var j=1;j<2;j++){
  var month = j
  print(month)
  var date_begin = ee.Date.fromYMD(year,month,1)
  var date_end = date_begin.advance(1,'month').advance(-1,'day')
  
  var imgcol_gldas = ee.ImageCollection(gldas_inputs_daily(date_begin, date_end, null));
  
  // define image name for asset !!!new version
  var days = date_end.difference(date_begin, 'day')
  var dateList = ee.List([])
  dateList = ee.List.sequence(0,days).map(function(i){
    i  = date_begin.advance(i, 'day').format('yyyy-MM-dd')
    return i
  }).getInfo()
  
  // judge exit imgs in the folder
  var indexes = ee.ImageCollection("projects/seibertli602/assets/GLDAS_daily").aggregate_array('system:index').getInfo();
  
  function contains(xs, x) {
      for (var i = 0; i < xs.length; i++) {
          if (xs[i] === x) return true;
      }
      return false;
  }
  
  // Batch export 
  var img;
  var n = dateList.length;
  var geometry = geometry_large
  
  for (var i = 0; i < n; i++) {
      var date = dateList[i];
      var task = date;
      
      if (contains(indexes, task)) continue; // if exist then next
  
      img = imgcol_gldas.filterDate(date).first().toFloat();
      Export.image.toAsset({
        image: img, 
        description: task,
        assetId:"projects/seibertli602/assets/GLDAS_daily/"+task,
        region: geometry,
        scale:27830,
        crs: "EPSG:4326",
        maxPixels: 1e13
      }); 
  }
}