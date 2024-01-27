/**** Start of imports. If edited, may not auto-convert in the playground. ****/
var point = /* color: #d63000 */ee.Geometry.Point([116.8859, 35.719148]),
    imageCollection = ee.ImageCollection("users/seibertli602/05Global/ERA"),
    geometry_large = /* color: #98ff00 */ee.Geometry.Polygon(
        [[[-180, 85],
          [-180, -60],
          [180, -60],
          [180, 85]]], null, false);
/***** End of imports. If edited, may not auto-convert in the playground. *****/
var start = 1990,end = 2020;
var years = ee.List.sequence(start,end,1).map(function(i){
  var startDate_mean = ee.Date.fromYMD(i,1,1)
  var endDate_mean = startDate_mean.advance(1,'year')
  var ET_scale = ee.ImageCollection('IDAHO_EPSCOR/TERRACLIMATE')
                  .filter(ee.Filter.date(startDate_mean, endDate_mean))
                  .select(['aet','pet'])
                  .sum()
                  .multiply(0.1)
  return ET_scale.set('system:time_start',startDate_mean.millis())
})
years = ee.ImageCollection(years)
print(years)
/*aet pet*/
var band = 'aet'
years = years.select(band).toBands()

print(years)
Export.image.toDrive({
  image: years.toFloat(), 
  description: "Terra_"+band,
  folder:"ERA5_Land",
  region: geometry_large,
  scale:27830,//11132(1deg),
  crs: "EPSG:4326",
  maxPixels: 1e13
}); 
