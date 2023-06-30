// ****************************************************************************************************************** //
// ****************************** Utilities for the ee-synthetic-earth package ************************************** //
// ****************************************************************************************************************** //

/**
 * Performs the NBR2 threshold algorithm approach to detect bare pixels in a collection of images.
   Method described in Dvorakova et al., 2023
 * @param {ImageCollection} col: The input rainfall image collection
 * @param {Date} start_date: A timestamp corresponding to the start of the time interval
 * @param {Date} end_date: A timestamp corresponding to the end of the time interval
 * @param {Geometry} geom: The geometry of the area of interest to retrieve the dry periods
 * @param {Float} diff_tol: The number of days tolerance for counting dry spells
                            (defaults to 0, which means any dry spell longer than a day is extracted)
 * @returns {Dictionary} dry_period_dates: Dictionary of dry period date tuples ([start_date, end_date])
                                           in order of length (longest to shortest)
 * @ignore
*/
exports.getDryPeriod = function(col, start_date, end_date, geom, options){

  var diff_tol = options.diff_tol || 0;

  // Filter data
  var col_filt = col.filterDate(start_date, end_date)
                    .sort('system:time_start');

  var col_arr = col_filt.toArray();
  col_arr = col_arr.arrayReshape(col_arr.arrayLength(0).toArray(), 1);

  // Compute the forward difference
  var difference = _forwardDifference(col_arr);
  var ones = ee.Image(ee.Array([1]))
  difference = ones.addBands(difference).toArray(0); // Prepend a 1 to the differences.

  // Make an index array sized for the length of the data in each pixel.
  var maxSize = col_filt.size();
  var indexes = ee.Image.constant(ee.Array(ee.List.sequence(0, maxSize)));
  indexes = indexes.arraySlice(0, 0, col_arr.arrayLength(0));

  var runStarts = indexes.arrayMask(difference.abs().gte(diff_tol));
  var runValues = col_arr.arrayMask(difference.abs().gte(diff_tol));

  // Append an extra index to mark the end of the last run.
  var runLengths = runStarts.addBands(col_arr.arrayLengths()).toArray(0);
  runLengths = _forwardDifference(runLengths.multiply(-1));

  // Sort lengths of dry periods, and sort the dry spree start in the same order
  var dryRunLengths = runLengths.arrayMask(runValues.eq(0)).arraySort();
  var dryRunStarts = runStarts.arrayMask(runValues.eq(0)).arraySort(dryRunLengths);

  // Extract indices for the first three longest dry runs
  var maxIndex = dryRunLengths.arrayArgmax().arrayGet(0);
  var maxIndex2 = dryRunLengths.arraySlice(0, 0, -3).arrayArgmax().arrayGet(0);
  var maxIndex3 = dryRunLengths.arraySlice(0, 0, -6).arrayArgmax().arrayGet(0);

  // Get start and end of each dry period
  var dry_period_start = dryRunStarts.arrayGet(maxIndex);
  var dry_period_end = dryRunStarts.arrayGet(maxIndex).add(dryRunLengths.arrayGet(maxIndex)).add(10).clamp(0, 365);
  var dry_period_start2 = dryRunStarts.arrayGet(maxIndex2);
  var dry_period_end2 = dryRunStarts.arrayGet(maxIndex2).add(dryRunLengths.arrayGet(maxIndex2)).add(10).clamp(0, 365);
  var dry_period_start3 = dryRunStarts.arrayGet(maxIndex3);
  var dry_period_end3 = dryRunStarts.arrayGet(maxIndex3).add(dryRunLengths.arrayGet(maxIndex3)).add(10).clamp(0, 365);

  // Extract the dry period dates for the first three longest dry periods for the entire area as a dictionary
  var dry_period_dates = ee.Dictionary(dry_period_start.rename('start')
                         .addBands(dry_period_end.rename('end'))
                         .addBands(dry_period_start2.rename('start2'))
                         .addBands(dry_period_end2.rename('end2'))
                         .addBands(dry_period_start3.rename('start3'))
                         .addBands(dry_period_end3.rename('end3'))
                         .reduceRegion({reducer: ee.Reducer.mode(),
                                        geometry: geom,
                                        scale: 100,
                                        crs: 'EPSG:4326',
                                        tileScale: 4,
                                        maxPixels: 1e13}));
  return dry_period_dates;
}

// computes the foward difference of an array image.
function _forwardDifference(image) {
  var left = image.arraySlice(0, 0, -1)
  var right = image.arraySlice(0, 1)
  return left.subtract(right)
}


/**
 * A small utility to add a 'doy' attribute to images in a collection
 * @param {Image} img: The input Image
 * @returns {Image}: An image with a set 'doy' attribute corresponding to the day of year
 * @ignore
*/
exports.addDate = function(img){
      var doy = img.date().getRelative('day', 'year');
      return img.set('doy', doy);
    };

/**
 * Calculate Mann Kendall's S statistic.
 This function returns the Mann Kendall's S statistic, assuming that n is
 less than 40. The significance of a calculated S statistic is found in
 table A.30 of Nonparametric Statistical Methods, second edition by
 Hollander & Wolfe.
 * @param {imageCollection} collection: Input image collection for which to calculate the mann-kendall S statistic.
 * @returns {Image} mk_stat: The Mann Kendall S statistic.
*/
exports.mannKendall = function(collection){
  var afterFilter = ee.Filter.lessThan({
    leftField: 'system:time_start',
    rightField: 'system:time_start'
  });

  var joined = ee.ImageCollection(ee.Join.saveAll('after').apply({
    primary: collection,
    secondary: collection,
    condition: afterFilter
  }));

  var sign = function(i, j) { // i and j are images
    //return ee.Image(j).neq(i) // Zero case
    //    .multiply(ee.Image(j).subtract(i).clamp(-1, 1)).int();
      var concordant = ee.Image(i).lt(j).rename('concordant');
      var discordant = ee.Image(i).gt(j).rename('discordant');
      return concordant.addBands(discordant);
  };

  var mk = ee.ImageCollection(joined.map(function(current) {
    var afterCollection = ee.ImageCollection.fromImages(current.get('after'));
    return afterCollection.map(function(image) {
      // The unmask is to prevent accumulation of masked pixels that
      // result from the undefined case of when either current or image
      // is masked.  It won't affect the sum, since it's unmasked to zero.
      return ee.Image(sign(current, image)).unmask(0);
    });
    // Set parallelScale to avoid User memory limit exceeded.
  }).flatten()).reduce('sum', 4);

  var mk_stat = mk.select('concordant_sum').subtract(mk.select('discordant_sum'));
  return mk_stat.toFloat()
}

/**
 * Generate a significance mask from the mann-kendall test results
 * @param {Image} img: Input image for which to generate a statistical significance mask.
 * @param {Image} mk_trend: Mann-Kendall S statistic image.
 * @returns {Image}: The Statistical significance mask for the 90%, 95% and 99% confidence intervals.
*/
exports.signifMask = function(img, mk_trend, period){
  // Define Kendall parameter values for a significance of 0.05
  //var period = end_year.get('year').subtract(start_year.get('year')).add(1);
  var kendall90 = ee.Number(_kendallCoefficient(period, 90));
  var kendall95 = ee.Number(_kendallCoefficient(period, 95));
  var kendall99 = ee.Number(_kendallCoefficient(period, 99));
  // Create final productivity trajectory output layer. Positive values are
  // significant increase, negative values are significant decrease.
  return ee.Image(-32768)
        .where(img.gt(0).and(mk_trend.abs().gte(kendall90)), 1)
        .where(img.gt(0).and(mk_trend.abs().gte(kendall95)), 2)
        .where(img.gt(0).and(mk_trend.abs().gte(kendall99)), 3)
        .where(img.lt(0).and(mk_trend.abs().gte(kendall90)), -1)
        .where(img.lt(0).and(mk_trend.abs().gte(kendall95)), -2)
        .where(img.lt(0).and(mk_trend.abs().gte(kendall99)), -3)
        .where(mk_trend.abs().lte(kendall90), 0);
}

/**
 * Hard-coded Kendall Coefficients as look-up dictionary
 * @param {Image} n: The kendall coefficient number to look up.
 * @param {Image} level: The confidence interval to use (90%, 95% or 99%).
 * @returns {Image}: The looked-up coefficient value.
*/
function _kendallCoefficient(n, level){
    // The minus 4 is because the indexing below for a sample size of 4
    n = n.subtract(4);
    var coefs = {90: ee.List([4, 6, 7, 9, 10, 12, 15, 17, 18, 22, 23, 27, 28, 32, 35, 37, 40, 42,
                  45, 49, 52, 56, 59, 61, 66, 68, 73, 75, 80, 84, 87, 91, 94, 98, 103,
                  107, 110, 114, 119, 123, 128, 132, 135, 141, 144, 150, 153, 159,
                  162, 168, 173, 177, 182, 186, 191, 197, 202]),
               95: ee.List([4, 6, 9, 11, 14, 16, 19, 21, 24, 26, 31, 33, 36, 40, 43, 47, 50, 54,
                    59, 63, 66, 70, 75, 79, 84, 88, 93, 97, 102, 106, 111, 115, 120,
                    126, 131, 137, 142, 146, 151, 157, 162, 168, 173, 179, 186, 190,
                    197, 203, 208, 214, 221, 227, 232, 240, 245, 251, 258]),
               99: ee.List([6, 8, 11, 18, 22, 25, 29, 34, 38, 41, 47, 50, 56, 61, 65, 70, 76, 81,
                    87, 92, 98, 105, 111, 116, 124, 129, 135, 142, 150, 155, 163, 170,
                    176, 183, 191, 198, 206, 213, 221, 228, 236, 245, 253, 260, 268,
                    277, 285, 294, 302, 311, 319, 328, 336, 345, 355, 364])}
    return coefs[level].get(n);
}


exports.stretchViz = function(viz_layer, viz_img, geom, property, options){

  var visPct = viz_img.reduceRegion({
    reducer: ee.Reducer.percentile([2,98]).setOutputs(['min','max']),
    geometry: geom,
    scale: 100,
    tileScale: 4
  });

  var viz_params = viz_layer.getVisParams();
  viz_params['min'] = visPct.getNumber(property+'_min').multiply(100).toInt16().divide(100);
  viz_params['max'] = visPct.getNumber(property+'_max').multiply(100).toInt16().divide(100);

  if (options.symmetry === 1) {
    viz_params['min'] = viz_params['min'].abs().max(viz_params['max'].abs()).multiply(-1);
    viz_params['max'] = viz_params['min'].abs().max(viz_params['max'].abs());
  }

  ee.Dictionary(viz_params).evaluate(function(params){
      viz_layer.setVisParams(params);
  });

  return [viz_layer, visPct.getNumber(property+'_min'), visPct.getNumber(property+'_max')]
};

/**
 * Load Annual Land Cover data from WorldCover, Dynamic World and GLAD Land Cover Change 2000-2020
 * @param {String} year: The year for which to load the land cover data
 * @param {FeatureCollection} geom: Geometry to stratify the land cover datasets
 * @returns {Array} landcover_data: an array containing World Cover-related data,
                                    Dynamic World-related data, and GLAD Land Cover Change data.
 * @ignore
 */
exports.landCoverDatasets = function(year, geom){

    // World Cover Dataset
    var world_cover = _worldCover(year);

    // Dynamic World Cover Dataset
    var dw = _dynamicWorld(year, geom);
    var dynamic_world = dw[0];
    var proba_hillshade = dw[1];

    // GLAD land cover change dataset 2000-2020
    var landmask = ee.Image("projects/glad/landBuffer4").mask();
    var LCLUC2020 = ee.Image('projects/glad/GLCLU2020/LCLUC_2020').updateMask(landmask);
    var wetland_cover = world_cover.eq(7).or(LCLUC2020.gte(100).and(LCLUC2020.lt(150))).selfMask();
    var tree_cover = LCLUC2020.gte(25).and(LCLUC2020.lte(100))
                     .or(LCLUC2020.gte(125).and(LCLUC2020.lte(200)))
                     .selfMask().rename('constant');

    var glad_lc_change = _gladLCChange(landmask);

    // Conflate the World Cover 2021 data with GLAD LULC2020 permanent water extent, which is more complete.
    world_cover = world_cover.where(LCLUC2020.gt(100).and(LCLUC2020.lt(150)), 6);

    var forest_change = ee.Image("UMD/hansen/global_forest_change_2020_v1_8");

    var landcover_data = [world_cover, tree_cover, wetland_cover,
                          dynamic_world, proba_hillshade, glad_lc_change, forest_change];

    return landcover_data
}

exports.makeLegend = function(title, palette, class_names, class_length){
  // Create a legend for the different crop types
  // set position of panel
  var legend = ui.Panel({
    style: {
      position: 'bottom-left',
      padding: '12px 15px'
    }
  });

  // Create legend title
  var legendTitle = ui.Label({
    value: title,
    style: {
      fontWeight: 'bold',
      fontSize: '18px',
      margin: '0 0 4px 0',
      padding: '0'
      }
  });

  legend.add(legendTitle);

  // Creates and styles 1 row of the legend.
  var makeRow = function(color, name) {
        // Create the label that is actually the colored box.
        var colorBox = ui.Label({
          style: {
            backgroundColor: color,
            // Use padding to give the box height and width.
            padding: '8px',
            fontSize: '12px',
            margin: '0 0 4px 0'
          }
        });

        // Create the label filled with the description text.
        var description = ui.Label({
          value: name,
          style: {margin: '0 0 4px 6px'}
        });

        // return the panel
        return ui.Panel({
          widgets: [colorBox, description],
          layout: ui.Panel.Layout.Flow('horizontal')
        });
  };

  // Add color and and names
  for (var i = 0; i <= class_length; i++) {
    legend.add(makeRow(palette[i], class_names[i]));
    }

  return legend
}

// Function to populate the color palette legends for the app layers
exports.populateLegend = function(legend_name, viz_params, add_char_min, add_char_max, options){

    // Create a legend for the different crop types
    // set position of panel
    var legend = ui.Panel({
      style: {
        position: 'bottom-left',
        padding: '12px 15px'
      }
    });

    // Create legend title
    var legend_title = ui.Label({
      value: legend_name,
      style: {
      fontWeight: 'bold',
      fontSize: '18px',
      margin: '0 0 0 0',
      padding: '0',
      //width: '115px'g
      }
      });

    legend.add(legend_title);

    // create the legend image
    var lon = ee.Image.pixelLonLat().select('latitude');
    var gradient = lon.multiply(ee.Number(viz_params.max).subtract(viz_params.min).divide(100)).add(viz_params.min);
    //var viz_palette =  viz_params['palette'].slice(0).reverse();
    //viz_params['palette'] = viz_palette;
    var legend_image = options.legend_image || gradient.visualize(viz_params);

    var max_label = ui.Label();
    // create text on top of legend
    var legend_panel_max = ui.Panel({
      widgets: [
      max_label
      //ui.Label(viz_params['max'] + add_char_max)
      ],
      });

    legend.add(legend_panel_max);

    // create thumbnail from the image
    var thumbnail = ui.Thumbnail({
      image: legend_image,
      params: {bbox: '0,0,10,100', dimensions:'10x25'},
      style: {padding: '1px', position: 'bottom-center', fontSize: '18px'}
      });

    legend.add(thumbnail);

    var min_label = ui.Label();
    // create text on top of legend
    var legend_panel_min = ui.Panel({
      widgets: [
      min_label
      ],
      });

    ee.Dictionary(viz_params).evaluate(function(params){
      max_label.setValue(params.max + add_char_max)
      min_label.setValue(params.min + add_char_min)
    });

    legend.add(legend_panel_min);

    return legend
};

/**
 * Nested function to load the World Cover data
 * @param {String} year: The year for which to load the land cover data
 * @returns {Image} world_cover: World Cover output
 * @ignore
*/
function _worldCover(year){
    var world_cover;
    if (year > '2020') {
        world_cover = ee.Image("ESA/WorldCover/v200/2021")
                      //.filterBounds(country.geometry()).mosaic()
                      // Remapping from world cover classes to a harmonized land cover nomenclature
                      .remap([10,20,30,40,50,60,70,80,90,95,100],
                             [3,4,5,2,1,8,9,6,7,3,5])
                      .rename('constant');
    } else {
        world_cover = ee.Image("ESA/WorldCover/v200/2020")
                      //.filterBounds(country.geometry()).mosaic()
                      // Remapping from world cover classes to a harmonized land cover nomenclature
                      .remap([10,20,30,40,50,60,70,80,90,95,100],
                             [3,4,5,2,1,8,9,6,7,3,5])
                      .rename('constant');
    }

    return world_cover
}

/**
 * Nested function to load the Dynamic World Data, aggregated for a specific year
 * @param {String} year: The year for which to load the land cover data
 * @param {FeatureCollection} geom: Geometry to stratify the land cover datasets
 * @returns {Image} world_cover: World Cover output
 * @ignore
*/
function _dynamicWorld(year, geom){
    // Alternative to WorldCover, which provides a yearly land cover dataset
    // (to be further tested whether it is applicable across geographies)
    // Temporal aggregation process tailored for the Sahelian context

    var dw_class_names = ['water', 'trees', 'grass', 'flooded_vegetation', 'crops',
                                    'shrub_and_scrub', 'built', 'bare', 'snow_and_ice'];

    var dw_col = ee.ImageCollection('GOOGLE/DYNAMICWORLD/V1')
                  .filterBounds(geom.geometry())
                  .filterDate(year+'-01-01', year+'-12-31T23:59:59');

    var dw_col_filt = dw_col.map(function(img){
      var img_mask = img.select(dw_class_names).updateMask(img.select(dw_class_names).gte(0.4));
      var img_mask_bs = img.select('bare').updateMask(img.select('bare').gte(0.6));
      img_mask = img_mask.addBands(img_mask_bs, ['bare'], true);
      return img_mask.addBands(img.select('label').updateMask(img_mask.reduce(ee.Reducer.max()).gt(0)));
    });

    // Get the most likely class probability.
    var top1Prob = dw_col.select(dw_class_names).reduce(ee.Reducer.mean(), 4).regexpRename('_mean', '');
    var top1Tree = dw_col_filt.select('trees').reduce(ee.Reducer.count(), 4).regexpRename('_count', '');
    var top1Grass = dw_col_filt.select('grass').reduce(ee.Reducer.count(), 4).regexpRename('_count', '');
    var top1Shrub = dw_col_filt.select('shrub_and_scrub').reduce(ee.Reducer.count(), 4).regexpRename('_count', '');
    var top1Label = dw_col_filt.select('label').reduce(ee.Reducer.mode(), 4).regexpRename('_mode', '');
    var top1All = dw_col.select('label').reduce(ee.Reducer.count(), 4).regexpRename('_count', '');

    // Create a hillshade of the most likely class probability on [0, 1];
    var proba_hillshade = ee.Terrain.hillshade(top1Prob.reduce(ee.Reducer.max()).multiply(100)
                                                 .reproject({crs:ee.Image(dw_col.first()).projection(), scale:10}))
                                                 .rename('proba');

    // Create an RGB image of the label (most likely class) on [0, 1].
    var dynamic_world = top1Label.remap([0, 1, 2, 3, 4, 5, 6, 7, 8], [6, 3, 5, 7, 2, 4, 1, 8, 9]).rename('label')
                       .where(top1Tree.divide(top1All).lt(0.6)
                              .and(top1Prob.select('crops').gt(top1Prob.select('shrub_and_scrub')))
                              .and(top1Prob.select('crops').gt(top1Prob.select('grass')))
                              .and(top1Prob.select('crops').gt(top1Prob.select('bare')))
                              .and(top1Prob.select('crops').gt(top1Prob.select('water')))
                              .and(top1Prob.select('crops').gt(top1Prob.select('flooded_vegetation')))
                              .and(top1Prob.select('crops').gt(top1Prob.select('built')))
                              .and(top1Prob.select('crops').gt(top1Prob.select('snow_and_ice'))), 2)
                       .where(top1Tree.divide(top1All).lt(0.6)
                              .and(top1Prob.select('crops').subtract(top1Prob.select('bare')).abs().lt(0.2))
                              .and(top1Prob.select('shrub_and_scrub').gt(top1Prob.select('grass')))
                              //.and(top1Prob.select('shrub_and_scrub').gt(top1Prob.select('crops'))))
                              .and(top1Shrub.gte(1))
                              .and(top1Grass.lt(top1Shrub)), 4)
                       .where(top1Tree.divide(top1All).lt(0.6)
                              .and(top1Prob.select('crops').subtract(top1Prob.select('bare')).abs().lt(0.2))
                              .and(top1Prob.select('grass').gt(top1Prob.select('shrub_and_scrub')))
                              //.and(top1Prob.select('grass').gt(top1Prob.select('crops')))
                              .and(top1Grass.gte(1))
                              .and(top1Shrub.lt(top1Grass)), 5)
                       .unmask(8);

    return [dynamic_world, proba_hillshade]
}

/**
 * Nested function to load the GLAD Land Cover Change
 * @param {String} landmask: A land surface mask to mask out the coastal water
 * @returns {Image} change_map: Reclassified GLAD Change data output
 * @ignore
*/
function _gladLCChange(landmask){

    // This data is used to refine land cover data, especially wrt the water extents
    var change = ee.Image('projects/glad/GLCLU2020/LCLUC').updateMask(landmask).rename('constant');

    var change_map = change.gt(50).and(change.lt(100)).multiply(7) // forest gain
                   .add(change.gt(150).and(change.lt(200)).multiply(8)) // wetland gain
                   .add(change.eq(240).multiply(1)) // forest --> short vegetation
                   .add(change.eq(247).multiply(4)) // short vegetation --> cropland
                   .add(change.eq(245).multiply(2)) // forest --> cropland
                   .add(change.eq(251).multiply(6)) // built-up gain
                   .add(change.eq(249).multiply(5)) // cropland --> short vegetation
                   .add(change.eq(246).multiply(3)) // wetland --> cropland
                   .add(change.eq(209).multiply(9)) // water loss
                   .add(change.eq(210).multiply(10)) // water gain
                   .add(change.eq(211).multiply(11)) // Ephemeral water
                   .selfMask()
                   .rename('constant');

    return change_map;
}

//Define function to extract time intervals to use to generate the temporal composites from Sentinel collections
exports.extractTimeRanges = function(start, end, agg_interval){
    /*
    Extract the time range data from the received time range and aggregation interval e.g.,
    input time interval: time_interval = ['2019-01-01','2020-01-01'], agg_interval: 60 days
    generate the following time intervals:
    time_range = [("2019-01-01T00:00:00Z", "2019-03-01T00:00:00Z"),
                ("2019-03-01T00:00:00Z", "2019-05-01T00:00:00Z"),
                ("2019-05-01T00:00:00Z", "2019-07-01T00:00:00Z"),
                ("2019-07-01T00:00:00Z", "2019-09-01T00:00:00Z"),
                ("2019-09-01T00:00:00Z", "2019-11-01T00:00:00Z"),
                ("2019-11-01T00:00:00Z", "2020-01-01T00:00:00Z")
    */

    var start_date = ee.Date(start);
    var end_date = ee.Date(end);

    // Number of intervals in the given "time_range" based on the specified "agg_interval" period
    var interval_no = ee.Date(end).difference(ee.Date(start), 'day').divide(agg_interval).round();
    var month_check = ee.Algorithms.If(ee.Number(30.4375 / agg_interval).round().gt(0),
                                       ee.Number(30.4375 / agg_interval).round(),
                                       ee.Number(1)); // The number of aggregation intervals within a month

    // Compute the relative date delta (in months) to add to each preceding period to compute the new one
    var rel_delta = ee.Number(end_date.difference(start_date, 'day'))
                    .divide(ee.Number(30.4375).multiply(interval_no)).ceil(); // 30.4375 days = average month length

    // Compute the first time interval end date by adding the relative date delta (in months) to the start date
    end_date = start_date.advance(start_date.advance(rel_delta, 'month')
                                  .difference(start_date, 'day')
                                  .divide(month_check), 'day')
                                  .advance(-1, 'second');

    var time_intervals = ee.List([ee.List([start_date, end_date])]);
    time_intervals = ee.List(ee.List.sequence(1, interval_no.subtract(1)).iterate(function(x,previous){
        start_date = ee.Date(ee.List(ee.List(previous).reverse().get(0)).get(1))
                     .advance(1, 'second'); //end_date of last element
        end_date = start_date
                   .advance(start_date.advance(rel_delta, 'month')
                   .difference(start_date, 'day')
                   .divide(month_check), 'day')
                   .advance(-1, 'second');

        return ee.List(previous).add(ee.List([start_date, end_date]));
    }, time_intervals));

    return time_intervals;
}

/**
 * This function adds a time band to the image
 * @param {Image} img: Input image for which to create time band
 * @returns {Image}: The original image conflated with a time band expressed in years
 * @ignore
*/
exports.createTimeBand = function(img) {
  // Scale milliseconds by a large constant to avoid very small slopes
  // in the linear regression output.
  return img.addBands(img.metadata('system:time_start').divide(3.154e10).add(1970));
};
