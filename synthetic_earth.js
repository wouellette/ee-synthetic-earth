// ****************************************************************************************************************** //
// ******************************** Implementation of synthetic earth functions ************************************* //
// ****************************************************************************************************************** //

// Import Dependencies
var utils = require('users/soilwatchtech/synthetic_earth:utils.js');

/**
 * Generate a Synthetic Composite Reflectance products (soil and vegetation)
 * @param {ImageCollection} collection: The input Image Collection
 * @param {ImageCollection} bs_collection: The input bare soil-masked Image Collection
 * @param {List} bands: List of bands to use for the synthetic reflectance composites
 * @param {List} agg_band: The band to use for aggregation
 * @param {List} algo: Bare soil pixel identification algorithm to use. options are 'GEOS3' or 'HISET'
 * @param {List} time_intervals: Nested list of date pairs ([DateStart0, DateEnd0], ..., [DateStartN, DateEndN])
                                 over which to compute the composites
 * @param {String} geom: The geometry to take into account for compute the histogram splitting approach (HISET)
 * @param {Boolean} dry_flag: Whether to generate the bare soil composite for the driest periods only.
                             It uses precipitation data to find time periods of absence of rainfall.
                             This may be useful if the bare soil composites are used for soil carbon assessment
                             that requires minimization of soil moisture conditions to increase reliability of modelling
 * @param {Integer} n_iter: Number of iterations to apply for the Weighted Geometric Median (default is 10).
 * @param {Integer} loss_func_param: Loss function parameter (defaults to -1) for the weighted geometric median.
 * @returns {ImageCollection} scr_ts: An image collection containing 4 temporal synthetic composites:
    - greenest earth: The small NDVI integral computed for the given time period (Eklundh & Jönsson, 2004)
      This is more representative than a median or max NDVI approach,
      as it incorporates both length and amplitude information of photosynthetically active part of season
    - barest earth: the barest earth condition computed using the weighted geometric median (Roberts et al., 2019)
    - bare earth: Bare soil composite with non-bare areas masked out (Heiden et al., 2022 for HISET
    or Dematte et al., 2019 for GEOS3)
 * @ignore
*/
exports.SCR = function(col, bs_col, bands, agg_band, algo, time_intervals, geom, options){

  var dry_flag = options.dry_flag || 0;
  var n_iter = options.n_iter || 10;
  var loss_func_param = options.loss_func_param || -1;
  var diff_tol = options.diff_tol || 7;

  var scr_ts = ee.ImageCollection(ee.List(time_intervals).map(function(time_interval){

    time_interval = ee.List(time_interval);

    var start_date = ee.Date(time_interval.get(0));
    var end_date = ee.Date(time_interval.get(1));
    // Create a intermediate timestamp for the aggregated images
    var mid_date = start_date.advance(start_date.difference(end_date, 'month').abs().divide(2), 'month').millis();

    // Compute the small integral for the time interval (Jönsson & Eklundh, 2004)
    var ndvi_image = exports.smallIntegral(col, bands, agg_band, start_date, end_date, mid_date);

    // Compute the bare soil composite from the bare soil collection
    var bs_image = exports.bareSoilComposite(bs_col, bands, ndvi_image.bandNames(), dry_flag, diff_tol, start_date, end_date, mid_date, geom);

    // Compute the weighted geometric median for the barest earth composite (Roberts et al., 2019)
    var barest_image = exports.weightedGeometricMedian(col, bands, agg_band, n_iter, loss_func_param, start_date, end_date);

    // Compile the greenest earth, bare earth and barest earth synthetic images into a single image
    return ndvi_image.addBands(bs_image)
                     .addBands(barest_image)
                     .set('system:time_start', mid_date);
    }))

    return scr_ts;
}

/**
 * Compute the Small integral of the NDVI temporal signal as per Jönsson & Eklundh, 2004
 * @param {ImageCollection} col: The input Image Collection
 * @param {List} bands: The list of bands to select for the aggregation
 * @param {List} agg_band: The band to use for aggregation
 * @param {Date} start_date: A timestamp corresponding to the start of the time interval for aggregation
 * @param {Date} end_date: A timestamp corresponding to the end of the time interval for aggregation
 * @param {Date} mid_date: A timestamp corresponding to the middle of the time interval for
                           which the small integral is calculated
 * @returns {Image} ndvi_image: An image containing the greenest earth composite,
                                i.e The small NDVI integral computed for the given time period (Eklundh & Jönsson, 2004)
 * @ignore
*/
exports.smallIntegral = function(col, bands, agg_band, start_date, end_date, mid_date){

    var ndvi_col = col.filterDate(start_date, end_date);

    // Compute 5th and 95th percentiles as pseudo min and max values.
    // This is done to avoid picking up outlier pixels (residual cloudy pixels or abnormally dark pixels)
    var ndvi_min_max = ndvi_col.select(agg_band).reduce(ee.Reducer.percentile([5, 95]), 4)
                             .set('system:time_start', mid_date);

    // Compute the amplitude of the time series
    var amplitude = ndvi_min_max.select(1).subtract(ndvi_min_max.select(0));

    // Generate the synthetic vegetation composite using the small integral calculation
    // (as per TIMESAT definition from Jönsson & Eklundh, 2004)
    // The calculation is a simplification of the actual integral calculation, as it simply takes the mean of the values
    // falling inside the small integral
    var ndvi_image = ndvi_col.map(function(img){
      return img.select(bands).updateMask(img.select(agg_band).gte(amplitude.divide(2).add(ndvi_min_max.select(0))));
    })
    .reduce(ee.Reducer.mean(), 4)
    .regexpRename('_mean', '_syvi')
    .set('system:time_start', mid_date);

    return ndvi_image
}

/**
 * Compute the Bare soil composite from a bare soil-masked collection using the geometric median
 * @param {ImageCollection} col: The input Image Collection
 * @param {List} bands: The list of bands to select for the aggregation
 * @param {List} band_names: Band names to rename the output bands to
 * @param {Boolean} dry_flag: Whether to generate the bare soil composite for the driest periods only.
                             It uses precipitation data to find time periods of absence of rainfall.
                             This may be useful if the bare soil composites are used for soil carbon assessment
                             that requires minimization of soil moisture conditions to increase reliability of modelling
 * @param {Integer} diff_tol: Number of days tolerance to take into account for the dry periods calculation.
 * @param {Date} start_date: A timestamp corresponding to the start of the time interval for aggregation
 * @param {Date} end_date: A timestamp corresponding to the end of the time interval for aggregation
 * @param {Date} mid_date: A timestamp corresponding to the middle of the time interval for aggregation
 * @param {Geometry} geom: The geometry of the area of interest to retrieve the dry periods
 * @returns {Image} bs_image: An image containing a bare soil composite, with non-bare areas masked out
                              (Heiden et al., 2022 for HISET or Dematte et al., 2019 for GEOS3)

 * @ignore
*/
exports.bareSoilComposite = function(col, bands, band_names, dry_flag, diff_tol, start_date, end_date, mid_date, geom){

    var bs_col = col.filterDate(start_date, end_date).map(utils.addDate).select(bands);

    // Calculate a list of dry period sprees (in days) in the course of the monitoring period longer than 7 consecutive days
    var dry_period_dates = utils.getDryPeriod(ee.ImageCollection("UCSB-CHG/CHIRPS/DAILY"),
                                              start_date, end_date, geom, {diff_tol: diff_tol});

    var bs_image;
    if (dry_flag === 1) {

        // if dry_flag is on, then calculate the bare soil geometric median of the all the bare pixels identified
        // within the detected dry sprees (take only the 3 longest dry sprees)
        var dry_filter = ee.Filter.or(ee.Filter.rangeContains('doy', dry_period_dates.get('start'),
                                                                     dry_period_dates.get('end')),
                                      ee.Filter.rangeContains('doy', dry_period_dates.get('start2'),
                                                                     dry_period_dates.get('end2')),
                                      ee.Filter.rangeContains('doy', dry_period_dates.get('start3'),
                                                                     dry_period_dates.get('end3')));

        // If there are fewer than 3 timestamps falling within the dry sprees periods, then fall back on
        // calculating the geometric median of all bare pixels in the time series for the monitoring period
        bs_image = ee.Image(ee.Algorithms.If(bs_col.filter(dry_filter).size().gte(3),
                                             bs_col.filter(dry_filter)
                                             .reduce(ee.Reducer.geometricMedian(bands.length), 4)
                                             .rename(band_names).regexpRename('syvi', 'sysi')
                                             .set('system:time_start', mid_date),
                                             bs_col.reduce(ee.Reducer.geometricMedian(bands.length), 4)
                                             .rename(band_names).regexpRename('syvi', 'sysi')
                                             .set('system:time_start', mid_date)
                                             )
                           );
    } else {
        // If dry_flag is off, calculate the geometric median for all pixels in the monitoring period
        bs_image = bs_col.reduce(ee.Reducer.geometricMedian(bands.length), 4)
                         .rename(band_names).regexpRename('syvi', 'sysi')
                         .set('system:time_start', mid_date);
    }

    return bs_image;
}

/**
 * Compute the Weighted Geometric Median on an image collection as per Roberts et al., 2019
 * @param {ImageCollection} col: The input Image Collection
 * @param {List} bands: The list of bands to select for the aggregation
 * @param {List} agg_band: The band to use for aggregation
 * @param {Integer} n_iter: The number of iterations to run to converge to the weighted geometric median
 * @param {Integer} loss_func_param: The loss function parameter for the Weighted Geometric Median computation (default -1)
 * @param {Date} start_date: A timestamp corresponding to the start of the time interval for aggregation
 * @param {Date} end_date: A timestamp corresponding to the end of the time interval for aggregation
 * @returns {Image} barest_image: An image containing the barest earth composite,
                                  i.e the barest earth condition as per Roberts et al., 2019
 * @ignore
*/
exports.weightedGeometricMedian = function(col, bands, agg_band, n_iter, loss_func_param, start_date, end_date){

  var bare_col = col.filterDate(start_date, end_date);

  // Method to calculate the Weighted Geometric Median according to Roberts et al., 2019
  // Multiply the NDVI values by the loss function parameter
  var barest_weight = bare_col.map(function(img){return img.addBands(img.select(agg_band)
                                                           .multiply(loss_func_param)
                                                           .rename('weights'))});
  // Sum up the weights
  var barest_weight_sum = barest_weight.select(['weights']).map(function(img){return img.exp()})
                                                           .reduce(ee.Reducer.sum());
  // Normalize the output
  var barest_weight_norm = barest_weight.map(function(img){return img.updateMask(img.select('weights').exp()
                                                                     .divide(barest_weight_sum))});
  // Calculate the mean of the normalized sum
  var barest_mean = barest_weight_norm.select(bands)
                                      .reduce(ee.Reducer.mean())
                                      .regexpRename('_mean', '_barest')
                                      .updateMask(ee.Image(1));

  // Iterate starting with the mean for n_iter times to converge to the weighted geometric median
  var weighted_gm = ee.Image(ee.List.sequence(0, n_iter).iterate(function(val, prev){
      prev = ee.Image(prev);
      var barest_weight_iter = barest_weight.map(function(img){
        var norm = img.select(bands).subtract(prev).pow(2).reduce(ee.Reducer.sum()).sqrt();
        var new_weights = img.select('weights').divide(norm.max(ee.Image(1e-6))).exp();
        return img.addBands(new_weights, ['weights'], true);
      });
      var barest_weight_sum = barest_weight_iter.select(['weights']).reduce(ee.Reducer.sum());
      barest_weight_iter = barest_weight_iter.map(function(img){
        return img.updateMask(img.select('weights').divide(barest_weight_sum))
      });

      return barest_weight_iter.select(bands)
                               .reduce(ee.Reducer.mean())
                               .regexpRename('_mean', '_barest')
                               .updateMask(ee.Image(1));
  }, barest_mean));

    return weighted_gm
}

/**
 * Performs the HISET approach to detect bare pixels in an collection of images. Method described in Heiden et al., 2022
 * @param {ImageCollection} col: The input Image Collection
 * @param {Image} lc_crop: A binary map indicating locations of cropland areas. Typically taken from a
                           recent high-resolution land cover data product like WorldCover or Dynamic World
 * @param {Image} lc_nat: A binary map indicating locations of other natural surfaces.
                          Typically, this can either be deciduous forest or grassland taken from a land cover product.
                          Heiden et al., 2022 yields improved results when using grassland rather than deciduous forest
 * @param {Image} lc_bu: A binary map indicating locations of built-up areas, typically taken from a
                         recent high-resolution land cover data product like WorldCover or Dynamic World
 * @param {Geometry} geom: The geometry of the area of interest to retrieve the histograms
                           The original paper retrieves one histogram per Sentinel-2 tile,
                           while here an histogram is generated for the area of interest
 * @param {Date} start_date: A timestamp corresponding to the start of the time interval for applying HISET
 * @param {Date} end_date: A timestamp corresponding to the end of the time interval for applying HISET
 * @param {Float} min_bucket_width: The minimum bucket width to set for the histogram splitting logic of HISET
                                    (default 0.0075)
 * @param {Integer} scale: The scale of the spatial reducer applied for the histogram computation (default 100m)
 * @returns {Array}: An array containing:
                     - the bare soil-masked image collection
                     - An active crop soil mask
                     - A permanent/fallow  soil mask
 * @ignore
*/
exports.addHISETMask = function(col, lc_crop, lc_nat, lc_bu, geom, options) {

  var start_date = options.start_date || 0;
  var end_date = options.end_date || 0;
  var min_bucket_width = options.min_bucket_width || 0.0075;
  var scale = options.scale || 100;

  // Filter collection temporally if start and end dates provided
  var col_red;
  if (start_date !== 0 && end_date !== 0){
    col_red = col.filterDate(start_date, end_date);
  } else {
    col_red = col;
  }

  // Compute the index of choice for computing HISET. Normalized Burn Ratio 2 has proven to yield best results
  // However, if NBR2 not computable because of absence of SWIR bands, fall back on NDVI
  var pv_min_max = col_red.map(function(img){
      var pv = ee.Image(ee.Algorithms.If(img.bandNames().contains('B11'),
                                         img.normalizedDifference(['B11', 'B12']),// Normalized Burn Ratio 2
                                         img.normalizedDifference(['B8', 'B4'])// Normalized Difference Vegetation Index
                                        )).rename('pv');
      return pv;
  })
  // Extract the minimum and maximum values using 5th and 95th percentile to avoid picking out noisy pixels
  .reduce(ee.Reducer.percentile([5, 95]), 4).rename(['pv_min', 'pv_max']);

  // Prepare the dictionary with the histogram inputs
  var lc_dict = ee.Dictionary({'t_max': ee.List([lc_bu, lc_crop, pv_min_max.select(1).rename('t_max')]),
                               't_min': ee.List([lc_crop, lc_nat, pv_min_max.select(0).rename('t_min')])
                });

  // Apply the HISET algorithm steps for minimum and maximum conditions
  var lc_dict_thresh = lc_dict.map(function(key, val){

      var lc1 = ee.List(val).get(0);
      var lc2 = ee.List(val).get(1);
      var pv_band = ee.Image(ee.List(val).get(2));

      // Compute histogram for the first land cover class
      var hist1 = ee.Array(
        pv_band.updateMask(lc1).reduceRegion({
        reducer: ee.Reducer.autoHistogram({minBucketWidth: min_bucket_width}),
        geometry: geom.geometry(),
        scale: scale,
        crs: 'EPSG:4326',
        maxPixels: 1e13,
        tileScale: 4
      }).get(key));

      // Compute histogram for the second land cover class
      var hist2 = ee.Array(
        pv_band.updateMask(lc2).reduceRegion({
        reducer: ee.Reducer.autoHistogram({minBucketWidth: min_bucket_width}),
        geometry: geom.geometry(),
        scale: scale,
        crs: 'EPSG:4326',
        maxPixels: 1e13,
        tileScale: 4
      }).get(key));

      var hist_range = hist1.slice(1, 0, 1).toList()
                       .cat(hist2.slice(1, 0, 1).toList())
                       .flatten().sort();
      var hmin = hist_range.get(0);
      var hmax = hist_range.get(-1);

      var hist1_count = hist1.slice(1, 1).reduce(ee.Reducer.sum(), [0]);
      var hist2_count = hist2.slice(1, 1).reduce(ee.Reducer.sum(), [0]);

      // Iteratively try to find the threshold by going over each histogram bin
      var hiset_thresh = ee.List.sequence(hmin, hmax, min_bucket_width).iterate(function(thresh, prev){

        var hiset = ee.Number(ee.List(prev).get(0));
        var thresh_min = ee.Number(ee.List(prev).get(1));

        var left_prop1 = hist1.slice(1, 1).mask(hist1.slice(1, 0, 1).lt(ee.Number(thresh)))
                                          .reduce(ee.Reducer.sum(), [0]).divide(hist1_count);
        var left_prop2 = hist2.slice(1, 1).mask(hist2.slice(1, 0, 1).lt(ee.Number(thresh)))
                                          .reduce(ee.Reducer.sum(), [0]).divide(hist2_count);
        var left_min = left_prop1.min(left_prop2).get([0, 0]);

        var right_prop1 = hist1.slice(1, 1).mask(hist1.slice(1, 0, 1).gt(ee.Number(thresh)))
                                           .reduce(ee.Reducer.sum(), [0]).divide(hist1_count);
        var right_prop2 = hist2.slice(1, 1).mask(hist2.slice(1, 0, 1).gt(ee.Number(thresh)))
                                           .reduce(ee.Reducer.sum(), [0]).divide(hist2_count);
        var right_min = right_prop1.min(right_prop2).get([0, 0]);

        thresh_min = ee.Number(ee.Algorithms.If(left_min.max(right_min).lt(hiset),
                                                thresh,
                                                thresh_min));
        hiset = ee.Number(ee.Algorithms.If(left_min.max(right_min).lt(hiset),
                                           left_min.max(right_min),
                                           hiset));

        return ee.List([hiset, thresh_min]);

      }, ee.List([10000, 0]));

      return hiset_thresh; // Return bare soil pixel stack
  });

  var mask_max = pv_min_max.select(1).gt(ee.Number(ee.List(lc_dict_thresh.get('t_max')).get(1)));
  var mask_min = pv_min_max.select(0).lt(ee.Number(ee.List(lc_dict_thresh.get('t_min')).get(1)));
  // Collateral active crop soil and permanent/fallow soil mask outputs
  var bs_crop_mask = mask_max.and(mask_min).selfMask().rename('bs_crop_mask');
  var bs_perm_mask = mask_max.not().and(mask_min).selfMask().rename('bs_perm_mask');

  // Go over the original Image Collection and apply the bare soil mask
  var col_masked = col_red.map(function(img){

      var pv = ee.Image(ee.Algorithms.If(img.bandNames().contains('B11'),
                                     img.normalizedDifference(['B11', 'B12']), // Normalized Burn Ratio 2
                                     img.normalizedDifference(['B8', 'B4']) // Normalized Difference Vegetation Index
                       )).rename('pv');

      return img.updateMask(pv.lt(ee.Number(ee.List(lc_dict_thresh.get('t_min')).get(1))))
                .addBands(img.select(0).neq(0).rename('nomask'));
  });

  return [col_masked, bs_crop_mask, bs_perm_mask];
};

/**
 * Performs the GEOS3 approach to detect bare pixels in a collection of images. Method described in Dematte et al., 2019
 * @param {Image} img: The input Image
 * @param {Boolean} rescale_flag: Whether to rescale the image to the [0,1] range or not
 * @param {Array} ndvi_thres: The upper and lower bounds of the NDVI threshold to use for the computation
 * @param {Array} nbr_thres: The upper and lower bounds of the NBR2 threshold to use for the computation
 * @param {Float} vnsir_thres: The upper bound of the VNSIR threshold to use for the computation
 * @returns {Image} geos3: An image where non-bare soil pixels are masked out
 * @ignore
*/
exports.addGEOS3Mask = function(img, options) {

  var rescale_flag = options.rescale_flag || 0;
  var ndvi_thres = options.ndvi_thres || [-0.25, 0.25];
  var nbr_thres = options.nbr_thres || [-0.3, 0.1];
  var vnsir_thres = options.vnsir_thres || 0.9;

  if (rescale_flag === 1) {
    img = img.divide(10000); // rescale to [0,1] reflectance.
  }

  var ndvi = img.normalizedDifference(['B8', 'B4']); // Normalized vegetation index
  var nbr2 = img.normalizedDifference(['B11', 'B12']); // Normalized Burn Ratio 2

  // Visible-toshortwave-infrared tendency index
  var vnsir = ee.Image(1)
              .subtract(ee.Image(2).multiply(img.select('B4'))
                                             .subtract(img.select('B3')).subtract(img.select('B2'))
              .add(ee.Image(3).multiply(img.select('B12').subtract(img.select('B8')))));

  // GEOS3 equation
  var geos3 = ndvi.gte(ndvi_thresh[0]).and(ndvi_thres[1])
              .and(nbr2.gte(nbr_thres[0]).and(nbr2.lte(nbr_thres[1])))
              .and(vnsir.lte(vnsir_thres)).rename('GEOS3');

  return geos3 // Return bare soil pixel stack
};

/**
 * Performs the NBR2 threshold algorithm approach to detect bare pixels in a collection of images.
   Method described in Dvorakova et al., 2023
 * @param {Image} img: The input Image
 * @param {Boolean} rescale_flag: Whether to rescale the image to the [0,1] range or not
 * @param {Array} nbr_thres: The upper bound of the NBR2 threshold to use for the computation
 * @returns {Image} geos3: An image where non-bare soil pixels are masked out
 * @ignore
*/
exports.addNormalizedNBR2Mask = function(img, options) {

  var rescale_flag = options.rescale_flag || 0;
  var nbr_thres = options.nbr_thres || 0.05;

  if (rescale_flag === 1) {
    img = img.divide(10000); // rescale to [0,1] reflectance.
  }

  var mean = img.reduce(ee.Reducer.sum(), 4).divide(img.bandNames().length()); // Mean Summed Reflectance
  var norm = img.divide(mean); // Normalized Reflectance

  var nbr2 = norm.normalizedDifference(['B11', 'B12']); // Normalized Burn Ratio 2

  return nbr2.lte(nbr_thres).rename('NBR2_norm'); // Return bare soil pixel stack
};

/**
 * Run the landtrendr algorithm on an image collection of annual images
 * @param {ImageCollection} col: Input Image Collection of annual images
 * @param {Integer} distDir: A factor to determine the direction of the indicators (-1 or 1)
 * @param {Dictionary} run_params: LandTrendr parameter dictionary
 * @returns {Image} bigDeltaImg: LandTrendr outputs for the detected change segment of greatest amplitude,
                                 including following bands:
                                 - 'yod': day of year band
                                 - 'endYr': end year band
                                 - 'startVal': time series value at first timestamp
                                 - 'endVal': time series value at end timestamp
                                 - 'mag': magnitude of the detected change segment
                                 - 'dur': duration of the detected change segment
                                 - 'rate': rate of change of the detected change segment
                                 - 'count': Number of changes detected
 * @ignore
*/
exports.landTrendr = function(col, distDir, run_params){

    var distDir =
    run_params.timeSeries = col; // add LT collection to the segmentation run parameter object
    // run LandTrendr spectral temporal segmentation algorithm
    var lt = ee.Algorithms.TemporalSegmentation.LandTrendr(run_params).select('LandTrendr');

    //------GET LandTrendR data
    var vertexMask = lt.arraySlice(0, 3, 4); // slice out the 'Is Vertex' row - yes(1)/no(0)
    var vertices = lt.arrayMask(vertexMask); // use the 'Is Vertex' row as a mask for all rows
    var left = vertices.arraySlice(1, 0, -1); // slice out the vertices as the start of segments
    var right = vertices.arraySlice(1, 1, null); // slice out the vertices as the end of segments
    var startYear = left.arraySlice(0, 0, 1); // get year dimension of LT data from the segment start vertices
    var startVal = left.arraySlice(0, 2, 3); // get spectral index dimension of LT data from the segment start vertices
    var endYear = right.arraySlice(0, 0, 1); // get year dimension of LT data from the segment end vertices
    var endVal = right.arraySlice(0, 2, 3); // get spectral index dimension of LT data from the segment end vertices

    // subtract the segment start year from the segment end year to calculate the duration of segments
    var dur = endYear.subtract(startYear);
    // substract the segment start index value from the segment end index value to calculate the delta of segments
    var mag = endVal.subtract(startVal);
    // calculate the rate of spectral change
    var rate = mag.divide(dur);

    var segInfo = ee.Image.cat([startYear.add(1), endYear, startVal, endVal, mag, dur, rate])
                  .toArray(0).mask(vertexMask.mask());
    // need to flip the delta here, since arraySort is working by ascending order
    var sortByThis = segInfo.arraySlice(0,4,5).toArray(0).multiply(-1);
    // sort the array by magnitude
    var segInfoSorted = segInfo.arraySort(sortByThis);
    // get the first segment in the sorted array (greatest magnitude vegetation loss segment)
    var bigDelta = segInfoSorted.arraySlice(1, 0, 1);

    // Convert the array Image into a multi-band image
    var bigDeltaImg = ee.Image.cat(bigDelta.arraySlice(0,0,1).arrayProject([1]).arrayFlatten([['yod']]),
       bigDelta.arraySlice(0,1,2).arrayProject([1]).arrayFlatten([['endYr']]),
       bigDelta.arraySlice(0,2,3).arrayProject([1]).arrayFlatten([['startVal']]).multiply(distDir),
       bigDelta.arraySlice(0,3,4).arrayProject([1]).arrayFlatten([['endVal']]).multiply(distDir),
       bigDelta.arraySlice(0,4,5).arrayProject([1]).arrayFlatten([['mag']]).multiply(distDir),
       bigDelta.arraySlice(0,5,6).arrayProject([1]).arrayFlatten([['dur']]),
       bigDelta.arraySlice(0,6,7).arrayProject([1]).arrayFlatten([['rate']]).multiply(distDir),
       segInfoSorted.arraySlice(0,4,5).arrayReduce(ee.Reducer.count(), [1]).arrayProject([1]).arrayFlatten([['count']]))

    return bigDeltaImg;
}
