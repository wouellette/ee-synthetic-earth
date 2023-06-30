// ****************************************************************************************************************** //
// ********************************* Synthetic Earth Example for Burkina Faso *************************************** //
// ****************************************************************************************************************** //

// Global parameters to define
// Masked out Sentinel-2 and Landsat parameters, but these could be used instead of nicfi to run the same analysis
var START_YEAR = '2016'; // Start year for the Dynamic World Dataset
var END_YEAR = '2022'; // End year for the derivation of temporal indicators and current year indicators
//var S2_BANDS = ['B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8', 'B8A', 'B11', 'B12']; // S2 bands to use in the workflow
//var S2_COLLECTION = 'COPERNICUS/S2_SR_HARMONIZED'; // S2 collection to use in the workflow
//var LS_BANDS = ['B2', 'B3', 'B4', 'B8', 'B11', 'B12']; // Landsat bands to use in the workflow
var NICFI_BANDS = ['B', 'G', 'R', 'N'];
var SURFACE_WATER = ee.Image("JRC/GSW1_2/GlobalSurfaceWater").select('max_extent').eq(1); // Surface Water Mask
//var CLD_PRB_THRESH = 20; // Cloud probability threshold to mask clouds. 40% is the default value of s2cloudless
//var CLD_FILTER = 40; // Threshold on sentinel-2 Metadata field determining whether cloud pixel percentage in image
//var NIR_DRK_THRESH = 0.15; // A threshold that determines when to consider a dark area a cloud shadow or not
//var CLD_PRJ_DIST = 10; // The distance (in no of pixels) in which to search from detected cloud to find cloud shadows
//var CLD_BUFFER = 200; // The cloud buffer (in meters) to use around detected cloud pixels to mask additionally
//var SNW_THRESH = 0.4; // The snow buffer (in meters) to use around detected cloud pixels to mask additionally
//var MASK_RES = 60; // resolution at which to generate and apply the cloud/shadow mask. 60m instead of 10m to speed up
var AGG_INTERVAL = 365; // Number of days for which to compute temporal aggregates
var ALGO = 'HISET'; // Bare soil retrieval algorithm to use. Either HISET (Heiden et al., 2022) or GEOS3 (Dematte et al., 2019)
// Builtup Mask. Default is World Settlement Footprint 2019.
var BUILTUP_DATA = ee.ImageCollection("projects/sat-io/open-datasets/WSF/WSF_2019").mosaic().unmask(0).eq(255);

// List of Countries for each continent where NICFI Basemaps are available
var africa_list_nicfi = ['Abyei', 'Angola', 'Benin', 'Botswana', 'Burkina Faso', 'Burundi', 'Cameroon',
                   'Cape Verde', 'Central African Republic', 'Chad', 'Comoros', 'Congo', "C te d'Ivoire",
                   "Côte d'Ivoire", 'Democratic Republic of the Congo', 'Djibouti', 'Eritrea', 'Ethiopia',
                   'Equatorial Guinea', 'Gabon', 'Gambia', 'Ghana', 'Guinea', 'Guinea-Bissau', 'Kenya', 'Lesotho',
                   'Liberia', 'Madagascar', 'Malawi', 'Mali', 'Mauritius', 'R union', "Réunion",
                   'Mozambique', 'Namibia', 'Niger', 'Nigeria', 'Rwanda', 'Sao Tome and Principe', 'Senegal',
                   'Seychelles', 'Sierra Leone', 'Somalia', 'South Africa', 'South Sudan', 'Sudan', 'Swaziland',
                   'Togo', 'Uganda', 'United Republic of Tanzania', 'Zambia', 'Zimbabwe'];

var america_list_nicfi = ['Mexico', 'Honduras', 'Guatemala', 'Belize', 'El Salvador', 'Nicaragua', 'Costa Rica',
                          'Panama', 'Colombia', 'Cuba', 'Bahamas', 'Venezuela', 'Guyana', 'French Guiana', 'Ecuador',
                          'Suriname', 'Uruguay', 'Paraguay', 'Bermuda', 'Anguilla', 'Netherlands Antilles', 'Peru',
                          'Bolivia', 'Argentina', 'Chile', 'Brazil', 'United States Virgin Islands',
                          'Falkland Islands (Malvinas)', 'Pitcairn', 'Turks and Caicos islands',
                          'Montserrat', 'Cayman Islands', 'Aruba', 'Bird Island', 'Puerto Rico', 'Navassa Island',
                          'South Georgia and the South Sandwich Islands', 'Bouvet Island', 'Haiti', 'Barbados',
                          'Saint Vincent and the Grenadines', 'Saint Lucia', 'Trinidad and Tobago',
                          'Saint Kitts and Nevis', 'Jamaica', 'Grenada', 'Antigua and Barbuda', 'Dominican Republic',
                          'Dominica', 'Guadeloupe', 'Martinique', 'Heard Island and McDonald Islands'];

var asia_list_nicfi = ['Papua New Guinea', 'Indonesia', 'Malaysia', 'Brunei Darussalam', 'Vanuatu', 'Singapore',
                       'Thailand', "Lao People's Democratic Republic", 'Viet Nam', 'Cambodia', 'Myanmar', 'Philippines',
                       'Timor-Leste', 'Sri Lanka', 'India', 'Bangladesh', 'Bhutan', 'Sri Lanka', 'Maldives',
                       'French Polynesia', 'Wallis and Futuna', 'Nepal', 'British Indian Ocean Territory',
                       'Cocos (Keeling) Islands', 'Christmas Island', 'Jammu and Kashmir', 'Arunachal Pradesh',
                       'Aksai Chin', 'China/India', 'Spratly Islands', 'Paracel Islands', 'Palau', 'Guam',
                       'Northern Mariana Islands', 'Solomon Islands', 'New Caledonia', 'Fiji', 'Tonga',
                       'American Samoa', 'Tuvalu', 'Tokelau', 'Niue', 'Cook Islands', 'Kiribati', 'Nauru',
                       'Micronesia', 'Marshall Islands'];

// Specify country for which to load NICFI data
var adm0_name = 'Burkina Faso';
// Select an ADMIN1 region of the country
var AOI = ee.FeatureCollection("projects/sat-io/open-datasets/geoboundaries/CGAZ_ADM1")
          .filter(ee.Filter.eq('shapeName', 'Centre-Nord'));

// Load AOI and zoom in on it
Map.centerObject(AOI);
Map.addLayer(AOI, {}, 'AOI');

// Dependencies Import
var palettes = require('users/gena/packages:palettes');
var utils = require('users/soilwatchtech/synthetic_earth:utils.js');
var syntheticEarth = require('users/soilwatchtech/synthetic_earth:synthetic_earth.js');

// Load latest available land cover datasets (World Cover, Dynamic World)
var landcover_data = utils.landCoverDatasets(END_YEAR, AOI);
var land_cover = landcover_data[0];
var tree_cover = landcover_data[1];
var wetland_cover = landcover_data[2];
var dynamic_world = landcover_data[3];
var proba_hillshade = landcover_data[4];
var glad_change = landcover_data[5];
var forest_change = landcover_data[6];

// Load NICFI data for the right continent depending on your area of interest
if (africa_list_nicfi.indexOf(adm0_name) >= 0) {
    var nicfi_col = ee.ImageCollection("projects/planet-nicfi/assets/basemaps/africa")
      .filterBounds(AOI.geometry())
      .map(function(img){return img.divide(10000).addBands(img.normalizedDifference(['N', 'R'])
                                   .rename('NDVI'))
                                   .rename(['B2', 'B3', 'B4', 'B8', 'NDVI'])
                                   .copyProperties(img, ['system:time_start'])

      });
} else if (america_list_nicfi.indexOf(adm0_name) >= 0) {
    var nicfi_col = ee.ImageCollection("projects/planet-nicfi/assets/basemaps/americas")
      .filterBounds(AOI.geometry())
      .map(function(img){return img.divide(10000).addBands(img.normalizedDifference(['N', 'R'])
                                   .rename('NDVI'))
                                   .rename(['B2', 'B3', 'B4', 'B8', 'NDVI'])
                                   .copyProperties(img, ['system:time_start'])
      });
} else if (asia_list_nicfi.indexOf(adm0_name) >= 0) {
    var nicfi_col = ee.ImageCollection("projects/planet-nicfi/assets/basemaps/asia")
      .filterBounds(AOI.geometry())
      .map(function(img){return img.divide(10000).addBands(img.normalizedDifference(['N', 'R'])
                                   .rename('NDVI'))
                                   .rename(['B2', 'B3', 'B4', 'B8', 'NDVI'])
                                   .copyProperties(img, ['system:time_start'])
      });
}

// Load Sentinel-2 collection if you prefer to run it on Sentinel-2 than NICFI Basemaps
/*var s2_col = utils.S2CloudMasked(S2_COLLECTION, S2_BANDS, START_YEAR, END_YEAR, SURFACE_WATER,
                                        CLD_FILTER, NIR_DRK_THRESH, CLD_PRJ_DIST, CLD_PRB_THRESH, CLD_BUFFER,
                                        SNW_THRESH, MASK_RES, AOI)
                                        .map(utils.applyBRDF('s2'));*/

// Apply bare soil masking algorithm, whether HISET or a simply NDVI threshold
if (ALGO == 'HISET'){
    var nicfi_bs_outputs = syntheticEarth.addHISETMask(nicfi_col,
                                                       land_cover.eq(2), // Use crop class
                                                       // Use natural surfaces classes, i.e. shrub and scrub and grass
                                                       land_cover.eq(4).or(land_cover.eq(5)),
                                                       BUILTUP_DATA.or(land_cover.eq(1)), // Use builtup data
                                                       AOI,
                                                       {start_year: START_YEAR+'-01-01', end_year: END_YEAR+'-12-31',
                                                        scale: 100, min_bucket_width: 0.0075});
    var nicfi_bs_col = nicfi_bs_outputs[0];
    var nicfi_bs_crop_mask = nicfi_bs_outputs[1];
    var nicfi_bs_perm_mask = nicfi_bs_outputs[2];
} else {
    nicfi_bs_col = nicfi.map(function(img){
        var img_rs = img.divide(10000); // rescale to [0,1] reflectance.
        var ndvi = img_rs.normalizedDifference(['B8', 'B4']); // Normalized vegetation index

        return img.addBands(ndvi.rename('NDVI')).updateMask(ndvi.gte(-0.25).bitwiseAnd(ndvi.lte(0.25)));
    });
    var nicfi_bs_crop_mask = ee.Image(0);
    var nicfi_bs_perm_mask = ee.Image(0);
}

// Generate time intervals of 1 year to extract temporal aggregates for each year
var time_intervals = utils.extractTimeRanges(START_YEAR+'-01-01', END_YEAR+'-12-31', AGG_INTERVAL);

// Synthetic Composite Reflectance for greenest earth (SYVI), barest earth and bare soil (SYSI)
var nicfi_scr = syntheticEarth.SCR(nicfi_col, nicfi_bs_col, ['B2', 'B3', 'B4', 'B8'], 'NDVI', ALGO, time_intervals, AOI,
                                   {n_iter: 10, loss_func_param: -1, dry_flag: 0, diff_tol: 0})
                .map(function(img){
                    var brightness_index = img.select('B4_barest').pow(2).add(img.select('B3_barest').pow(2))
                                           .divide(2).sqrt().rename('bi');
                    return img.addBands(brightness_index).multiply(10000).toInt16().copyProperties(img, ['system:time_start'])
                })
                .map(utils.createTimeBand);

// Extract the latest output (i.e latest year) and the respective greenest earth, bare soil and barest earth aggregates
var scr_bright = nicfi_scr.sort('system:time_start', false).first();
var ndvi_image_nicfi = scr_bright.select('.*_syvi');
var bs_image_nicfi = scr_bright.select('.*_sysi');
var barest_image_nicfi = scr_bright.select('.*_barest');

// Same bare soil mask applied with the dry_flag, i.e. will take bare soil pixels correspond to dry periods (i.e. no rain)
var nicfi_scr_dry = syntheticEarth.SCR(nicfi_col, nicfi_bs_col, ['B2', 'B3', 'B4', 'B8'], 'NDVI', ALGO, time_intervals, AOI,
                                       {n_iter: 10, loss_func_param: -1, dry_flag: 1, diff_tol: 0})
                    .map(function(img){
                        var brightness_index = img.select('B4_barest').pow(2).add(img.select('B3_barest').pow(2))
                                               .divide(2).sqrt().rename('bi');
                        return img.addBands(brightness_index).multiply(10000).toInt16().copyProperties(img, ['system:time_start'])
                    })
                    .map(utils.createTimeBand);

// Extract the latest output (i.e latest year) and the respective greenest earth, bare soil and barest earth aggregates
var scr_bright_dry = nicfi_scr_dry.sort('system:time_start', false).first();
var ndvi_image_nicfi_dry = scr_bright_dry.select('.*_syvi');
var bs_image_nicfi_dry = scr_bright_dry.select('.*_sysi');
var barest_image_nicfi_dry = scr_bright_dry.select('.*_barest');

// Mask out the areas of the bare soil aggregate with the crop soil and permanent soil masks
nicfi_bs_crop_mask = nicfi_bs_crop_mask.updateMask(bs_image_nicfi.select(0).neq(0));
nicfi_bs_perm_mask = nicfi_bs_perm_mask.updateMask(bs_image_nicfi.select(0).neq(0));
bs_image_nicfi = bs_image_nicfi.updateMask(nicfi_bs_crop_mask.unmask(0).or(nicfi_bs_perm_mask.unmask(0)));
bs_image_nicfi_dry = bs_image_nicfi_dry.updateMask(nicfi_bs_crop_mask.unmask(0).or(nicfi_bs_perm_mask.unmask(0)));

// Load data in map viewer
Map.addLayer(ndvi_image_nicfi.clip(AOI.geometry()),
             {bands:['B4_syvi', 'B3_syvi', 'B2_syvi'], min: 0, max:2000},
             'PlanetScope NICFI SYVI RGB' + END_YEAR);
Map.addLayer(ndvi_image_nicfi.clip(AOI.geometry()),
             {bands:['B8_syvi', 'B4_syvi', 'B3_syvi'], min: 0, max:4000},
             'PlanetScope NICFI SYVI False RGB' + END_YEAR);
Map.addLayer(bs_image_nicfi.clip(AOI.geometry()),
             {bands:['B4_sysi', 'B3_sysi', 'B2_sysi'], min: 0, max:2000},
             'PlanetScope NICFI SYSI RGB' + END_YEAR);
Map.addLayer(bs_image_nicfi_dry.clip(AOI.geometry()),
             {bands:['B4_sysi', 'B3_sysi', 'B2_sysi'], min: 0, max:2000},
             'PlanetScope NICFI DRY SYSI RGB' + END_YEAR);
Map.addLayer(bs_image_nicfi.clip(AOI.geometry()),
             {bands:['B8_sysi', 'B4_sysi', 'B3_sysi'], min: 0, max:4000},
             'PlanetScope NICFI SYSI False RGB' + END_YEAR);
Map.addLayer(barest_image_nicfi.clip(AOI.geometry()),
             {bands:['B4_barest', 'B3_barest', 'B2_barest'], min: 0, max:2000},
             'PlanetScope NICFI barest earth RGB' + END_YEAR);
Map.addLayer(barest_image_nicfi.clip(AOI.geometry()),
             {bands:['B8_barest', 'B4_barest', 'B3_barest'], min: 0, max:4000},
             'PlanetScope NICFI barest earth False RGB' + END_YEAR);
Map.addLayer(nicfi_bs_crop_mask.unmask(0).add(nicfi_bs_perm_mask.unmask(0).multiply(2))
                              .selfMask().clip(AOI.geometry()),
             {palette: ['pink', 'brown'], min: 1, max: 2},
             'NICFI Bare Soil Mask');

// Generate a pretty multi-temporal NICFI basemap composite
var nicfi_mt = nicfi_col.filterDate(END_YEAR+'-01-01', END_YEAR+'-12-31T23:59:59')
                    .sort('system:time_start', false).qualityMosaic('NDVI').clip(AOI.geometry())
                    .select('NDVI')
     .addBands(nicfi_col.filterDate(ee.String(ee.Number.parse(END_YEAR).subtract(1)).cat('-01-01'),
                                ee.String(ee.Number.parse(END_YEAR).subtract(1)).cat('-12-31T23:59:59'))
                    .sort('system:time_start', false).qualityMosaic('NDVI').clip(AOI.geometry())
                    .select('NDVI')
                    )
     .addBands(nicfi_col.filterDate(ee.String(ee.Number.parse(END_YEAR).subtract(2)).cat('-01-01'),
                                ee.String(ee.Number.parse(END_YEAR).subtract(2)).cat('-12-31T23:59:59'))
                    .sort('system:time_start', false).qualityMosaic('NDVI').clip(AOI.geometry())
                    .select('NDVI')
                    )
     .multiply(10000).toInt16();

Map.addLayer(nicfi_mt,
             {bands:['NDVI', 'NDVI_1', 'NDVI_2'], min: 0, max:7000},
             'PlanetScope NICFI NDVI composite (R:' + END_YEAR +
             ', G: ' + ee.String(ee.Number.parse(END_YEAR).subtract(1)).getInfo() +
             ', B: ' + ee.String(ee.Number.parse(END_YEAR).subtract(2)).getInfo() + ')');

// Compute the Mann-Kendall Significance test for the NICFI bare soil brightness composite time series
var mk_trend_scr_nicfi = utils.mannKendall(nicfi_scr_dry.select('bi')).rename('mannkendall');
// Compute the Temporal trend of bare soil brightness from the computed bare soil brightness index ('bi')
var scr_nicfi_ols = nicfi_scr_dry.select(['system:time_start', 'bi'])
                            .reduce(ee.Reducer.linearFit(), 4).select('scale').addBands(mk_trend_scr_nicfi);
// Generate a significance mask from the mann-kendall test results
var signif_scr_nicfi = utils.signifMask(scr_nicfi_ols.select('scale'), mk_trend_scr_nicfi, nicfi_scr_dry.size())
                       .rename('significance');
// Add significance mask to trend outputs
scr_nicfi_ols = scr_nicfi_ols.addBands(signif_scr_nicfi);

// Load temporal color palettes
var landtrendr_palette = palettes.colorbrewer.RdYlGn[11].slice(0);
var landtrendr_palette_rev = palettes.colorbrewer.RdYlGn[9].slice(0).reverse();

// Load brightness trend masked for bare soil pixels only, and with the mann-kendall significance for the 95% confidence interval
var nicfi_brightness_layer = ui.Map.Layer(scr_nicfi_ols.select('scale')
                                          .updateMask(signif_scr_nicfi.neq(0)
                                                      .and(signif_scr_nicfi.abs().gte(2))
                                                      .and(bs_image_nicfi_dry.select(0).gt(0)))
                                          .clip(AOI.geometry()),
                                          {min: -200, max: 200, palette: landtrendr_palette_rev},
                                          'NICFI Barest Brightness trend');
Map.add(nicfi_brightness_layer);

var ltParams = {maxSegments:            6,
              spikeThreshold:         0.9,
              vertexCountOvershoot:   3,
              preventOneYearRecovery: true,
              recoveryThreshold:      0.25,
              pvalThreshold:          0.05,
              bestModelProportion:    0.75,
              minObservationsNeeded:  3};

// Filter the composite time series for SYVI bands (greenest earth) only
nicfi_scr_dry = nicfi_scr_dry.map(function(img){
                return img.addBands(img.normalizedDifference(['B8_syvi', 'B4_syvi'])
                                    .multiply(10000).toShort().rename('NDVI_syvi'))
                       .select(['NDVI_syvi', 'B8_syvi', 'B4_syvi', 'B3_syvi', 'B2_syvi'])
                       .copyProperties(img, ['system:time_start'])
            });

// Compute the Mann-Kendall Significance test for the NICFI greenest earth composite time series
var mk_trend_landtrendr_nicfi = utils.mannKendall(nicfi_scr_dry.select('NDVI_syvi')).rename('mannkendall');
// Compute the landtrendr algorithm from the greenest earth composite
var landtrendr_nicfi = syntheticEarth.landTrendr(nicfi_scr_dry, 1, ltParams).addBands(mk_trend_landtrendr_nicfi);
// Generate a significance mask from the mann-kendall test results
var signif_landtrendr_nicfi = utils.signifMask(landtrendr_nicfi.select('mag'), mk_trend_landtrendr_nicfi, nicfi_scr_dry.size())
                       .rename('significance');
// Add significance mask to trend outputs
landtrendr_nicfi = landtrendr_nicfi.addBands(signif_landtrendr_nicfi);

// Generate legend to interpret the trend information
var landtrendr_legend = utils.populateLegend('Landtrendr/BI trends',
                                       {min: -1000, max: 1000, palette: landtrendr_palette},
                                       " -", " +", {});
Map.add(landtrendr_legend);

var bs_mask_legend = utils.makeLegend('Bare Soil Mask', ['pink', 'brown'], ['active crop soil', 'permanent/fallow soil'], 1);
Map.add(bs_mask_legend);

// Load landtrendr magnitude with the mann-kendall significance for the 95% confidence interval
var nicfi_landtrendr_layer = ui.Map.Layer(landtrendr_nicfi.updateMask(signif_landtrendr_nicfi.neq(0)
                                                          .and(signif_landtrendr_nicfi.abs().gte(2)))
                                                          .clip(AOI.geometry()),
                                          {min: -1000, max: 1000, bands: ['mag'], palette: landtrendr_palette},
                                          'NICFI SYVI Landtrendr Magnitude of Change');
Map.add(nicfi_landtrendr_layer);

// Export generated composite results.
Export.image.toAsset({
  image: ndvi_image_nicfi_dry.clip(AOI.geometry()),
  description:'greenest_earth_2022',
  assetId: 'users/soilwatchtech/BurkinaFaso/greenest_earth_2022',
  region: AOI.geometry(),
  crs: 'EPSG:4326',
  scale: 4.5,
  maxPixels:1e13
});

// Export generated composite results.
Export.image.toAsset({
  image: bs_image_nicfi_dry.clip(AOI.geometry()),
  description:'bare_soil_2022',
  assetId: 'users/soilwatchtech/BurkinaFaso/bare_soil_2022',
  region: AOI.geometry(),
  crs: 'EPSG:4326',
  scale: 4.5,
  maxPixels:1e13
});

// Export generated composite results.
Export.image.toAsset({
  image: barest_image_nicfi_dry.clip(AOI.geometry()),
  description:'barest_earth_2022',
  assetId: 'users/soilwatchtech/BurkinaFaso/barest_earth_2022',
  region: AOI.geometry(),
  crs: 'EPSG:4326',
  scale: 4.5,
  maxPixels:1e13
});

// Export generated composite results.
Export.image.toAsset({
  image: nicfi_bs_crop_mask.unmask(0).add(nicfi_bs_perm_mask.unmask(0).multiply(2))
        .selfMask().clip(AOI.geometry()),
  description:'soil_mask_2022',
  assetId: 'users/soilwatchtech/BurkinaFaso/soil_mask_2022',
  region: AOI.geometry(),
  crs: 'EPSG:4326',
  scale: 4.5,
  maxPixels:1e13
});

// Export generated composite results.
Export.image.toAsset({
  image: scr_nicfi_ols.clip(AOI.geometry()),
  description:'bi_trend_2022',
  assetId: 'users/soilwatchtech/BurkinaFaso/bi_trend_2016_2022',
  region: AOI.geometry(),
  crs: 'EPSG:4326',
  scale: 4.5,
  maxPixels:1e13
});

// Export generated composite results.
Export.image.toAsset({
  image: landtrendr_nicfi.clip(AOI.geometry()),
  description:'landtrendr_2016_2022',
  assetId: 'users/soilwatchtech/BurkinaFaso/landtrendr_2016_2022',
  region: AOI.geometry(),
  crs: 'EPSG:4326',
  scale: 4.5,
  maxPixels:1e13
});