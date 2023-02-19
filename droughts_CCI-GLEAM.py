# Input ..................................................................
# Directories
path_ccism = './data/inputs/ESACCI-SOILMOISTURE-L3S-SSMV-COMBINED-fv06.1.nc'   # ESA CCI dataset of soil moisture
path_gleam = './data/inputs/SMroot_GLEAM_v3.5a.nc'                             # GLEAM SMroot dataset

# Parameters for calculation
limit_pearson = 0.4            # Minimum value of Pearson coefficient to merge data
t_frame = [1990, 2010]         # Period of interest [first_year, last_year]
upscale = True                 # Spatial resolution of outcome. If True : 30arcmin, False : 15arcmin
perc_vtlm = 0.2                # Percentile that defines the threshold value to identify a drought in dates with runoff for the Variable Threshold Level method
perc_cdpm = 0.8                # Percentile that defines the threshold value to identify a drought in dates without runoff for the Consecutive Dry Period method
n_cores = 50                   # Number of cores to be used during the multiprocessing

# Output
path_out  = './data/dhe/'

# Packages ...............................................................
import xarray, pandas, datetime, sys, numpy, extremes, multiprocess
pandas.options.mode.chained_assignment = None

# Functions ..............................................................
'''Based on functions in extremes.py'''
# Importing NetCDF
def import_dataset(ds_path, upscale):
    # Importing NetCDF and defining raster characteristics
    xa = xarray.open_dataset(ds_path, decode_times=False)
    var = list(xa.keys())[0]
    if upscale == True: xa = xa.coarsen(lon=2, lat=2).mean()
    
    # Defining space related properties
    lons, lats = xa['lon'].values, xa['lat'].values
    coords = [str((lon,lat)) for lat in lats for lon in lons]
    
    # Defining time related properties
    times = xa['time'].values
    dates = pandas.to_datetime(xarray.decode_cf(xa).time.values)
    
    # Converting xarray information to 2D array and data frame
    da = xa[var].values
    ar = da.reshape(da.shape[0], da.shape[1] * da.shape[2])
    df = pandas.DataFrame(ar, index=dates, columns=coords)
    
    # Dropping empty pixels
    df.dropna(how='all', axis=1, inplace=True)
    
    return df, lons, lats

# Exporting NetCDF
def export_dataset(xa, var, lats, lons, dates, perc, path_out, t_frame, limit_pearson):
    # Defining time properties
    year, month, day = dates[0].year, dates[0].month, dates[0].day
    date_ref = dates[0].toordinal()
    times = [t.toordinal()-date_ref for t in dates]
    
    # Setting-up dataset
    ar = xarray.DataArray(
        data = xa,
        dims = ('time', 'lat', 'lon'),
        coords = {'time': times, 'lat': lats, 'lon': lons})
    
    # Assigning attributes
    ar.time.attrs['standard_name'] = 'time'
    ar.time.attrs['units'] = f'days since {year}-{month}-{day}'
    ar.time.attrs['calendar'] = 'standard'
    ar.lat.attrs['standard_name'] = 'latitude'
    ar.lat.attrs['units'] = 'degrees_north'
    ar.lon.attrs['standard_name'] = 'longitude'
    ar.lon.attrs['units'] = 'degrees_east'
    
    # Converting to dataset
    ds = xarray.Dataset({var: ar})
    
    # Exporting dataset
    res = int(abs(lons[0]-lons[1]) * 60)
    perc = int(perc*100)
    encode = {var: {'dtype':'float32', 'zlib':True, 'complevel':9}}
    ds.to_netcdf(f'{path_out}\\{var}_CCI-GLEAM_{res}arcmin_monthly_{t_frame[0]}-{t_frame[1]}_p={perc}%_pearson={int(limit_pearson*100)}%.nc', encoding=encode)

# Pre-processing .........................................................
print('> Pre-processing started...')

# Opening datasets 
df_cci, lons, lats = import_dataset(path_ccism, upscale)
print("  '-> CCI dataset loaded")
df_glm,  _  ,  _   = import_dataset(path_gleam, upscale)
print("  '-> GLEAM dataset loaded")

# Filtering list of common pixels
pixels = list(set(df_cci.columns.tolist()) & set(df_glm.columns.tolist()))
print("  '-> workable pixels selected")

# Creating empty data arrays
dates = pandas.date_range(f'{t_frame[0]}-01-01', f'{t_frame[1]}-12-31', freq='MS')
df = numpy.full([len(dates), len(pixels)], numpy.nan)
sm = numpy.full([len(dates), len(lats), len(lons)], numpy.nan)
x, y, pxls = lons.tolist(), lats.tolist(), []

# Combining datasets: filling CCIsm gaps with GLEAMsmroot data
for pixel in pixels:
    # Creating temporal dataframe with the two timeseries of interest
    ts1, ts2 = df_cci.loc[:,pixel], df_glm.loc[:,pixel]
    df_tmp = pandas.concat([ts1,ts2], axis=1)
    df_tmp.columns = ['CCI','GLM']
    
    # Calculating Pearson correlation
    tmp = df_tmp.dropna(how='any', axis=0)
    if tmp.shape[0] > 3:
        corr = tmp.CCI.corr(tmp.GLM, method='pearson')
        
    # Combining datasets if data agrees
        if corr >= limit_pearson:
            tmp = df_tmp.groupby(pandas.Grouper(freq='MS')).mean()
            tmp = tmp.loc[dates,:]
            tmp = numpy.where(pandas.isna(tmp.CCI), tmp.GLM, tmp.CCI)
            loc = pixels.index(pixel)
            df[:, loc] = tmp
            
    # Assigning timeseries to 3D dataset
            lon, lat = eval(pixel)[0], eval(pixel)[1]
            lon, lat = x.index(lon), y.index(lat)
            sm[:, lat, lon] = tmp
    
    # Checking progress
            pxls.append(pixel)
    sys.stdout.write(f"\r  '-> evaluation progress: {round((pixels.index(pixel) + 1)/len(pixels)*100,3)}%   ")
    sys.stdout.flush()

df = pandas.DataFrame(df, index=dates, columns=pixels)
df.dropna(how='all', axis=1, inplace=True)
pixels = pxls.copy()
sys.stdout.write("\r  '-> evaluation complete\n  '-> finished")

# Processing .............................................................
print('\n> Processing started...')

# Identifying agricultural droughts
pool = multiprocess.Pool(processes = n_cores)
results = [pool.apply_async(extremes.droughts_identification, args=[df.iloc[:,loc].to_frame(), 'MS', 'MS', perc_vtlm, perc_cdpm, 0, None, None]) for loc in range(df.shape[1])]
output = ([p.get() for p in results])
droughts = numpy.empty([df.shape[0], df.shape[1]])
for loc in range(df.shape[1]):
    droughts[:, loc] = output[loc]
droughts = pandas.DataFrame(droughts, index=df.index, columns=df.columns)
print("  '-> agricultural droughts identified")

# Insert identified agricultural droughts timeseries to global-scale array
xa = numpy.full([df.shape[0], len(lats), len(lons)], numpy.nan)
for pixel in pixels:
    # Extracting timeseries of identified droughts
    ts = droughts.loc[:, pixel].to_numpy()
    lon, lat = eval(pixel)[0], eval(pixel)[1]
    lon, lat = x.index(lon), y.index(lat)
    xa[:, lat, lon] = ts
            
    # Printing progress
    sys.stdout.write(f"\r  '-> droughts allocation - progress: {round((pixels.index(pixel)+1)/len(pixels)*100,2)}%  ")
    sys.stdout.flush()
print(f"\n  '-> finished")

# Post-processing ........................................................
print('> Post-processing started...')
# Creating and exporting datasets
export_dataset(xa, 'droughts' , lats, lons, dates, perc_vtlm, path_out, t_frame, limit_pearson)
print("  '-> finished\nDone!")
