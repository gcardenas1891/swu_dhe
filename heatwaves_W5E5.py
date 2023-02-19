# Input ..................................................................
# Directories
path_temp = './inputs/tasmax_W5E5v2.0_19790101-20191231.nc'
t_frame = [1990,2010]             # Time-frame to evaluate: [initial_year, final_year]
t_scale = 'MS'                    # Time aggregation: [None, 'MS', 'AS']. Where: {None=daily, MS=monthly, AS=yearly}

# Heatwaves analysis properties
days_cons = 3                     # Number of consecutive days to be considered a heatwaves
days_window = 5                   # If: t_scale == None, number of days to evaluate the threshold value. As the main day will be placed in the middle of the time window, the values must be odd numbers
perc_vtlm = 0.90                  # Percentile that defines the threshold value to identify a heatwave for the Variable Threshold Level method
n_cores = 50                      # Number of cores to be used during the multiprocessing

# Output
path_out = './data/dhe/'

# Packages ...............................................................
import extremes, pandas, xarray, numpy, multiprocess
pandas.options.mode.chained_assignment = None

# Pre-processing .........................................................
print('  Pre-processing started...')
# Importing NetCDF and defining raster characteristics
xa = xarray.open_dataset(path_temp)
var = list(xa.keys())[0]

# Defining space related properties
lons, lats = xa['lon'].values, xa['lat'].values
coords = [str((lon,lat)) for lat in lats for lon in lons]

# Defining time related properties
times = xa.time.values
to, tf = numpy.datetime64(f'{t_frame[0]}-01-01'), numpy.datetime64(f'{t_frame[1]}-12-31')
to, tf = times[numpy.where(times==to)][0], times[numpy.where(times==tf)][0]

# Converting xarray information to 2D array and data frame
xa = xa.sel(time=slice(to, tf))
da = xa[var].values
ar = da.reshape(da.shape[0], da.shape[1] * da.shape[2])

dates = pandas.to_datetime(xa.time.values)
df = pandas.DataFrame(ar, index=dates, columns=coords)

# Aggregating results by time
if t_scale != None: dates = pandas.date_range(to, tf, freq=t_scale)
print("   '-> finished")

# Processing .............................................................
print('  Processing started...')

pool = multiprocess.Pool(processes = n_cores)
results = results = [pool.apply_async(extremes.heatwaves_identification, args=[df.iloc[:,loc].to_frame(), perc_vtlm, days_cons, days_window, t_scale, None]) for loc in range(df.shape[1])]
output = ([p.get() for p in results])
rows = output[0].shape[0]
heatwaves = numpy.empty([rows, df.shape[1]])
for loc in range(df.shape[1]):
    heatwaves[:, loc] = output[loc]
heatwaves = heatwaves.reshape(rows, da.shape[1], da.shape[2]).astype('float32')

print(f"  '-> finished")

# Post-processing ........................................................
print('  Post-processing started')

# Defining time properties
year, month, day = dates[0].year, dates[0].month, dates[0].day
date_ref = dates[0].toordinal()
times = [t.toordinal()-date_ref for t in dates]

# Setting-up dataset
ar = xarray.DataArray(
    data = heatwaves,
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
ds = xarray.Dataset({'heatwaves': ar})
ds.attrs['Developer'] = 'Gabriel A. Cardenas B.'
ds.attrs['Institution'] = 'Department of Physical Geography, Faculty of Geosciences, Utrecht University'

# Exporting dataset
time_dict = {'D':'daily', 'MS':'monthly', 'AS':'yearly'}
if t_scale == None: t_scale = 'D'
perc = str(int(perc_vtlm*100)).zfill(2)
encode = {'heatwaves': {'dtype':'float32', 'zlib':True, 'complevel':9}}
ds.to_netcdf(f'{path_out}\\heatwaves_W5E5_30arcmin_{time_dict[t_scale]}_{t_frame[0]}-{t_frame[1]}_p={perc}%.nc', encoding=encode)
print("   '-> finished\nDone!")
