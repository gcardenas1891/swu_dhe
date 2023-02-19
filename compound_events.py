# Input ..................................................................
# Directories
path_heat = './data/dhe/heatwaves_W5E5_30arcmin_monthly_1990-2010_p=90%.nc'                       # Netcdf with heatwaves identified
path_dhyd = './data/dhe/droughts_GSIM-GRDC_30arcmin_monthly_1990-2010_p=20%_gaps=25%.nc'          # Netcdf with hydrological droughts identified
path_dagr = './data/dhe/droughts_CCI-GLEAM_30arcmin_monthly_1990-2010_p=20%_pearson=40%.nc'       # Netcdf with agricultural droughts identified

# Output
path_out = './data/dhe/'

# Packages ...............................................................
import extremes, pandas, xarray, numpy, os
pandas.options.mode.chained_assignment = None

# Processing .............................................................
print('  Processing started...')

# Opening datasets
ds1 = xarray.open_dataset(path_heat)
da1 = ds1['heatwaves'].values.astype(int)
ds2 = xarray.open_dataset(path_dhyd)
da2 = ds2['droughts'].values
if path_dagr != None:
    ds3 = xarray.open_dataset(path_dagr)
    da3 = ds3['droughts'].values
    tag = '-CCI-GLEAM'
else:
    da3, tag = da1 * 0, ''

# Assigning compound ID (Hydrological droughts: 1, Agricultural drought: 2, Heatwave: 4) and combining datasets
heat, drt1, drt2 = da1 * 4, da2 * 1, da3 * 2
compounds = heat + drt1 + drt2

print(f"  '-> finished")

# Post-processing ........................................................
print('  Post-processing started')

# Defining spatial-temporal dataset properties
lons, lats = ds1['lon'].values, ds1['lat'].values
dates = pandas.to_datetime(ds1['time'].values)

# Defining output dataset time properties
year, month, day = dates[0].year, dates[0].month, dates[0].day
date_ref = dates[0].toordinal()
times = [t.toordinal()-date_ref for t in dates]

# Setting-up dataset
ar = xarray.DataArray(
    data = compounds,
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
ds = xarray.Dataset({'compounds': ar})
ds.attrs['Developer'] = 'Gabriel A. Cardenas B.'
ds.attrs['Institution'] = 'Department of Physical Geography, Faculty of Geosciences, Utrecht University'

# Exporting dataset
time_dict = {'D':'daily', 'MS':'monthly', 'AS':'yearly'}
t_scale = extremes.time_step(dates)
t_frame = [dates[0].year, dates[-1].year]

encode = {'compounds': {'dtype':'float32', 'zlib':True, 'complevel':9}}
ds.to_netcdf(os.path.join(path_out, f'compounds_W5E5-GSIM-GRDC{tag}_30arcmin_{time_dict[t_scale]}_{t_frame[0]}-{t_frame[1]}.nc'), encoding=encode)
print("   '-> finished\nDone!")
