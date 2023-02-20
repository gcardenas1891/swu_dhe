# Input ..................................................................
# Directories
path_dis  = './data/inputs/gsim-grdc_discharges.csv'         # Discharge dataset: GRDC and GSIM combined
path_gsim = './data/inputs/GSIM_metadata/'                   # Complete GSIM dataset directory

# Time-series and spatial parameters
t_frame = [1990,2010]     # Period of analysis: [starting_year. ending_year]
t_scale = 'MS'            # Temporal scale: {'D':daily, 'MS':monthly, 'AS':annually}
t_gaps  = 0.25            # Percentage of temporal gaps that each discharge time series can have
res = 30                  # Spatial resolution of output dataset in arcmin

# Droughts identification parameters
perc_vtlm = 0.2           # Percentile that defines the threshold value to identify a drought in dates with runoff for the Variable Threshold Level method
perc_cdpm = 0.8           # Percentile that defines the threshold value to identify a drought in dates without runoff for the Consecutive Dry Period method
n_cores = 50              # Number of cores to be used during the multiprocessing

# Output
path_out = './data/dhe/'

# Packages ...............................................................
import pandas, numpy, datetime, os, sys, xarray, multiprocess, extremes, geopandas

# Functions ..............................................................
# Resampling dataframe to daily scale
def resampler(df, scale, t_scale):
    tmp = df.shift(1, freq=scale).iloc[-1,:].to_frame().T
    df = pandas.concat([df,tmp], axis=0)
    df = df.resample('D').ffill()
    df.drop(tmp.index, inplace=True)
    return df

# Exporting NetCDF
def export_dataset(xa, var, lats, lons, dates, perc, path_out, t_frame, t_scale, t_gaps):
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
    time_dict = {'D':'daily', 'MS':'monthly', 'AS':'yearly'}
    res = int(abs(lons[0]-lons[1]) * 60)
    perc = int(perc*100)
    encode = {var: {'dtype':'float32', 'zlib':True, 'complevel':9}}
    ds.to_netcdf(f'{path_out}/{var}_GSIM-GRDC_{res}arcmin_{time_dict[t_scale]}_{t_frame[0]}-{t_frame[1]}_p={perc}%_gaps={int(t_gaps*100)}%.nc', encoding=encode)

# Pre-processing .........................................................
print('  Pre-processing started...')
# Opening discharge dataframe
df = pandas.read_csv(path_dis, sep=';', index_col=0).T
df.index = pandas.to_datetime(df.index)

# Opening index discharge dataframes and formatting
idx1 = pandas.read_csv(f'{path_gsim}/GSIM_catalog/GSIM_metadata.csv', sep=',', index_col='gsim.no')
idx2 = pandas.read_csv(f'{path_gsim}/GSIM_catalog/GSIM_catchment_characteristics.csv', sep=',', index_col='gsim.no')
idx = pandas.concat([idx1, idx2], axis=1)

# Creating base grid
dx = dy = res/60
lons, lats = numpy.arange(-180+dx/2, 180, dx).tolist(), numpy.arange(90-dy/2, -90, -dy).tolist()
coords = [(lon, lat) for lon in lons for lat in lats]
data = pandas.DataFrame(coords, columns=['lon','lat'])
grd = geopandas.GeoDataFrame(data, geometry=geopandas.points_from_xy(data.lon, data.lat), crs='EPSG:4326')
print("  '-> finished")

# Processing .............................................................
print('  Processing started...')
# Time framing discharge data frame
df_dis = df[(df.index >= datetime.datetime(t_frame[0],1,1)) & (df.index <= datetime.datetime(t_frame[1],12,31))]

# Filtering station with allowed number of data-gaps
stations = [station for station in df_dis.columns if df_dis[station].count() >= (int((1-t_gaps) * df_dis.shape[0]))]
df_dis = df_dis.loc[:,stations]

# Identifying hydrological droughts
scale = 'MS' if t_scale == 'D' else t_scale
pool = multiprocess.Pool(processes = n_cores)
results = [pool.apply_async(extremes.droughts_identification, args=[df_dis.iloc[:,loc].to_frame(), 'MS', scale, perc_vtlm, perc_cdpm, 0, None, None]) for loc in range(df_dis.shape[1])]
output = ([p.get() for p in results])
droughts = numpy.empty([df_dis.shape[0], df_dis.shape[1]])
for loc in range(df_dis.shape[1]):
    droughts[:, loc] = output[loc]
droughts = pandas.DataFrame(droughts, index=df_dis.index, columns=df_dis.columns)
if t_scale == 'D':
    df_dis   = resampler(df_dis, scale, t_scale)
    droughts = resampler(droughts, scale, t_scale)
print("  '-> hydrological droughts identified")

# Filtering index discharge data frame
df_idx = idx.loc[df_dis.columns,:].reset_index(drop=False).rename(columns={'index':'gsim_id'})
df_idx.dropna(subset='area.est', axis=0, inplace=True)
df_idx.sort_values(by='area.est', ascending=False, inplace=True)

# Creating dictionary of stations
df_idx.index = range(1, df_idx.shape[0]+1)
stations = df_idx.gsim_id.tolist()
idx_dict, dis_dict = dict(zip(df_idx.index, stations)), dict(zip(stations, df_idx.index))

# Defining influence area of each station: catchments
da = numpy.full([len(lats), len(lons)], numpy.nan)
for station in stations:
    filename = f'{path_gsim}/GSIM_catchments/{station.lower()}.shp'
    
    # Importing basin shapefile and framing base grid extent
    if os.path.exists(filename):
        bsn = geopandas.read_file(filename).set_index('Id')
        bnd = bsn.bounds
        frm = grd[(grd.lon >= bnd.loc[0,'minx'] - dx) & (grd.lon <= bnd.loc[0,'maxx'] + dx) & (grd.lat >= bnd.loc[0,'miny'] - dy) & (grd.lat <= bnd.loc[0,'maxy'] + dy)]
        
    # Intersecting geodataframes and rasterizing results
        tmp = geopandas.overlay(frm, bsn, how='intersection')
        for pixel,coord in tmp.iterrows():
            lon, lat = lons.index(coord['lon']), lats.index(coord['lat'])
            da[lat, lon] = int(dis_dict[station])
    
    # Printing progress
    sys.stdout.write(f"\r  '-> catchment delineation - progress: {round((stations.index(station)+1)/len(stations)*100,2):.2f}%  ")
    sys.stdout.flush()
da = numpy.where(da > 0, da, numpy.nan)
sys.stdout.write(f"\r  '-> catchment delineation - progress: finished\n")

# Insert identified hydrological droughts timeseries to catchment array
dates = df_dis.index
xa = numpy.full([len(dates), len(lats), len(lons)], numpy.nan)
ds = xa.copy()
for ind in range(int(numpy.nanmin(da)), int(numpy.nanmax(da)+1)):
    # Extracting timeseries of identified droughts
    ts = droughts.loc[:, idx_dict[ind]].to_numpy()
    ro = df_dis.loc[:, idx_dict[ind]].to_numpy()
    
    # Locating pixel to insert data
    pixels = numpy.where(da == ind)
    for i in range(len(pixels[0])):
        lat, lon = pixels[0][i], pixels[1][i]
        xa[:, lat, lon] = ts
        ds[:, lat, lon] = ro
    
    # Printing progress
    sys.stdout.write(f"\r  '-> droughts allocation - progress: {round((ind)/(numpy.nanmax(da)+1)*100,2):.2f}%")
    sys.stdout.flush()
sys.stdout.write("\r  '-> droughts allocation - progress: finished\n  '-> finished")

# Post-processing ........................................................
print('\n  Post-processing started...')
export_dataset(xa, 'droughts' , lats, lons, dates, perc_vtlm, path_out, t_frame, t_scale, t_gaps)
print("  '-> finished\nDone!")
