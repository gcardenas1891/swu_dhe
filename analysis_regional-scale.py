# Input ..................................................................
path_swu = './data/inputs/datasets_regional-scale.xlsx'                              # Location of cooiling water or net electricity generation datasets from EIA
dataset  = 'Water'                                                                   # Type of dataset to evaluate: {'Water':cooling water, 'Electricity':net electricity production}
path_dhe = './data/dhe/compounds_W5E5-GSIM-GRDC_30arcmin_monthly_1990-2019.nc'       # Location of identified extreme events (hydrological droughts and heatwaves)
path_out = './data/analysis/'                                                        # Location of outputs

# Packages ...............................................................
import xarray, pandas, numpy, sys, warnings, os
from scipy import stats
warnings.filterwarnings("ignore")

# Functions ..............................................................
def detrending_series(ts):
    dtrend = ts.copy().to_frame()
    
    # Linear fitting
    annual = dtrend.groupby(pandas.Grouper(freq='AS')).mean()
    annual = annual.dropna().reset_index()
    x, y = annual.iloc[:,0].apply(lambda x: x.toordinal()), annual.iloc[:,1]
    m, b = numpy.polyfit(x, y, 1)
    
    # Detrending
    dtrend['dates'] = dtrend.index
    dtrend.dates = dtrend.dates.apply(lambda x: x.toordinal())
    dtrend['detrended'] = dtrend['ele'] - (dtrend.dates * m + b)
    resid = dtrend.loc[:,'detrended'].rename('detrend')
    return resid

# Pre-processing .........................................................
print('> Pre-processing started')
# Opening SWU data
df_coor = pandas.read_excel(path_swu, sheet_name='Coords', index_col='ID')
df_data = pandas.read_excel(path_swu, sheet_name=dataset, index_col=0)

# Opening DHE dataset
ds_dhe = xarray.open_dataset(path_dhe)
ds_dhe = ds_dhe.sel(lat=slice(50,20), lon=slice(-170,-50))
dates = pandas.to_datetime(xarray.decode_cf(ds_dhe).time.values)

# Filtering SWU results
df_data = df_data.T[df_data.max() > 0].T
df_data = df_data.where(df_data > 0, numpy.nan)
stations = df_data.columns.tolist()
print("  '-> finished")

# Processing .............................................................
print('> Processing started')
# Retrieving information
df_out = pandas.DataFrame()

# Extracting time series from dataset pixel
for station in stations:    
    # SWU data
    ts_swu = df_data.loc[:,station].squeeze().dropna().rename('ele')
    ts_swu = ts_swu.reindex(pandas.date_range(ts_swu.index[0], ts_swu.index[-1], freq='MS'))
    detrend = detrending_series(ts_swu)
    
    # DHE data
    lon, lat, fuel, csys = df_coor.loc[station, 'Lon'], df_coor.loc[station, 'Lat'], df_coor.loc[station, 'Fuel'], df_coor.loc[station, 'Sys']
    ts_dhe = pandas.Series(ds_dhe.sel(lon=lon, lat=lat, method='nearest').compounds.values, index=dates, name='dhe')
    
    # Combining data and formatting
    df_pixel = pandas.concat([ts_swu, detrend, ts_dhe], axis=1)
    df_pixel = df_pixel.loc[ts_swu.index,:]
    df_pixel.dropna(subset='dhe', inplace=True)
    
# Analysing dataframes
    # Creating a 'month' column
    df_pixel['months'] = df_pixel.index.month

    # Filtering dataframe per time series, extreme event, occurrance and input dataset
    events = ['none','dhy','dag','dhy+dag','htw','dhy+htw','dag+htw','cmp']
    corr, head = [lon,lat,fuel,csys], ['lon','lat','fuel','sys']
    
    for i in [0,1,4,5]:              # DHE time-series: Types of events {'none','dhy','htw','dhy+htw'}
        event = events[i]
        df_tmp = df_pixel[df_pixel['dhe'] == i]
                
    # Calculating rank in percentile
        df_tmp = df_tmp.loc[:,['months','detrend']].dropna(how='any', axis=0)
        perc = []
        
        for row,col in df_tmp.iterrows():
            value, month = col['detrend'], row.month
            tmp = df_pixel['detrend']
            p_tmp = stats.percentileofscore(tmp, value, kind='strict')
            perc.append(p_tmp)
            
    # Calculating median percentile value per event, sector and pixel
        if len(perc) != 0: p = numpy.median(perc)/100
        else: p = numpy.nan
        corr.append(p), head.append(f'{event}-detrend')
    
# Compiling results
    df_corr = pandas.DataFrame(corr, index=head).T
    df_out = pandas.concat([df_out, df_corr], axis=0)
    
# Printing progress
    sys.stdout.write(f"\r    '-> Progress: {stations.index(station)+1}/{len(stations)}")
    sys.stdout.flush()
print('')

# Exporting results
tag, dim = ('cooling','withd') if dataset == 'Water' else ('electricity','prod')
t1, t2 = str(max(df_data.index[0].year, dates[0].year)), str(min(df_data.index[-1].year, dates[-1].year))
df_out = df_out.set_index(keys=['lon','lat']).dropna(how='all', axis=0).reset_index(drop=False)
df_out.to_csv(os.path.join(path_out,f'regional_{dim}_elec-{tag}_{t1}-{t2}.csv'), sep=';', index=False)

print('Done!')
