# Input ..................................................................
# SWU data
path_swu = './Huang_2018_v2/'  # Location where all Huang et al. (2018) sectoral water use datasets are stored
file_swu = None                # To speed up the process, choose one file to evaluate: {None, 1:cons_dom, 2:cons_elec, 3:cons_irr, 4:cons_liv, 5:cons_mfg, 6:withd_dom, 7:withd_elec, 8:withd_irr, 9:withd_liv, 10:withd_mfg}

# DHE datasets
path_dhe1 = './data/dhe/compounds_W5E5-GSIM-GRDC-CCI-GLEAM_30arcmin_monthly_1990-2010.nc'    # NetCDF file with compound events identified, including agricultural droughts
path_dhe2 = './data/dhe/compounds_W5E5-GSIM-GRDC_30arcmin_monthly_1990-2010.nc'              # NetCDF file with compound events identified, excluding agricultural droughts

# Parameters
period = [1990,2010]           # Time frame: [initial_year,final_year]

# Output
path_out = './data/analysis/'  # Output folder directory

# Packages ...............................................................
import xarray, pandas, datetime, numpy, os, sys, extremes, warnings
from dateutil.relativedelta import relativedelta
from scipy import stats
warnings.filterwarnings("ignore")

# Dictionaries ...........................................................
swu_dict = {'dom':'Domestic', 'elec':'Thermoelectric', 'irr':'Irrigation', 'liv':'Livestock', 'mfg':'Manufacture',
            'cons':'Consumption', 'withd':'Withdrawal'}
events = ['none','dhy','dag','dhy+dag','htw','dhy+htw','dag+htw','cmp']
seasons = {0:('djf',[12,1,2]), 1:('mam',[3,4,5]), 2:('jja',[6,7,8]), 3:('son',[9,10,11])}

# Functions ..............................................................
def import_swu(sector, path_swu, period):
    # Importing dataset and extracting properties
    path = os.path.join(path_swu,f'{sector}.nc')
    xa = xarray.open_dataset(path, decode_times=False)
    var = list(xa.keys())[0]
    lons, lats, times = xa.lon.values, xa.lat.values, xa.month.values
    
    # Extracting values and creating dataframe
    da = xa[var].values
    dates = pandas.date_range('01/01/1971', '31/12/2010', freq='MS')
    coords = list(zip(lons,lats))
    df = pandas.DataFrame(da, index=dates, columns=coords)
    
    # Time-framing dataset
    df = df[(df.index >= datetime.datetime(period[0],1,1)) & (df.index <= datetime.datetime(period[1],12,31))]
    df = df.T[df.T.max(axis=1) > 0].T   # To remove pixels where all values are zero
    return df

def detrending_series(ts):
    # Formatting
    dtrend = ts_swu.copy().dropna()
    dates = pandas.date_range(dtrend.index[0], dtrend.index[-1], freq='MS')
    dtrend = dtrend.reindex(dates, method='nearest').to_frame()
    
    # Linear fitting
    annual = dtrend.groupby(pandas.Grouper(freq='AS')).mean()
    annual.reset_index(inplace=True)
    x, y = annual.iloc[:,0].apply(lambda x: x.toordinal()), annual.iloc[:,1]
    m, b = numpy.polyfit(x, y, 1)
    
    # Detrending
    dtrend['dates'] = dtrend.index
    dtrend.dates = dtrend.dates.apply(lambda x: x.toordinal())
    dtrend['detrended'] = dtrend[sector] - (dtrend.dates * m + b)
    resid = dtrend.loc[:,'detrended'].rename('detrend')
    return resid

# Pre-processing .........................................................
print('> Pre-processing started')

## Selecting SWU datasets
swu_files = os.listdir(path_swu)
swu_files = [file[:-3] for file in swu_files if '.nc' in file and '_min' not in file]
if file_swu != None:
    swu_files = swu_files[file_swu-1:file_swu]

## Opening DHE relate datasets
_, _, df_dhe1, _, _, _, _ = extremes.import_dataset(path_dhe1, 'compounds', period, None, None)
_, _, df_dhe2, _, _, _, _ = extremes.import_dataset(path_dhe2, 'compounds', period, None, None)

print("  '-> finished")

# Processing .............................................................
print('> Processing started')

## [SWU] Sectoral Water Use
# Iterating per sector/dimension
for sector in swu_files:
    dim, sec = swu_dict[sector.split('_')[0]], swu_dict[sector.split('_')[1]]
    df_out = pandas.DataFrame()
    
# Selecting suitable DHE dataset for evaluation
    if sector.split('_')[1] == 'irr': df_dhe = df_dhe1
    else: df_dhe = df_dhe2
    
# Opening dataset and selecting possible coordinates to evaluate
    df = import_swu(sector, path_swu, period)
    coords = [eval(str(x)) for x in df.columns]
    
# Extracting time series from dataset pixel
    for pixel in coords:
        lon, lat = pixel
        ts_swu = df.loc[:,[pixel]].squeeze().rename(sector)
        df_pxl = pandas.DataFrame()
        
# Decomposing time series
        detrend = detrending_series(ts_swu)
        df_pixel = pandas.DataFrame([ts_swu, detrend]).T
        
## [DHE] Drought-Heatwave Events
# Fetching information
        dhe = df_dhe.loc[:,str(pixel)].rename('dhe')
        df_pixel = pandas.concat([df_pixel,dhe], axis=1).dropna(subset=['detrend','dhe'], how='any')
        if df_pixel.shape[0] > 0:           # Checking if pixel has both SWU and DHE information
            df_pixel['months'] = df_pixel.index.month
            
## Responses analysis
# Filtering dataframe per extreme event, occurrence and input dataset
            for i in range(len(events)):    # DHE time-series: Types of events {'none','dhy','dag','dhy+dag','htw','dhy+htw','dag+htw','cmp'}
                event = events[i]
                df_tmp = df_pixel[df_pixel.dhe == i]
                    
# Calculating rank in percentile
                if df_tmp.shape[0] > 0:     # Checking if pixel still has both SWU and DHE information after filtering by extreme event
                    df_tmp = df_tmp.loc[:,['months','detrend']].dropna(how='any', axis=0)
                    perc, perc50, name50 = [], [], []
                    
                    for row,col in df_tmp.iterrows():
                        value, month = col.detrend, row.month
                        tmp = df_pixel.detrend
                        p_tmp = stats.percentileofscore(tmp, value, kind='strict')
                        perc.append(p_tmp)
                    ts_prc = pandas.DataFrame([df_tmp.months.tolist(), perc], index=['months', f'{event}']).T
                    
# Calculating median percentile value per event, sector, pixel and season
                    p50 = ts_prc[f'{event}'].quantile(0.5, interpolation='nearest')
                    perc50.append(p50), name50.append(f'{event}')
                    
                    for j in seasons.keys():
                        season, months = seasons[j]
                        ts_sea = ts_prc[(ts_prc.months == months[0]) | (ts_prc.months == months[1]) | (ts_prc.months == months[2])]
                        
                        if ts_sea.shape[0] > 0:  # Checking if pixel still has SWU and DHE inforamtion after filtering by extreme event and season
                            p50 = ts_sea[f'{event}'].quantile(0.5, interpolation='nearest')
                        else:
                            p50 = numpy.nan
                        perc50.append(p50), name50.append(f'{event}-{season}')
                else: perc50, name50 = [numpy.nan] * 5, [event] + [f'{event}-{seasons[s][0]}' for s in seasons.keys()]
                
# Compiling results
                df_perc = pandas.DataFrame(perc50, index=name50, columns=[(lon,lat)]).T
                df_pxl = pandas.concat([df_pxl, df_perc], axis=1)
            df_out = pandas.concat([df_out, df_pxl], axis=0)
            
# Printing progress
        sys.stdout.write(f"\r  '-> {dim} - {sec}: {round((coords.index(pixel)+1)/len(coords)*100,3):.3f}%")
        sys.stdout.flush()
    print('')
    
## Exporting
    df_out[['lon','lat']] = df_out.index.tolist()
    df_out = df_out.loc[:,['lon','lat'] + df_out.columns[:-2].tolist()].reset_index(drop=True)
    df_out.to_csv(os.path.join(path_out,f'global_{sector}_{period[0]}-{period[1]}.csv'), sep=';', index=False)
    
print('Done!')