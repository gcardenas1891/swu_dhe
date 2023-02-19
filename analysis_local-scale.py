# Input ..................................................................
# Datasets
path_swu = './data/swu/datasets_local-scale.xlsx'                                           # Location of local-scale summary dataset
path_htw = './data/dhe/heatwaves_W5E5_30arcmin_monthly_1990-2019_p=90%.nc'                  # Location of identified heatwaves
path_drt = './data/dhe/droughts_GSIM-GRDC_30arcmin_monthly_1990-2019_p=20%_gaps=25%.nc'     # Location of identified hydrological droughts
path_out = './data/analysis/'                                                               # Location of outputs

# Packages ...............................................................
import xarray, pandas, numpy, os
from scipy import stats

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
    dtrend['detrended'] = dtrend['swu'] - (dtrend.dates * m + b)
    resid = dtrend.loc[:,'detrended'].rename('detrend')
    
    return resid

# Pre-processing .........................................................
# Opening SWU dataset
df_coor = pandas.read_excel(path_swu, sheet_name='coords', index_col='Location')
df_data = pandas.read_excel(path_swu, sheet_name='swu', index_col=0, header=None).T

# Opening DHE datasets
ds_htw, ds_drt = xarray.open_dataset(path_htw), xarray.open_dataset(path_drt)
var_htw, var_drt = list(ds_htw.keys())[0], list(ds_drt.keys())[0]

# Formatting datasets
df_data = df_data.sort_values(by=['Sector','Location','City']).set_index(['Sector','Location','City'])
dates = pandas.to_datetime(xarray.decode_cf(ds_htw).time.values)
dhe_dict = {'htw':(ds_htw, var_htw), 'dhy':(ds_drt, var_drt)}

# Processing .............................................................
df_out = pandas.DataFrame()

# Evaluating per city
for row,col in df_data.iterrows():
    sector, loc, city = row
        
    # Extracting coordinates for location (lon,lat)
    ts_swu = df_data.loc[row,:].squeeze().dropna().rename('swu')
    ts_swu = ts_swu.reindex(pandas.date_range(ts_swu.index[0], ts_swu.index[-1], freq='MS'))
    detrend = detrending_series(ts_swu)
    df_pixel = pandas.DataFrame([ts_swu, detrend]).T
        
    # Fetching information from datasets and formatting
    for event in ['htw','dhy']:
        lon, lat = df_coor.loc[city,f'Lon_{event}'], df_coor.loc[city,f'Lat_{event}']
        ts_dhe = pandas.Series(dhe_dict[event][0].sel(lon=lon, lat=lat, method='nearest')[dhe_dict[event][1]].values, index=dates, name=event)
        ts_dhe = ts_dhe / ts_dhe.max()

        df_pixel = pandas.concat([df_pixel, ts_dhe], axis=1)
    df_pixel = df_pixel.loc[ts_swu.index,:].dropna(subset='htw')   # Time-framing to common time SWU and DHE
        
    # Identifying compound events and dropping duplicates
    df_pixel['dhy+htw'] = df_pixel.htw + df_pixel.dhy
    
    df_pixel['dhy+htw'] = numpy.where(df_pixel['dhy+htw'] == 2, 1, 0)
    df_pixel['htw']     = numpy.where(df_pixel['dhy+htw'] == 1, 0, df_pixel['htw'])
    df_pixel['dhy']     = numpy.where(df_pixel['dhy+htw'] == 1, 0, df_pixel['dhy'])
    df_pixel['none']    = numpy.where((df_pixel['htw'] != 1) & (df_pixel['dhy'] != 1) & (df_pixel['dhy+htw'] != 1), 1, 0)
    
    # Analysing percentiles
    df_pixel['months'] = df_pixel.index.month
    corr, head = list(row), ['Sector','Location','City']
    for event in ['none','dhy','htw','dhy+htw']:
        df_tmp = df_pixel[df_pixel[event] == 1]
        df_tmp = df_tmp.loc[:,['months','detrend']].dropna(how='any', axis=0)
        perc = []

        for y,x in df_tmp.iterrows():                           
            value, month = x['detrend'], y.month       
            tmp = df_pixel.detrend
            p_tmp = stats.percentileofscore(tmp, value, kind='strict')
            perc.append(p_tmp)

        if len(perc) != 0: p = numpy.median(perc)/100
        else: p = numpy.nan
        corr.append(p), head.append(f'{event}')
    df_per = pandas.DataFrame(corr, index=head).T
    df_out = pandas.concat([df_out, df_per], axis=0)

# Exporting results
df_out.to_csv(os.path.join(path_out,f'local_withd_dom-mfg.csv'), sep=';', index=False)
print('Done!')