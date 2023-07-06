# Input ..................................................................
# SWU data
path_swu = './inputs/datasets_global-scale_aquastat.csv'   # File location of recorded water use data: datasets_global-scale_aquastat.csv, datasets_country-scale_usgs.csv
path_shp = './inputs/shp/aquastat_countries_30arcmin.shp'  # Shapefile location of regions: aquastat_countries_30arcmin.shp, usgs_counties_30arcmin.shp

# DHE datasets
path_dhyd = './dhe/droughts_GSIM-GRDC_30arcmin_monthly_1990-2019_p=20%_gaps=25%.nc'      # NetCDF file with hydrological droughts identified
path_dagr = './dhe/droughts_CCI-GLEAM_30arcmin_monthly_1990-2019_p=20%_pearson=40%.nc'   # NetCDF file with agricultural droughts identified
path_heat = './dhe/heatwaves_W5E5_30arcmin_monthly_1990-2019_p=90%.nc'                   # NetCDF file with heatwaves identified

# Parameters
threshold_drought  = [20,3]   # [Minimum % of country/county area affected by a drought , Minimum number of months in a year with a drought]
threshold_heatwave = [20,3]   # [Minimum % of country/county area affected by a heatwave, Minimum number of months in a year with a heatwave]

# Output
path_out  = '.data/analysis/'

# Packages ...............................................................
import pandas, geopandas, xarray, numpy, os, sys, warnings
from scipy import stats
warnings.filterwarnings('ignore')

# Functions ..............................................................
def open_dataset(path_dhe, df):
    # Opening dataset and extracting time series
    ds  = xarray.open_dataset(path_dhe)
    var = list(ds.keys())[0]
    da = ds[var].values
    ar = da.reshape(da.shape[0], da.shape[1] * da.shape[2])

    dates = pandas.to_datetime(ds.time)
    lons, lats = ds['lon'].values, ds['lat'].values
    coords = [(lon,lat) for lat in lats for lon in lons]

    pixels = list(df.pixel.unique())
    df_dhe = pandas.DataFrame(ar, index=dates, columns=coords)
    df_dhe = df_dhe.loc[:,pixels]
    return df_dhe

def formatting_dhe(dhe, threshold, months, label):
    # Filtering DHE's dataframe to annual scale
    dhe
    dhe.index = pandas.to_datetime(dhe.index)
    dhe.dropna(how='all', axis=1, inplace=True)
    dhe.loc[:,:] = numpy.where(dhe >= threshold/100, 1, 0)
    dhe = dhe.groupby(pandas.Grouper(freq='AS')).sum()
    dhe.loc[:,:] = numpy.where(dhe >= months, 1, 0)
    
    # Formatting datetime
    dhe['Year'] = dhe.index.year
    dhe = dhe.set_index('Year').T
    
    # Melting dataset
    out = pandas.DataFrame()
    for location in dhe.index:
        tmp = dhe.loc[location,:].rename(label).to_frame()
        tmp.reset_index(inplace=True)
        tmp['ID'] = location
        out = pandas.concat([out,tmp], axis=0)
    out.set_index(['ID','Year'], inplace=True)
    return out

# Pre-processing .........................................................
print('Pre-processing started...')

# Formatting 
gdf = geopandas.read_file(path_shp)
df  = pandas.DataFrame(gdf)
df['pixel'] = list(zip(df.x, df.y))

# Opening dataset and extracting time series
dhyd = open_dataset(path_dhyd, df)
dagr = open_dataset(path_dagr, df)
heat = open_dataset(path_heat, df)
    
# Assigning extreme event to location
df_dhyd, df_dagr, df_heat = pandas.DataFrame(), pandas.DataFrame(), pandas.DataFrame()
locations = df.ID.unique().tolist()

for loc in locations:
    df_tmp   = df[df.ID == loc]
    area_tot = df_tmp.area.sum()
    headers  = df_tmp.pixel.tolist()
    
    for event in ['dhyd','dagr','heat']:
        df_dhe = eval(event)
        dhe = df_dhe.loc[:,headers]

        for row,col in df_tmp.iterrows():
            pixel, area = col.pixel, col.area
            tmp = dhe[pixel]
            dhe[pixel] = numpy.where(tmp == 1, area/area_tot, 0)
    
        dhe_loc = dhe.sum(axis=1)
        globals()[f'df_{event}'] = pandas.concat([globals()[f'df_{event}'],dhe_loc.rename(str(loc))], axis=1)
    
    # Printing progress
    sys.stdout.write(f"\r  '-> locations evaluated: {locations.index(loc)+1}/{len(locations)}")
    sys.stdout.flush()

# Defining swu dataset properties
events  = ['Droughts','Heatwaves','Compounds']
source = os.path.basename(path_swu).split('_')[-1].split('.')[0]
if   source == 'aquastat':
    sectors = ['Irrigation','Municipal','Industrial','Thermoelectric','Livestock']
elif source == 'usgs':
    sectors = ['Irrigation','Domestic','Manufacturing','Thermoelectric','Livestock']
    
# Formatting DHE dataframe
hyd  = formatting_dhe(df_dhyd, threshold_drought[0] , threshold_drought[1] , 'Hydrological')
agr  = formatting_dhe(df_dagr, threshold_drought[0] , threshold_drought[1] , 'Agricultural')
heat = formatting_dhe(df_heat, threshold_heatwave[0], threshold_heatwave[1], 'Heatwaves')

# Opening SWU dataset and merging
df = pandas.read_csv(path_swu, sep=';', index_col=['ID','Year'], dtype={'ID':str})
df = pandas.concat([df,hyd,agr,heat], axis=1)

# Filtering merged dataset
df = df.dropna(subset=['Hydrological','Agricultural','Heatwaves'], how='all', axis=0)
df = df.dropna(subset=sectors, how='all', axis=0)

# Processing .............................................................
print('\nProcessing started...')
# Creating dictionary with evaluation areas and associated region
gdf = gdf[~pandas.isna(gdf.region)]
dict_region = dict(zip(gdf.ID.astype(str),gdf.region))

# Formatting dataframe
df[['region','drop']] = df.index.tolist()
df.drop('drop', axis=1, inplace=True)

# Associating locations to regions
df.region.replace(dict_region, inplace=True)
for row,col in df.iterrows():
    if col.region not in gdf.region.unique():
        df.drop(row, inplace=True)

# Saving excel file
scale  = 'global' if source == 'aquastat' else 'country'
period = os.path.basename(path_heat).split('_')[4]
filename = f'{path_out}{scale}_{source}_withd_{period}.xlsx'
writer = pandas.ExcelWriter(filename)

# Analysing responses by sectors
for sector in sectors:
    df_out = pandas.DataFrame()
    
    # Framing to useful columns
    drought = 'Agricultural' if sector == 'Irrigation' else 'Hydrological'
    df_swu = df.loc[:,[sector, drought, 'Heatwaves','region']]
    df_swu.rename(columns={drought:'Droughts'}, inplace=True)
    df_swu.dropna(subset=[sector], how='any', axis=0, inplace=True)
    
    # Filtering by MESSAGE region
    for region in df_swu.region.unique():
        df_tmp = df_swu[df_swu.region == region]
    
        # Calculating compound events
        df_tmp.loc[:,'Heatwaves'] = numpy.where(df_tmp.Heatwaves==1, 2, df_tmp.Heatwaves)
        df_tmp.loc[:,'DHE'] = df_tmp.Droughts + df_tmp.Heatwaves
        
        # Calculating the score percentile
        perc = []
        for row,col in df_tmp.iterrows():
            value = col[sector]
            p = stats.percentileofscore(df_tmp[sector], value, kind='strict')
            perc.append(p)
        df_tmp['Percentile'] = perc
        
        # Evaluating if there are extreme events
        if (df_tmp.DHE.max() > 0) and (df_tmp.shape[0] > 2):
            out = df_tmp.loc[:,['region',sector,'Percentile','DHE']]
            df_out = pandas.concat([df_out,out], axis=0)
            
    # Printing progress
    sys.stdout.write(f"\r  '-> sectors evaluated: {sectors.index(sector)+1}/{len(sectors)}")
    sys.stdout.flush()
    
    # Export spreadsheet
    df_out.to_excel(writer, sheet_name=sector.lower(), index=False)
writer.close()
print('\nDone!')
