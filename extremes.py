# Packages ....................................................................
import xarray, pandas, numpy, time, sys, os, datetime, multiprocess

# Functions ...................................................................
## Defining time step ...............................
def time_step(dates):
    dt = (dates[1] - dates[0]).days
    if   dt == 1: t_step = 'D'
    elif (dt > 1) and (dt <= 31): t_step = 'MS'
    elif (dt > 31) and (dt <= 366): t_step = 'AS'
    return t_step

## Extracting data set properties ...................
def dataset_properties(ds_path):
    # Importing NetCDF and defining raster characteristics
    xa = xarray.open_dataset(ds_path, decode_times=False)
    
    # Defining space related properties
    try:
        lons, lats = xa['lon'].values, xa['lat'].values
    except:
        try:
            lons, lats = xa['longitude'].values, xa['latitude'].values
        except:
            lons, lats = xa['X'].values, xa['Y'].values
    nrow, ncol = len(lats), len(lons)
    
    # Defining time related properties
    times = xa['time'].values
    dates = pandas.to_datetime(xarray.decode_cf(xa).time.values)
    t_step = time_step(dates)
    
    return xa, lons, lats, times, dates, t_step

## Framing datetime .................................
def date_frame(times, dates, t_frame, t_step):
    # Indexing dates in a dataframe
    df_dates = pandas.DataFrame([times,dates], index=['xr','dt']).T
    
    # Framing to period of evaluation
    date_o, date_f = datetime.datetime(t_frame[0],1,1), datetime.datetime(t_frame[1],12,31)
    df_dates = df_dates[(df_dates.dt >= date_o) & (df_dates.dt <= date_f)]
    df_dates.reset_index(drop=True, inplace=True)
    
    # Formatting useful dates
    to, tf = df_dates.xr[df_dates.index[0]], df_dates.xr[df_dates.index[-1]]
    dates = pandas.date_range(date_o, date_f, freq=t_step)   # to correct dates format
    
    return dates, date_o, date_f, to, tf

## Importing NetCDF .................................
def import_dataset(ds_path, ds_variable, t_frame, t_scale, eval_type):
    # Importing NetCDF and defining raster characteristics
    xa, lons, lats, times, dates, t_step = dataset_properties(ds_path)
    
    # Framing period of evaluation
    dates, date_o, date_f, to, tf = date_frame(times, dates, t_frame, t_step)

    # Setting list of pixel locations (tuple of coordinates) for model dataset
    coords = [str((lon,lat)) for lat in lats for lon in lons]
    
    # Converting xarray information to 2D array and data frame
    xa = xa.sel(time=slice(to, tf))
    da = xa[ds_variable].values
    ar = da.reshape(da.shape[0], da.shape[1] * da.shape[2])
    df = pandas.DataFrame(ar, index=dates, columns=coords)
    
    # Aggregating results by time
    if t_scale != None:
        dates = pandas.date_range(date_o, date_f, freq=t_scale)
        if eval_type == 'drought':
            df = df.groupby(pandas.Grouper(freq=t_scale)).mean()
    
    return xa, da, df, dates, t_step, lons, lats

## Time window time series ..........................
def time_window(ts, days_window):
    dp = ts.copy()
    window = [int((days_window-1)/2) - dt for dt in range(days_window)]
    window.remove(0)
    for t in window:
        dp = pandas.concat([dp, ts.shift(t)], axis=1)
    return dp

## Variable Threshold Level method ..................
def vtl_method(ts, t_step, t_scale, days_window=5):
    # Formatting data frame
    station = ts.columns[0]
    ts_vtl = ts.copy()
    ts_vtl['year'] = ts_vtl.index.year
    ts_vtl['month_day'] = ts_vtl.index.strftime('%m-%d')
    ts_vtl = ts_vtl.pivot('month_day', 'year', station).reset_index().rename_axis(None, axis=1)

    # Calculating percentiles for daily time step
    if (t_step == 'D') and (t_scale == None):
        # Identifying location of 29/feb
        leap_day = datetime.datetime(2000,2,29).strftime('%m-%d')
        leap_day_row = ts_vtl[ts_vtl.month_day == leap_day].index[0]
        ts_vtl.drop('month_day', axis=1, inplace=True)
    
        # Identify the 'fake' leap years: make them "-1"
        leap_years = [year for year in ts_vtl.columns if ((year%400==0) or (year%100!=0)) and (year%4==0)]
        not_leap_year_columns = [year for year in ts_vtl.columns if year not in leap_years]
    
        # Generating the time window (dp: data point)
        dp = time_window(ts_vtl, days_window)
        dp_prctl = dp.rank(pct=True, axis=1, na_option='keep')
        prctl = dp_prctl.iloc[:, :len(ts_vtl.columns)]
    
        # Reshaping data frame to series
        prctl.loc[leap_day_row, not_leap_year_columns] = -1
        prctl = prctl.T.stack(dropna=False)
        prctl = prctl[prctl != -1]  # Remove 'fake' leap years
    
    # Calculating percentiles for monthly and yearly time steps
    else:
        ts_vtl.drop('month_day', axis=1, inplace=True)
        prctl = ts_vtl.rank(pct=True, axis=1, na_option='keep')
        prctl = prctl.T.stack(dropna=False)
    
    # Turning percentile results Numpy array type
    prctl = prctl.reset_index(drop=True)
    return prctl

## Consecutive Dry Period method ....................
def cdp_method(ts, perc_cdpm):
    # Getting number of consecutive days without discharge (=zero) ('nd' stands for 'no discharge')
    ts_cdp = ts.reset_index(drop=True).squeeze()
    ts_nd = numpy.where(ts_cdp == 0, 1, 0)
    ts_nd = pandas.Series(ts_nd)
    ts_nd_cum = ts_nd * (ts_nd.groupby((ts_nd != ts_nd.shift()).cumsum()).cumcount() + 1)
    
    # Calculating threshold value CDPM
    threshold = ts_nd_cum.quantile(perc_cdpm)
    prctl_cdpm = ts_nd_cum.rank(pct=True, axis=0)
    return ts_cdp, ts_nd, threshold

# Combined approach .................................
def comb_method(ts, t_step, t_scale, perc_vtlm, perc_cdpm, days_window):
    # Calling results from other methods
    prctl_vtl = vtl_method(ts, t_step, t_scale, days_window)
    ts_cdp, ts_nd, threshold = cdp_method(ts, perc_cdpm)
    
    # Getting number of drought days from VTLM
    prctl_vtl_pos = prctl_vtl.reset_index(drop=True)
    prctl_vtl_pos[ts_cdp == 0] = numpy.nan
    prctl_vtl_pos = pandas.Series(prctl_vtl_pos)
    ts_vtl_pos = numpy.where(prctl_vtl_pos <= perc_vtlm, 1, 0)
    ts_vtl_pos = pandas.Series(ts_vtl_pos)
    
    # Generating time series with either a [no-flow day] or a [tlvm drought]
    ts_combined = numpy.where((ts_nd == 1) | (ts_vtl_pos == 1), 1, 0)
    ts_combined = pandas.Series(ts_combined)
    ts_combined_cum = ts_combined * (ts_combined.groupby((ts_combined != ts_combined.shift()).cumsum()).cumcount() + 1)
    
    # Calculating percentiles for combined time series, rescaling and compiling together with VTLM percentiles
    prctl_comb = ts_combined_cum.rank(pct=True, axis=0, na_option='keep')
    prctl = prctl_vtl_pos.where(prctl_vtl_pos.notna(), 1-prctl_comb)
    return prctl

## Extremes severity classification .................
def severity_classification(severity, prctl):
    # Formatting input data
    severity.sort()
    extreme = prctl.copy()
    
    # Assigning severity value
    for i in range(len(severity)):
        perc = severity[i] / 100
        sev = len(severity) - i
        extreme[extreme < perc] = sev + 1      # this is done to prevent that percentiles == 1 are classified as a extreme event of severity == 1
    extreme[extreme <= 1] = 1
    extreme = extreme - 1                      # rescaling severity values: 0 for no extreme events identified
    return extreme

## Droughts identification ..........................
def droughts_identification(ts, t_step, t_scale, perc_vtlm, perc_cdpm, q_zero, days_window, severity):
    # Identifying pixels without NaN data
    if ts.isna().sum()[0] < ts.shape[0]:

    # Converting values lower than 'q_zero' equal to zero
        ts[ts <= q_zero] = 0
    
    # Droughts identification method selection
        if   ts.squeeze().quantile(0.05) > 0 : prctl = vtl_method( ts, t_step, t_scale, days_window)
        elif ts.squeeze().quantile(0.05) == 0: prctl = comb_method(ts, t_step, t_scale, perc_vtlm, perc_cdpm, days_window)
    
    # Assigning severity
        droughts = prctl.copy().to_numpy()
        if severity == None:
            droughts[droughts > perc_vtlm] = 0
            droughts[(droughts <= perc_vtlm) & (droughts > 0)] = 1
        elif type(severity) == list:
            droughts = severity_classification(severity, droughts)
    
    else: droughts = ts.squeeze().to_numpy()
    return droughts

## Consecutive days analysis ........................
def days_consecutive(prctl, perc_heat, days_cons):
    ts_heat = numpy.where(prctl >= perc_heat, 1, 0)
    ts_heat = pandas.Series(ts_heat)
    ts_heat_cum = ts_heat * (ts_heat.groupby((ts_heat != ts_heat.shift()).cumsum()).cumcount() + 1)
    heat_idx = ts_heat_cum[ts_heat_cum >= days_cons].index
    noheat_idx = ts_heat_cum[ts_heat_cum < days_cons].index
    return heat_idx, noheat_idx    

## Heatwaves identification independently ...........
def heatwaves_identification(ts, perc, days_cons, days_window, t_scale, severity):
    # Identifying only-NaN pixels
    if ts.isna().sum()[0] == 0:
    
    # Filtering to periods with high temperature;
        hot_limit = ts.squeeze().quantile(0.5)
        ts = ts.where(ts >= hot_limit, numpy.nan)
    
    # Fetching results from Variable Threshold Level method
        prctl = vtl_method(ts, 'D', None, days_window)
    
    # Identifying heatwaves according the number of consecutive days
        heat_idx, noheat_idx = days_consecutive(prctl, perc, days_cons)
    
    # Assigning severity
        heatwaves = 1 - prctl                          # to reverse the percentiles and to make them comparable with droughts format
        if severity == None:
            heatwaves.loc[heat_idx] = 1
        else:
            severity = [100-sev for sev in severity]   # to reverse the severity and to make it comparable with droughts format
            heatwaves = severity_classification(severity, heatwaves)
        heatwaves.loc[noheat_idx] = 0
    
    # Aggregating results by time
        if t_scale != None:
            heat = numpy.where(heatwaves > 0, 1, 0)
            ts_heat = pandas.Series(heat, index=ts.index)
            ts_heat = ts_heat.groupby(pandas.Grouper(freq=t_scale)).sum()
            ts_heat = numpy.where(ts_heat >= 3, 1, 0)   # '3' is the minimum number of days per month with identified heatwave to aggregate results
            
            ts_sev  = pandas.Series(heatwaves.to_numpy(), index=ts.index)
            ts_sev  = ts_sev.groupby(pandas.Grouper(freq=t_scale)).max()
            ts_sev  = ts_sev.to_numpy()
            
            heatwaves = ts_heat * ts_sev
    
    else: heatwaves = ts.squeeze().to_numpy()
    return heatwaves

## Compound events identification ....................
def compound_identification(ds1, ds2, ds3, grouping):
    # Opening datasets
    ds1 = xarray.open_dataset(ds1, decode_times=False)
    da1 = ds1['heatwaves'].values.astype(int)
    try:
        ds2 = xarray.open_dataset(ds2, decode_times=False)
        da2 = ds2['droughts'].values
    except: da2 = da1 * 0
    try:
        ds3 = xarray.open_dataset(ds3, decode_times=False)
        da3 = ds3['droughts'].values
    except: da3 = da1 * 0
    
    # Converting extreme event values to compound ID (Hydrological droughts: 100000, Agricultural drought: 200000, Heatwave: 400000)
    heat, drt1, drt2 = da1, da2, da3
    for sev in range(1,6):
        drt1 = numpy.where(drt1 == sev, int(1000000 + (sev * 10**(sev))), drt1)
        drt2 = numpy.where(drt2 == sev, int(2000000 + (sev * 10**(sev))), drt2)
        heat = numpy.where(heat == sev, int(4000000 + (sev * 10**(sev))), heat)
    
    # Combining the datasets 
    comp_all = heat + drt1 + drt2
    
    # Creating data array with extreme values only (not severity included) [Drought: 1, Heatwave: 2, Compound: 3]
    with numpy.errstate(invalid='print'):
        compounds = comp_all // 1000000
    
    if grouping == True:
        compounds = numpy.where((compounds == 1) | (compounds == 2) | (compounds == 3), 1, compounds)
        compounds = numpy.where((compounds == 4) , 2, compounds)
        compounds = numpy.where((compounds == 5) | (compounds == 6) | (compounds == 7), 3, compounds)
    
    # Creating data array with severity values only
    with numpy.errstate(invalid='print'):
        severity = comp_all % 1000000
    for sev in range(1,6):
        sev = 6 - sev
        severity = numpy.where(severity >= (sev * 10**(sev)), sev, severity)
    return compounds, severity

## Multiprocessing ...................................
def processing_timeseries(eval_type, df, da, t_step, t_scale, perc=0.2, perc_cdpm=0.8, q_zero=0, days_window=5, days_cons=3, severity=[20], n_cores=4):
    pool = multiprocess.Pool(processes = n_cores)
    if eval_type == 'drought':
        results = [pool.apply_async(droughts_identification,  args=[df.iloc[:,loc].to_frame(), t_step, t_scale, perc, perc_cdpm, q_zero, days_window, severity]) for loc in range(df.shape[1])]
    elif eval_type == 'heatwave':
        results = [pool.apply_async(heatwaves_identification, args=[df.iloc[:,loc].to_frame(), perc, days_cons, days_window, t_scale, severity]) for loc in range(df.shape[1])]
    output = ([p.get() for p in results])
    rows = output[0].shape[0]
    extremes_all = numpy.empty([rows, df.shape[1]])
    for loc in range(df.shape[1]):
        extremes_all[:, loc] = output[loc]
    extremes = extremes_all.reshape(rows, da.shape[1], da.shape[2]).astype('float32')
    return extremes

## Creating datasets .................................
def create_dataset(ds_output, lons, lats, dates, var_name):
    # Defining time properties
    year, month, day = dates[0].year, dates[0].month, dates[0].day
    date_ref = dates[0].toordinal()
    times = [t.toordinal()-date_ref for t in dates]
    
    # Setting-up dataset
    ar = xarray.DataArray(
        data = ds_output,
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
    ds = xarray.Dataset({var_name: ar})
    ds.attrs['Institution'] = 'Department of Physical Geography, Faculty of Geosciences, Utrecht University'
    return ds
 
## Exporting results .................................
def netcdf_export(ds_output, var_name, lons, lats, dates, perc, out_path, out_file='extremes'):
    # Creating dataset
    xa_output = create_dataset(ds_output, lons, lats, dates, var_name)
    
    # Setting file name
    dt_dict = {'D':'daily', 'MS':'monthly', 'AS':'yearly'}
    step = time_step(pandas.to_datetime(xarray.decode_cf(xa_output).time.values))
    res = str(int(abs(lons[0]-lons[1])*60)).zfill(2)
    if perc == None: 
        out_file = f'{var_name}_{out_file}_{res}arcmin_{dt_dict[step]}_{dates[0].year}-{dates[-1].year}'
    else:
        out_file = f'{var_name}_{out_file}_{res}arcmin_{dt_dict[step]}_{dates[0].year}-{dates[-1].year}_p={str(int(perc*100)).zfill(2)}%'
    path_out = os.path.join(out_path, out_file)
    
    # Exporting dataset
    encode = {var_name: {'dtype':'float32', 'zlib':True, 'complevel':9}}
    xa_output.to_netcdf(path_out+'.nc', encoding=encode)