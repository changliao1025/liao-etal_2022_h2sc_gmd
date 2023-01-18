######################################################
############# well i chosed ##########################
wellist = ['PP-01',
           'PP-02',#'PP-03',
           'PR-05',
           'PT-07','PR-07',
           'PR-10','PR-11']

wellat  = [-2.608976,
           -2.606521,#-2.605713,
           -2.604537,
           -2.604057,-2.602816,
           -2.601951,-2.601140]

wellon  = [-60.208901,
           -60.208150,#-60.207503,
           -60.207577,
           -60.207640,-60.207436,
           -60.207320,-60.207212]

###### well data observation ######
# 2.0) input the well data
year = '2004'
df_well = pd.read_csv('./data/1.0_well_obs_mon.csv', na_values='-9999', parse_dates=True, skiprows=0)
df_well.index = pd.to_datetime(df_well['date'])
#df_well = df_well[wellist]
# print(df_well)
well_obs_mon = df_well.loc[year]
well_obs_avg = well_obs_mon.mean(axis=0)
# print(well_obs_mon)
# print(well_obs_avg)

# 2.0) input the lat and lon
infil = '/qfs/people/lili400/compy/project/NGT/manaus/data/2d_data/k34_dem/k34_dem1_adj.nc'
ncfil =  xr.open_mfdataset(infil, parallel=True)
#print(ncfil)
ncdem = ncfil["Band1"][:]
dem0  = ncdem.values
lat   = ncfil.lat.values
lon   = ncfil.lon.values
#print(lat.shape)
print(ncdem.shape)
#print(ncdem.values)
#print(ncdem.sel(lat=wellat[0], lon=wellon[0], method='nearest').values)

# 2.2) find the ilat and ilon for well
welat_i = []
welon_i = []
for t in range(len(wellat)):
    welat_i.append( min(range(len(lat)), key=lambda i: abs(lat[i]-wellat[t])) )
    welon_i.append( min(range(len(lon)), key=lambda i: abs(lon[i]-wellon[t])) )
print('well ilat: ', welat_i,' , ',welon_i)
# print(ncdem[welat_i,welon_i])

#######################################################
#################### Parflow simulation ###############
expnm = 'par_35'
ymdat = '1905'

wtd_mod = np.load('./data/'+expnm+'_wtd.npy') # (nt, ny, nx)
olf_mod = np.load('./data/'+expnm+'_overland_flow.npy') # (nt, ny, nx)
olf_mod = olf_mod / 3600.    # m3/hr -> m3/s
print(wtd_mod.shape)  ### this is what you need ###

wsf_mod = np.load('./data/'+expnm+'_surface_storage.npy') # (nt, ny, nx)
wsb_mod = np.load('./data/'+expnm+'_subsurface_storage.npy') # (nt, ny, nx)
print(wsf_mod.shape)

# elm input
elm_pr = np.load('./data/'+expnm+'_elm_rain.npy') # (nt, ny, nx)
elm_et = np.load('./data/'+expnm+'_elm_et.npy') # (nt, ny, nx)
elm_pet= np.load('./data/'+expnm+'_elm_pet.npy') # (nt, ny, nx)
print(len(elm_pet))

## date ##
date_mod = pd.read_csv('./data/'+expnm+'_date.csv', na_values='-9999', parse_dates=True, skiprows=0)
date_mod.index = pd.to_datetime(date_mod['date'])
mod_doy = [pd.Period(ri,freq='D').dayofyear for ri in date_mod.index]
date_mod['doy'] = mod_doy
print(date_mod)

##### sim ELM #######
elm_day = pd.DataFrame() # mm/day
elm_day['doy'] = mod_doy
elm_day['pr'] = elm_pr
elm_day['et'] = elm_et
elm_day['pet'] = elm_pet
elm_day.index = pd.to_datetime(date_mod['date'])

elm_mon = elm_day.resample('1M', label='right').mean() # 1M, W-MON
#print(elm_mon)

##### sim well WTD #######
well_mod_day = pd.DataFrame()
for si in range(len(wellat)):
    well_mod_day[wellist[si]] = wtd_mod[:,welat_i[si],welon_i[si]]

well_mod_day.index = pd.to_datetime(date_mod['date'])
well_mod_day['doy'] = mod_doy
well_mod_mon = well_mod_day.resample('1M', label='right').mean() # 1M, W-MON
print(well_mod_day)
# print(well_mod_mon)

##### 2D average #####
date_idx = date_mod[ymdat]['step0'].tolist()
#print(date_idx)
sim_wtd_2D_avg = wtd_mod[date_idx,:,:].mean(axis=0)

######## select out the specific year to compare OBS #########
well_sim_day = well_mod_day[ymdat]
rive_sim_day = rive_mod_day[ymdat]

well_sim_mon = well_mod_mon[ymdat]
rive_sim_mon = rive_mod_mon[ymdat]
well_sim_avg = well_sim_mon.mean(axis=0)
