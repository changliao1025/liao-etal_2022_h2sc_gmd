import numpy as np

iFlag_obs = 1
iFlag_sim = 1

if iFlag_obs == 1:
    sFilename = sWorkspace_auxiliary + slash + 'situ' + slash + 'INPA-LBA_WellData_2001_2016.xlsx'
    xl = pd.ExcelFile(sFilename)
    aSheet = xl.sheet_names  # see all sheet names 14, last one is summary
    aElevation = np.array([59, 59, 60, 61, 60, 87,87, 81, 101, 96, 101,101, 101 ])
    aDate_host=list()
    nyear = iYear_end - iYear_start + 1
    for iYear in range(iYear_start, iYear_end + 1):
        for iMonth in range(1,13):
            dom = day_in_month(iYear, iMonth)
            for iDay in range(1, dom+1):
                dSimulation = datetime.datetime(iYear, iMonth, iDay)
                aDate_host.append( dSimulation )
                pass

    aDate_host=np.array(aDate_host)
    nobs_host = len(aDate_host)
    aData_host = np.full( (13,nobs_host), np.nan, dtype=float)
    # we skip some data because language is not in english
    #4 of them not used
    aFlag = np.full( 13, 0, dtype=int )
    for iSheet in np.arange(1,14, 1):
        sSheet = aSheet[iSheet-1]
        if sSheet == 'PZ_PT-07':
            continue
        if sSheet == 'PP03':
            continue
        if sSheet == 'PP2':
            continue
        if sSheet == 'PP01':
            continue

        print(sSheet)
        aFlag[iSheet -1 ] =1
        df = pd.read_excel(sFilename, \
                           sheet_name=sSheet, \
                           header=None, \
                           skiprows=range(5), \
                           usecols='A,E')
        df.columns = ['Date','WTD']
        dummy1 = df['Date']
        dummy2 = np.array(dummy1)
        dummy3 = dt2cal(dummy2)
        #aDate_obs = pd.to_datetime(np.array(dummy1))
        nobs =len(dummy3)
        aDate_obs = list()
        for iObs in range(nobs):
            d1=dummy3[iObs,0]
            d2=dummy3[iObs,1]
            d3= dummy3[iObs,2]
            dummy4= datetime.datetime(d1,d2 , d3 )
            aDate_obs.append( dummy4 )
            pass

        aDate_obs= np.array(aDate_obs)
        dummy5 = df['WTD']
        aWTD_obs_dummy = np.array(dummy5)  # mg/l

        #now fit the data inside the host

        #the existing data is outside limit,

        dummy_index = aDate_obs-aDate_host[0]
        #dummy_index1 = dummy_index[0].days
        dummy_index1 = [x.days for x in dummy_index ]
        dummy_index1 = np.array(dummy_index1)
        dummy_index2 = np.where( dummy_index1 < nobs_host )

        dummy_obs= aWTD_obs_dummy[dummy_index2]
        dummy_index3 = dummy_index1[dummy_index2]
        aData_host[iSheet-1, dummy_index3 ] = dummy_obs
        pass

    #need index

    #remove unused sites
    dummy = np.where( aFlag == 1)
    aElevation1= aElevation[dummy]

    aElevation2 , indices = np.unique(aElevation1, return_index=True)
    #aOrder  = np.argsort(aElevation)
    #aElevation_sort = np.sort(aElevation2)
    dummy2 = dummy[0][indices]
    aData_host2 = aData_host[dummy2, :  ]



    #get the average
    aWTD_obs = np.nanmean(aData_host2, axis=1)
    aWT_obs = aElevation2 - aWTD_obs