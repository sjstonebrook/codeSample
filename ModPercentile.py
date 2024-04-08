import os
import matplotlib.pyplot as plt
import numpy as np
import rasterio as rio
import xarray as xr
import rioxarray as rxr
from shapely.geometry import mapping, box
import geopandas as gpd
import glob
import warnings
import pandas as pd
import matplotlib.pyplot as plt

# this is just ignoring the warning that comes up when using np.nanpercentile when a slice has only Nan values
    # usually happens over water areas where there are no data values available
warnings.filterwarnings(action = 'ignore', message = 'All-NaN slice encountered')

mainPath = '/projectnb/modislc/users/sjstone/MODIS_LST_EVI/'

# function to compute the specified percentile for the specified product, variable, and year
# output: an array with dimensions 1200 x 1200 with a percentile value for each pixel
# input:
    # year: 2000 - 2022 (EVI)
    # year: 2002 - 2022 (LST)
    # product: can be 'MOD13A2.061' (EVI), 'MYD21A2.061' (temp)
    # selectedBands: for EVI ['1 km 16 days EVI'], for temperature ['LST_Day_1KM']
    # selectedPerc: 0-100, the desired percentile to compute 

def yearPercentile(year, product, selectedBands, selectedPerc):
    # selecting all of the files for the specified product and year
    yearFiles = glob.glob(mainPath + 'data/' + product + '/' + year + '*/*.hdf', recursive = True)
    
    prodList = []
    prodDate = []
    # reading in the selected files and appending them to a list
    for file in yearFiles:
        date = file.split('/')[8]
        fileRead = rxr.open_rasterio(file,
                                    masked = True,
                                    variable = selectedBands).squeeze()
        if(product == 'MYD21A2.061'):
            fileList = fileRead.to_array().values[0].flatten() * 0.02
        else:
            fileList = fileRead.to_array().values[0].flatten()
        
        prodList.append(fileList)
        prodDate.append(date)
        
    prodDF = pd.DataFrame(prodList) # generating a dataframe from the files
    prodDF.index = pd.to_datetime(prodDate, format = '%Y-%m-%d') # assigning the index as a date time
    prodDF = prodDF.sort_index() # sorting the index by the date
    
    prodInterp = prodDF.interpolate(method = 'time', limit_direction = 'both', axis = 0) # interpolating the NA values
    
    prodArrayList = []
    
    # creating an array for each of the dates 
    for row in list(range(0, prodInterp.shape[0])):
        prodSel = prodInterp.iloc[row]
        prodArray = np.array(prodSel).reshape(1200, 1200)
        
        prodArrayList.append(prodArray)
        
    prodStack = np.stack(prodArrayList) # generating a stack of the arrays from the different dates 
     
    percentile = np.nanpercentile(prodStack, selectedPerc, axis = 0) # taking the percentile 
    
    return(percentile)


# function that computed the fit coeficient for locations where there are no NA values 
# returns NA value where there are NA vlaues

def fit(x, yearVec):
    naSum = x.isna().sum()
    if(naSum > 0):
        coef = np.nan
    else:
        coef = np.polyfit(yearVec, x, 1)[0]
    return(coef)
    
infoFile = rxr.open_rasterio(glob.glob(mainPath + 'data/' + 'MYD21A2.061' + '/' + '2019' + '*/*.hdf', recursive = True)[0])

## LST Calculations

# medLSTPerc = []
selectedPerc = [10, 50, 90]


for perc in selectedPerc:
    medLSTPerc = []

    for yr in list(range(2003, 2022)):
        LSTperc = yearPercentile(str(yr), 'MYD21A2.061', ['LST_Day_1KM'], perc)
        medLSTPerc.append(LSTperc)

    medianLSTStack = np.stack(medLSTPerc)

    medianList = []
    for loc in list(range(0, medianLSTStack.shape[0])):
        flattened = medianLSTStack[loc, :, :].flatten()
        medianList.append(flattened)

    medianDF = pd.DataFrame(medianList)

    medianCoef = medianDF.apply(fit, axis = 0, yearVec = list(range(2003, 2022)))

    medianCoefarray = np.array(medianCoef).reshape(1200, 1200)
    
    rio.open(
        mainPath + 'data/coefficientArray/LSTCoefArray_' + str(perc) + '_percentile.tif',
        'w',
        height=medianCoefarray.shape[0],
        width=medianCoefarray.shape[1],
        count=1,
        dtype=medianCoefarray.dtype.name,
        crs=infoFile.rio.crs,
        transform=infoFile.rio.transform(),
        compress='lzw'
    ).write(medianCoefarray,1)
    print('LST percentile ' + str(perc) + ' done')
    
# EVI calculations    
# medEVIPerc = []
# selectedPerc = [10, 50, 90]

# for perc in selectedPerc:
#     medEVIPerc = []
    
#     for yr in list(range(2001, 2022)):
#         EVIperc = yearPercentile(str(yr), 'MOD13A2.061', ['1 km 16 days EVI'], perc)
#         medEVIPerc.append(EVIperc)
#         print(str(yr) + ' percentile done', end = ' ')

#     medianEVIstack = np.stack(medEVIPerc)

#     medianEVIlist = []
#     for loc in list(range(0, medianEVIstack.shape[0])):
#         flattened = medianEVIstack[loc, :, :].flatten()
#         medianEVIlist.append(flattened)

#     medianEVIdf = pd.DataFrame(medianEVIlist)

#     print('median data frame generated')

#     medianEVICoef = medianEVIdf.apply(fit, axis = 0, yearVec = list(range(2001, 2022)))

#     print('coefficients calculated')

#     medianEVICoefArray = np.array(medianEVICoef).reshape(1200, 1200)

#     print('coef array generated')
#     rio.open(
#         mainPath + 'data/coefficientArray/EVICoefArray_' + str(perc) + '_percentile.tif',
#         'w',
#         height=medianEVICoefArray.shape[0],
#         width=medianEVICoefArray.shape[1],
#         count=1,
#         dtype=medianEVICoefArray.dtype.name,
#         crs=infoFile.rio.crs,
#         transform=infoFile.rio.transform(),
#         compress='lzw'
#     ).write(medianEVICoefArray,1)
#     print('EVI percentile ' + str(perc) + ' done')
    
