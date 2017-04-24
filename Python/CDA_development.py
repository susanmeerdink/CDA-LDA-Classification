# Created By: Susan Meerdink
# 4/17/17
# This python code will develop canonical discriminant analysis coefficients/eigenvectors for spectral library
# --------------------------------------------------------------------------------------------------------------------
# Inputs
libLocation = 'F:\\Classification-Products\\2 - Statistical Analysis\\'  # Uses Log Transformed spectra
outLocation = 'F:\\Classification-Products\\3 - CDA Development\\'
dateTag = '140416'

from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn import linear_model
import numpy as np
import random
from sets import Set

# Read in spectra and metadata files
libSpecCalFile = libLocation + dateTag + '_transformed_spectral_library_spectra.csv'
libMetaCalFile = libLocation + dateTag + '_transformed_spectral_library_metadata.csv'
libSpecValFile = libLocation + dateTag + '_transformed_spectral_library_spectra.csv'
libMetaValFile = libLocation + dateTag + '_transformed_spectral_library_metadata.csv'
metaCal = np.loadtxt(libMetaCalFile, dtype=object, delimiter=',')
metaVal = np.loadtxt(libMetaValFile, dtype=object, delimiter=',')
spectraCal = np.loadtxt(libSpecCalFile, dtype=object, delimiter=',')  # Load in spectra - skips first line
metaSpecCal = spectraCal[:, 0:5]  # save first 5 columns of spectra separately
spectraCal = np.delete(spectraCal, [0, 1, 2, 3, 4], 1)  # remove the 5 columns of metadata in spectra
spectraCal = spectraCal.astype(np.double)  # convert from string array to double array
spectraCal = np.nan_to_num(spectraCal)
spectraVal = np.loadtxt(libSpecValFile, dtype=object, delimiter=',')  # Load in spectra - skips first line
metaSpecVal = spectraVal[:, 0:5]  # save first 5 columns of spectra separately
spectraVal = np.delete(spectraVal, [0, 1, 2, 3, 4], 1)  # remove the 5 columns of metadata in spectra
spectraVal = spectraVal.astype(np.double)  # convert from string array to double array
spectraVal = np.nan_to_num(spectraVal)
wavelengths = [365.913666,  375.577667,  385.246674,  394.919647,  404.596649,  414.278656,   423.964661,  433.654663,
                 443.349670,  453.049652,  462.752655,  472.460663,  482.173645,  491.890656,  501.611664,  511.337646,
                 521.067688, 530.801636,  540.540649,  550.283630,  560.031677,  569.783630,  579.539673,  589.300659,
                 599.065674,  608.835632,  618.608643,  628.387634,  638.169678,  647.957642,  657.748657,  667.544678,
                 655.475647,  665.282654, 675.084656,  684.881653,  694.672668,  704.459656,  714.240662,  724.016663,
                 733.786682,  743.552673,  753.312683,  763.067688,  772.816650,  782.561646,  792.300659,  802.035645,
                 811.764648,  821.487671,  831.206665, 840.919678,  850.627686,  860.330688,  870.028687,  879.720642,
                 889.407654,  899.090637,  908.766663,  918.438660,  928.104675,  937.766663,  947.422668,  957.072632,
                 966.718689,  976.358643,  985.994629,  995.624634,1005.253662, 1014.863647, 1024.483643, 1034.093628,
                 1043.693604, 1053.293701, 1062.883667, 1072.473633, 1082.063599, 1091.633667, 1101.213623, 1110.773682,
                 1120.343628, 1129.893677, 1139.443604, 1148.993652, 1158.533691, 1168.073608, 1177.603638, 1187.133667,
                 1196.653687, 1206.163696, 1215.673706, 1225.183716, 1234.683716, 1244.173706, 1253.663696, 1263.143677,
                 1253.353638, 1263.333618, 1273.303711, 1283.273682, 1293.243652,1303.213623, 1313.193604, 1323.163696,
                 1333.133667, 1343.103638, 1353.073608, 1363.043701, 1373.013672, 1382.983643, 1392.953613, 1402.923706,
                 1412.893677, 1422.863647, 1432.833618, 1442.793701, 1452.763672,1462.733643, 1472.703613, 1482.663696,
                 1492.633667, 1502.603638, 1512.573608, 1522.533691, 1532.503662, 1542.463623, 1552.433716, 1562.403687,
                 1572.363647, 1582.333618, 1592.293701, 1602.263672, 1612.223633,1622.183716, 1632.153687, 1642.113647,
                 1652.073608, 1662.043701, 1672.003662, 1681.963623, 1691.933716, 1701.893677, 1711.853638, 1721.813599,
                 1731.773682, 1741.733643, 1751.693604, 1761.663696, 1771.623657,1781.583618, 1791.543701, 1801.503662,
                 1811.453613, 1821.413696, 1831.373657, 1841.333618, 1851.293701, 1861.253662, 1871.213623, 1872.363647,
                 1866.843628, 1876.913696, 1886.963623, 1897.023682, 1907.083618,1917.133667, 1927.183716, 1937.233643,
                 1947.273682, 1957.313599, 1967.363647, 1977.393677, 1987.433716, 1997.463623, 2007.503662, 2017.523682,
                 2027.553711, 2037.583618, 2047.603638, 2057.623779, 2067.643555,2077.653564, 2087.663574, 2097.683594,
                 2107.683594, 2117.693604, 2127.693604, 2137.703613, 2147.693604, 2157.693604, 2167.693604, 2177.683594,
                 2187.673584, 2197.663574, 2207.643555, 2217.633545, 2227.613770,2237.583740, 2247.563721, 2257.533691,
                 2267.513672, 2277.483643, 2287.443604, 2297.413574, 2307.373779, 2317.333740, 2327.293701, 2337.243652,
                 2347.203613, 2357.153564, 2367.093750, 2377.043701, 2386.993652,2396.933594, 2406.873779, 2416.803711,
                 2426.743652, 2436.673584, 2446.603760, 2456.533691, 2466.453613, 2476.383545, 2486.303711, 2496.223633]

# Develop canonical discriminant variables
clf = LinearDiscriminantAnalysis()  # http://scikit-learn.org/stable/modules/generated/sklearn.discriminant_analysis.LinearDiscriminantAnalysis.html#sklearn.discriminant_analysis.LinearDiscriminantAnalysis
numIterations = 10
cdaValues = np.empty([24, 224])
cdaResults = np.empty([numIterations, 4])  # Array for intercept, slope, r2, and rmse
for i in range(0, numIterations):
    internalCalIndex = random.sample(xrange(0, spectraCal.shape[0]), int(round(0.7 * spectraCal.shape[0])))  # Select 70% of sample for internal calibration and development of canonical variables
    internalValIndex = list(Set(range(0, spectraCal.shape[0])).difference(internalCalIndex))
    calX = spectraCal[internalCalIndex, :]
    calY = np.array(metaCal[internalCalIndex, 14]).astype(np.int)
    valX = spectraCal[internalValIndex, :]
    valY = np.array(metaCal[internalValIndex, 14]).astype(np.int)
    clf.fit(calX, calY)
    predictResults = clf.predict(valX)
    cdaValues = np.add(cdaValues, clf.coef_)

    # Calculate results from CDA development
    regr = linear_model.LinearRegression()
    linearResults = regr.fit(valY.reshape(len(valY), 1), predictResults.reshape(len(predictResults), 1), sample_weight=None)  # http://scikit-learn.org/stable/modules/generated/sklearn.linear_model.LinearRegression.html
    cdaResults[i, 0] = linearResults.intercept_
    cdaResults[i, 1] = linearResults.coef_
    cdaResults[i, 2] = linearResults.score(valY.reshape(len(valY), 1), predictResults.reshape(len(predictResults), 1), sample_weight=None)
    cdaResults[i, 3] = np.mean((regr.predict(valY.reshape(len(valY), 1)) - predictResults.reshape(len(predictResults), 1)) ** 2)

avgCDA = cdaValues / numIterations
valResults = clf.predict(spectraVal)
regr = linear_model.LinearRegression()
linearResults = regr.fit(metaVal[:, 14].reshape(len(metaVal[:, 14]), 1), valResults.reshape(len(valResults), 1))  # http://scikit-learn.org/stable/modules/generated/sklearn.linear_model.LinearRegression.html
print linearResults.score(metaVal[:, 14].reshape(len(metaVal[:, 14]), 1), valResults.reshape(len(valResults), 1))
tempLib = np.dot(spectraVal, avgCDA.T)

# Save Canonical Variables to Text File
outCDA = file(outLocation + dateTag + '_CDA_spectral_library_spectra.csv', 'wb')
np.savetxt(outCDA, avgCDA, fmt='%s', delimiter=',')
outCDA.close()

# Save linear regression results from Canonical Variables to Text File
outCDA = file(outLocation + dateTag + '_CDA_development_results.csv', 'wb')
cdaWrite = np.vstack((np.array(['Intercept', 'Slope', 'R^2', 'RMSE']), cdaResults))
np.savetxt(outCDA, cdaWrite, fmt='%s', delimiter=',')
outCDA.close()