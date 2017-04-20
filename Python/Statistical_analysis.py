# Created By: Susan Meerdink
# 4/17/17
# This python code will read in spectral libraries and determine some basic statistics for future analysis
# --------------------------------------------------------------------------------------------------------------------
# Inputs
libLocation = 'F:\\Classification-Products\\1 - Spectral Library\\Combined Single Date\\'
outLocation = 'F:\\Classification-Products\\2 - Statistical Analysis\\'
dateTag = '140416'

# Import Modules
import numpy as np
import scipy.stats
import matplotlib.pyplot as plt
from statsmodels.stats.multicomp import pairwise_tukeyhsd
from statsmodels.stats.multicomp import MultiComparison

# Get metadata and spectra files and read them in for full library
libSpecFile = libLocation + dateTag + '_spectral_library_spectra.csv'
libMetaFile = libLocation + dateTag + '_spectral_library_metadata.csv'
metadata = np.loadtxt(libMetaFile, dtype=object, delimiter=',', skiprows=0)  # Load in metadata
# Right now the code above is not reading in first line of file
# headers = metadata[0, :]  # save headers separate of metadata
# headers = np.char.strip(headers.astype(str))  # remove whitespace from headers
# metadata = np.delete(metadata, 0, 0)  # remove the headers
spectra = np.loadtxt(libSpecFile, dtype=object, delimiter=',')  # Load in spectra - skips first line
spectraMeta = spectra[:, 0:5]  # save first 5 columns of spectra separately
spectra = np.delete(spectra, [0, 1, 2, 3, 4], 1)  # remove the 5 columns of metadata in spectra
spectra = spectra.astype(np.double)  # convert from string array to double array
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

# Loop through wavelengths and calculate whether sample belongs from normal distribution
normalStats = np.empty([224])
kurtosisStats = np.empty([224])
skewStats = np.empty([224])
for w in range(0, 224):
    # Determine if the sample comes from a normal distribution (pvalue < 0.05 data is not from normal distribution)
    stats = scipy.stats.normaltest(spectra[:, w])  # https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.normaltest.html#scipy.stats.normaltest
    kurtosis = scipy.stats.kurtosis(spectra[:, w], nan_policy='omit')  # https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.kurtosistest.html#scipy.stats.kurtosistest
    skew = scipy.stats.skew(spectra[:, w], nan_policy='omit')  # https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.skewtest.html#scipy.stats.skewtest
    normalStats[w] = stats.pvalue
    kurtosisStats[w] = kurtosis
    skewStats[w] = skew

# Transform data and calculate whether sample now belongs to normal distribution
with np.errstate(all='ignore'):  # zero values in spectra will produce error, ignore them
    spectraTransform = np.log(spectra)
    spectraTransform[np.isneginf(spectraTransform)] = 0
normalTransformStats = np.empty([224])
kurtosisTransformStats = np.empty([224])
skewTransformStats = np.empty([224])
for w in range(0, 224):
    # Determine if the sample comes from a normal distribution now that spectra is transformed
    stats = scipy.stats.normaltest(spectraTransform[:, w], nan_policy='omit')
    kurtosis = scipy.stats.kurtosis(spectraTransform[:, w], nan_policy='omit')
    skew = scipy.stats.skew(spectraTransform[:, w], nan_policy='omit')
    normalTransformStats[w] = stats.pvalue
    kurtosisTransformStats[w] = kurtosis
    skewTransformStats[w] = skew

# Read in spectra files for calibration and validation for Transformation
libSpecCalFile = libLocation + dateTag + '_spectral_library_calibration_spectra.csv'
libSpecValFile = libLocation + dateTag + '_spectral_library_validation_spectra.csv'
spectraCal = np.loadtxt(libSpecCalFile, dtype=object, delimiter=',')  # Load in spectra - skips first line
spectraMetaCal = spectraCal[:, 0:5]  # save first 5 columns of spectra separately
spectraCal = np.delete(spectraCal, [0, 1, 2, 3, 4], 1)  # remove the 5 columns of metadata in spectra
spectraCal = spectraCal.astype(np.double)  # convert from string array to double array
spectraVal = np.loadtxt(libSpecValFile, dtype=object, delimiter=',')  # Load in spectra - skips first line
spectraMetaVal = spectraVal[:, 0:5]  # save first 5 columns of spectra separately
spectraVal = np.delete(spectraVal, [0, 1, 2, 3, 4], 1)  # remove the 5 columns of metadata in spectra
spectraVal = spectraVal.astype(np.double)  # convert from string array to double array
with np.errstate(all='ignore'):  # zero values in spectra will produce error, ignore them
    spectraTransformCal = np.log(spectraCal)
    spectraTransformVal = np.log(spectraVal)
    spectraTransformCal[np.isneginf(spectraTransformCal)] = 0
    spectraTransformVal[np.isneginf(spectraTransformVal)] = 0

# Run multiple anova on new transformed dataset
dominant = metadata[:, 13]
listDominant = np.unique(dominant)  # Next code is based on the fact there are 24 classes
anovaResults = np.empty([224, 2])
tukeyResults = np.empty([0, 8])
outTukey = file(outLocation + dateTag + '_tukeyhsd_results.csv', 'wb')
for w in range(0, 224):
    if any(spectraTransform[:, w]):
        anovaResults[w, 0], anovaResults[w, 1] = scipy.stats.f_oneway(
                                        spectraTransform[np.where(dominant == listDominant[0])[0], w],
                                        spectraTransform[np.where(dominant == listDominant[1])[0], w],
                                        spectraTransform[np.where(dominant == listDominant[2])[0], w],
                                        spectraTransform[np.where(dominant == listDominant[3])[0], w],
                                        spectraTransform[np.where(dominant == listDominant[4])[0], w],
                                        spectraTransform[np.where(dominant == listDominant[5])[0], w],
                                        spectraTransform[np.where(dominant == listDominant[6])[0], w],
                                        spectraTransform[np.where(dominant == listDominant[7])[0], w],
                                        spectraTransform[np.where(dominant == listDominant[8])[0], w],
                                        spectraTransform[np.where(dominant == listDominant[9])[0], w],
                                        spectraTransform[np.where(dominant == listDominant[10])[0], w],
                                        spectraTransform[np.where(dominant == listDominant[11])[0], w],
                                        spectraTransform[np.where(dominant == listDominant[12])[0], w],
                                        spectraTransform[np.where(dominant == listDominant[13])[0], w],
                                        spectraTransform[np.where(dominant == listDominant[14])[0], w],
                                        spectraTransform[np.where(dominant == listDominant[15])[0], w],
                                        spectraTransform[np.where(dominant == listDominant[16])[0], w],
                                        spectraTransform[np.where(dominant == listDominant[17])[0], w],
                                        spectraTransform[np.where(dominant == listDominant[18])[0], w],
                                        spectraTransform[np.where(dominant == listDominant[19])[0], w],
                                        spectraTransform[np.where(dominant == listDominant[20])[0], w],
                                        spectraTransform[np.where(dominant == listDominant[21])[0], w],
                                        spectraTransform[np.where(dominant == listDominant[22])[0], w],
                                        spectraTransform[np.where(dominant == listDominant[23])[0], w])
        # If the anova turns back a pvalue < 0.05, do multicomparison to figure out what samples are different
        if anovaResults[w, 1] < 0.05:
            mc = MultiComparison(spectraTransform[:, w], dominant)  # http://statsmodels.sourceforge.net/0.6.0/_modules/statsmodels/stats/multicomp.html
            result = mc.tukeyhsd()  # http://statsmodels.sourceforge.net/devel/generated/statsmodels.sandbox.stats.multicomp.MultiComparison.tukeyhsd.html
            inResults = np.array([mc.groupsunique[mc.pairindices[0]], mc.groupsunique[mc.pairindices[1]], result.meandiffs, result.confint[:, 0], result.confint[:, 1], result.std_pairs, result.reject]).T
            inResults = np.column_stack((np.repeat(wavelengths[w], len(result.reject)), inResults))
            tukeyResults = np.vstack((tukeyResults, inResults))

# Set up csv file to output statistical results
outStats = file(outLocation + dateTag + '_statistical_analysis.csv', 'a')  # Opening in append mode
row1 = np.hstack(('normal distribution p value for original spectra', normalStats))
row2 = np.hstack(('kurtosis p value for original spectra', kurtosisStats))
row3 = np.hstack(('skew p value for original spectra', skewStats))
row4 = np.hstack(('normal distribution p value for transformed spectra', normalTransformStats))
row5 = np.hstack(('kurtosis p value for transformed spectra', kurtosisTransformStats))
row6 = np.hstack(('skew p value for transformed spectra', skewTransformStats))
row7 = np.hstack(('anova results for transformed spectra', anovaResults[:, 1]))
inRows = np.vstack((row1, row2, row3, row4, row5, row6, row7))
np.savetxt(outStats, inRows, fmt='%s', delimiter=',')
outStats.close()

# Set up csv file to output transform spectra
outSpec = file(outLocation + dateTag + '_transformed_spectral_library_spectra.csv', 'wb')
outSpecCal = file(outLocation + dateTag + '_transformed_spectral_library_calibration_spectra.csv', 'wb')
outSpecVal = file(outLocation + dateTag + '_transformed_spectral_library_validation_spectra.csv', 'wb')
headerOutSpec = 'Flightline, Date, PolygonName, X, Y,' + ','.join(map(str, wavelengths))
allSpec = np.hstack((spectraMeta, spectraTransform))
np.savetxt(outSpec, allSpec, header=headerOutSpec, fmt='%s', delimiter=",")
allSpecCal = np.hstack((spectraMetaCal, spectraTransformCal))
np.savetxt(outSpecCal, allSpecCal, header=headerOutSpec, fmt='%s', delimiter=",")
allSpecVal = np.hstack((spectraMetaVal, spectraTransformVal))
np.savetxt(outSpecVal, allSpecVal, header=headerOutSpec, fmt='%s', delimiter=",")
outSpec.close()
outSpecCal.close()
outSpecVal.close()

# Set up csv file to output Tukey Results
writeTukey = np.vstack((np.array(['Wavelength', 'Group1', 'Group2', 'MeanDiff', 'Confidence_Interval', 'Confidence_Interval', 'Std_Pairs', 'Reject']), tukeyResults))
np.savetxt(outTukey, writeTukey, fmt='%s', delimiter=',')
outTukey.close()