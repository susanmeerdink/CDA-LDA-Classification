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

# Get library files and read them in
libSpecFile = libLocation + dateTag + '_spectral_library_spectra.csv'
libMetaFile = libLocation + dateTag + '_spectral_library_metadata.csv'
metadata = np.loadtxt(libMetaFile, dtype=object, delimiter=',')  # Load in metadata
headers = metadata[0, :]  # save headers separate of metadata
headers = np.char.strip(headers.astype(str))  # remove whitespace from headers
metadata = np.delete(metadata, 0, 0)  # remove the headers
spectra = np.loadtxt(libSpecFile, dtype=object, delimiter=',')  # Load in metadata
spectraMeta = spectra[0:4, :]  # save first 5 columns of spectra separately
spectra = np.delete(spectra, [0, 1, 2, 3, 4], 1)  # remove the 5 columns of metadata in spectra
spectra = spectra.astype(np.double)  # convert from string array to double array

# Loop through wavelengths and calculate whether sample belongs from normal distribution
normalStats = np.empty([224])
kurtosisStats = np.empty([224])
skewStats = np.empty([224])
for w in range(0, 224):
    # Determine if the sample comes from a normal distribution (pvalue < 0.05 data is not from normal distribution)
    stats = scipy.stats.normaltest(spectra[:, w])  # https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.normaltest.html#scipy.stats.normaltest
    kurtosis = scipy.stats.kurtosistest(spectra[:, w], nan_policy='omit')  # https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.kurtosistest.html#scipy.stats.kurtosistest
    skew = scipy.stats.skewtest(spectra[:, w], nan_policy='omit')  # https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.skewtest.html#scipy.stats.skewtest
    normalStats[w] = stats.pvalue
    kurtosisStats[w] = kurtosis.pvalue
    skewStats[w] = skew.pvalue

# Transform data and calculate whether sample now belongs to normal distribution
with np.errstate(divide='ignore'):
    spectraTransform = np.log(spectra)
print spectraTransform[5,:]
normalTransformStats = np.empty([224])
kurtosisTransformStats = np.empty([224])
skewTransformStats = np.empty([224])
for w in range(0, 224):
    # Determine if the sample comes from a normal distribution now that spectra is transformed
    stats = scipy.stats.normaltest(spectraTransform[:, w], nan_policy='omit')
    kurtosis = scipy.stats.kurtosistest(spectraTransform[:, w], nan_policy='omit')
    skew = scipy.stats.skewtest(spectraTransform[:, w], nan_policy='omit')
    normalTransformStats[w] = stats.pvalue
    kurtosisTransformStats[w] = kurtosis.pvalue
    skewTransformStats[w] = skew.pvalue

# Set up csv file to output results
outStats = file(outLocation + dateTag + '_statistical_analysis.csv', 'a')  # Opening in append m
row1 = np.hstack(('normal distribution p value for original spectra', normalStats))
row2 = np.hstack(('kurtosis p value for original spectra', kurtosisStats))
row3 = np.hstack(('skew p value for original spectra', skewStats))
row4 = np.hstack(('normal distribution p value for transformed spectra', normalTransformStats))
row5 = np.hstack(('kurtosis p value for transformed spectra', kurtosisTransformStats))
row6 = np.hstack(('skew p value for transformed spectra', skewTransformStats))
inRows = np.vstack((row1, row2, row3, row4, row5, row6))
np.savetxt(outStats, inRows, fmt='%s', delimiter=',')
# pvalue is still small but distribution is closer to normal

# Run multiple anova on new transformed dataset
dominant = metadata[:, 13]
listDominant = np.unique(dominant)
for w in range(0, 224):
    if any(spectraTransform[:, w]):
        f, p = scipy.stats.f_oneway(spectraTransform[np.where(dominant == listDominant[0])[0], w],
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
        print ('One-way ANOVA')
        print ('=============')

        print ('F value:', f)
        print ('P value:', p, '\n')

        mc = MultiComparison(spectraTransform[:, w], dominant)
        result = mc.tukeyhsd()

        print(result)
        print(mc.groupsunique)


