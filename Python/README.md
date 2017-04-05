# CDA-LDA-Classification

Project related to CDA and LDA classification (using only python)

### Dependencies
---
glob: https://pymotw.com/2/glob/
gdal: https://pypi.python.org/pypi/GDAL
numpy: 


### Inputs
---
1. Polygons: These are polygons of plant species across the study area
2. Image: This is a flight of AVIRIS or AVIRIS+MASTER imagery from the study area

### Code Order
---
1. Developing Spectral Libraries
	a. Extract spectra from image files: FL02 - FL11.
	b. Statistics - anova
	c. Split into training and validation 
2. Develop Canonical Discriminant Analysis coefficients/egienvectors
	a. Run CDA on spectral library
	b. Run accuracy analysis
	c. Apply CDA coefficients to image
3. Run Linear Discriminant Analysis Classification
	a. Run LDA on imagery
	b. Run accuracy analysis/ classification results
