# <p> patternize - An R package for quantifying &#13;&#10;color pattern variation <img src="https://cloud.githubusercontent.com/assets/6349171/22620648/29ecb77e-eb08-11e6-8f7e-80d3a3807fda.png" alt="patternize" width="300" align="right"></p>

Quantifying variation in color patterns to study and compare the consistency of their expression necessitates the homologous alignment and color-based segmentation of images. Patternize is an R package that quantifies variation in color patterns as obtained from image data. Patternize defines homology between pattern positions across specimens either through fixed landmarks or image registration. Pattern identification is performed by categorizing the distribution of colors using either an RGB threshold or an unsupervised image segmentation. The quantification of the color patterns can be visualized as heat maps and compared between sets of samples.

```diff
Please do not hesitate to contact me with any questions or suggestions!
```

<b>Install patternize in R from CRAN:</b>

```
install.packages("patternize")
```

<b>Install patternize in R from GitHub using devtools:</b>

```
install.packages("devtools")
library(devtools)
install_github("StevenVB12/patternize")
# Morpho has committed a change that affects computeTransform used in
# the patternize functions patLanRGB and patLanK (not available in CRAN)
install_github("zarquon42b/Morpho")
```

<b>For examples see package examples or https://github.com/StevenVB12/patternize-examples</b>

<b>Workflow</b>

<img src="https://cloud.githubusercontent.com/assets/6349171/24803941/bf7e7274-1ba4-11e7-87e3-fdbc37b1c512.png" alt="workflow" width="500" align="center"></p>

<b>Setting landmarks</b>

I set landmarks using the Fiji distribution of ImageJ (https://fiji.sc/):

```
Fiji
> File > open (image)
> Multi-point (make sure points are set in the same order for each sample)
> Save As > XY coordinates
```

But you can also use (setting a resampleFactor > 0 speeds up things, but reduces resolution):

```
landmarkList <- sampleLandmarks(sampleList, resampleFactor = NULL, crop = c(0,0,0,0))
```

<b>Setting RGB value</b>

You can assign an RGB vector manually (e.g. RGB <- c(255,0,0) for red) or use:

```
RGB <- sampleRGB(image, resampleFactor = NULL, crop = c(0,0,0,0))
```

<b>Making cartoon for plotting (or masking)</b>

To plot a cartoon of the organism or trait of interest I use XY coordinates of an outline or lines obtained in the Fiji distribution of ImageJ (https://fiji.sc/). This is an annoying manual task, but if done precisely it can provide an outline or cartoon for any type of shape. The cartoon should be drawn for the reference (target) image when using image registration (patRegRGB or patRegK) or when using landmark transformation to a target image (patLanRGB or patLanK; transformRef = 'sample_ID'). When using transformRef = 'meanshape', the cartoon will also be transformed to the mean shape.

outline

```
Fiji
> File > open (image)
> Polygon selections (draw outline)
> Save As > XY coordinates
```

lines (same as setting landmarks)

```
Fiji
> File > open (image)
> Multi-point (draw lines by setting points)
> Save As > XY coordinates
```
