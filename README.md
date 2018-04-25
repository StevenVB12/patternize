# <p> patternize - An R package for quantifying &#13;&#10;color pattern variation <img src="https://cloud.githubusercontent.com/assets/6349171/22620648/29ecb77e-eb08-11e6-8f7e-80d3a3807fda.png" alt="patternize" width="300" align="right"></p>

Quantifying variation in color patterns to study and compare the consistency of their expression necessitates the homologous alignment and color-based segmentation of images. Patternize is an R package that quantifies variation in color patterns as obtained from image data. Patternize defines homology between pattern positions across specimens either through fixed landmarks or image registration. Pattern identification is performed by categorizing the distribution of colors using either an RGB threshold or an unsupervised image segmentation. The quantification of the color patterns can be visualized as heat maps and compared between sets of samples.

```diff
Please do not hesitate to contact me with any questions or suggestions!
```

<b>Install patternize in R from CRAN (not up to date at this time!!):</b>

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

<b>Installation errors</b>

Some people have noted platform specific installation errors (mostly Mac). If you don't find your solution here, please contact me.

```
error: unable to load shared object 
'/Library/Frameworks/R.framework/Versions/3.4/Resources/library/rgl/libs/rgl.so'

Solution: download XQuartz https://www.xquartz.org 
```




<b>For examples see package examples or https://github.com/StevenVB12/patternize-examples</b>

<b>Workflow</b>

<img src="https://user-images.githubusercontent.com/6349171/27639019-f1a91760-5c0c-11e7-9705-40a7a24700b1.png" alt="workflow" width="500" align="center"></p>

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

If you are working with <i>Heliconius</i> or related butterflies, consider using this landmark scheme:

<img src="https://user-images.githubusercontent.com/6349171/38930020-44607508-4306-11e8-8449-0363040adabe.png" alt="landmarks_Heliconius_FW" width="500" align="center"></p>
<img src="https://user-images.githubusercontent.com/6349171/38930021-44746eaa-4306-11e8-995b-34096c59dc4f.png" alt="landmarks_Heliconius_HW" width="500" align="center"></p>

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

<b>Calibrate images using ColorChecker</b>

The ColorChecker (Macbeth ColorChecker) is a color calibration target consisting of a framed arrangement of 24 squares of painted samples. The chart patches are selected to mimic those of natural objects and to have consistent color appearance under a variety of lighting conditions. By including this chart in an image, it can be used to calibrate images taken in different lightning conditions.

<img src="https://user-images.githubusercontent.com/6349171/37557345-0a13112c-29fb-11e8-9d55-ebf4ba21a0ea.png" alt="ColorChecker calibrated" width="500" align="center"></p>

You can calibrate images like this:

```
IDlist <- c('image1', 'image2')
extension <- '.jpg'
colorChecker(IDlist, extension)
```

The function is also available as a stand alone app: https://stevenvanbelleghem.shinyapps.io/shiny_colorchecker/ (see https://github.com/StevenVB12/shiny_ColorChecker)

This calculates a second order polynomial regression between the observed and expected RGB values of the 24 ColorChecker patches and performs the calibration of the image. The function will ask you to define the corners of the ColorChecker like this (also works for images taken with a skew angle):

<img src="https://user-images.githubusercontent.com/6349171/37557207-35132aee-29f9-11e8-84b7-04b31850cea1.png" alt="ColorChecker corners" width="200" align="center"></p>

