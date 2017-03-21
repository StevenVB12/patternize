# <p> patternize - An R package for quantifying &#13;&#10;color pattern variation <img src="https://cloud.githubusercontent.com/assets/6349171/22620648/29ecb77e-eb08-11e6-8f7e-80d3a3807fda.png" alt="patternize" width="300" align="right"></p>

Quantifying variation in color patterns to study and compare the consistency of their expression necessitates the homologous alignment and color-based segmentation of images. Patternize is an R package that quantifies variation in color patterns as obtained from image data. Patternize defines homology between pattern positions across specimens either through fixed landmarks or image registration. Pattern identification is performed by categorizing the distribution of colors using either an RGB threshold or an unsupervised image segmentation. The quantification of the color patterns can be visualized as heat maps and compared between sets of samples.

<b>Install patternize in R from CRAN:</b>

```
install.packages("patternize")
```

<b>Install patternize in R from GitHub using devtools:</b>

```
install.packages("devtools")
library(devtools)
install_github("StevenVB12/patternize")
# Morpho has commited a change that affects computeTransform used in
# the patternize functions patLanRGB and patLanK (not available in CRAN)
install_github("zarquon42b/Morpho")
```

<b>For examples see package examples or https://github.com/StevenVB12/patternize-examples</b>
