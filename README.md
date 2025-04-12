# Exponential-average-growth
This project is used for modification of the standard reaction rate in finite element package Febio.
Original diffusion like kernel has been non-invariant formulation like
$$\frac{d\rho}{dt}=\sum{\frac{d\rho_{i}}{dt}exp(-\frac{d}{D})}$$.
Indexing is done for initial neighorhood elements.
Small d is used for element centroid distance and big D as the material parameter.
For an invariant formulation, this operation can be rewritten with an integral operation using the initial element volume like
$$\frac{d\rho}{dt}=(\sum{\frac{d\rho_{i}}{dt}exp(-\frac{d}{D})}V_{i})/\sum{V_{i}}$$.
For build, the plugin should be using CMake to specify Febio SDK folder that is placed in a similar folder with Febio Studio.
For using the plugin as a default, you should add DLL name to febio.xml in bin folder.
This reaction rate is replaces the original in bone growth material with the same material properties.
Example plate used plane strain formulation with X varied load.
(example images%plate_fixed.png)
![example](https://github.com/Panavia-Tornado/Exponential-average-growth/blob/main/example%20images/plate_fixed.png)
Result density distribution without parameter D and included over a long period of time(left and right part of image accordingly).
![res](https://github.com/Panavia-Tornado/Exponential-average-growth/blob/main/example%20images/result.jpg)
