# LeArEst 1.0.0

* improved algorithm for length and area estimation
* replaced all occurrences of "variance" to "standard deviation" through the whole package
* added support for Student error distribution (T1,..., T5)

# LeArEst 0.2.0

* added support for multithreaded calculation of area
* added support for estimating object as a circle
* added new slicing method ("star") in area estimation functions
* new imports: doParallel, foreach, parallel
* new parameters in areaest function (parallel = FALSE, 
  slicing = c("hv", "star"), representation = c("ellipse", "circle"))
* new options in area estimation web interface ("Slicing", "Parallelization",
  "Represent object as")

# LeArEst 0.1.5

* added support for opencpu 2.x.x

# LeArEst 0.1.4

* initial version