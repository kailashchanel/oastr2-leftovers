import Pkg
Pkg.add("Cosmology")
Pkg.add("Measurement")
Pkg. add("Unitful")
using Cosmology
using Measurements
using Unitful

#define cosmological model. For this example I will use the Planck 2015 
#cosmological parameters but this can be easily modified. 

c = cosmology(h=0.7, OmegaM=0.3, OmegaR=0)

age(c, 0)

age(c, 1.42) #z = 1.42


