# TrafficAssignment
improved Traffic Assignment by Paired Alternative Segments (iTAPAS) algorithm for static multi-class traffic assignment problem with generalized link cost function.

**Generalized link cost function: cᵐᵢⱼ = fᵐ(vᵢⱼ) * tᵢⱼ**
-   cᵐᵢⱼ : generalized link cost for link ij, vehicle class m
-   tᵢⱼ  : travel time on link ij
-   vᵢⱼ  : travel speed on link ij
-   fᵐ(v): a polynomial rate of consumption function on v for vehicle class m
           (analogous to rate of fuel consumption or rate of pollutant emission)
           
**Required properties of the generalized cost function**
-   Strictly positive
-   Monotonically non-decreasing
-   Continuously differentiable

**Arguments**
-   networkName : network
-   tol         : tolerance level for relative gap convergence
-   maxIters    : maximum number of iterations
-   maxRunTime  : maximum wall clock run time
-   log         : presents results for every iteration if log is on

**DataFiles** (available at: https://github.com/anmol1104/TrafficAssignment/tree/main/src/Network)
-   class   : Enlists coefficients of fᵐ(v) for each class
-   network : Details the topology of the network
-   demand  : Enlists OD pairs and corresponding demand for each class in passenger car equivalent (PCE)

**IO Units**
-   length  : miles
-   time    : hour
-   volume  : litre (of fuel)
-   mass    : kg (of emissions)
-   cost    : $

For more details on the iTAPAS algorithm refer to:
-   Bar-Gera, H., 2010. Traffic assignment by paired alternative segments. Transp. Res. Part B Methodological 44, 1022-1046.
-   Xie, J., Xie, C., 2016. New insights and improvements of using paired alternative segments for traffic assignment. Transp. Res. Part B Methodological 93, 406-424.
