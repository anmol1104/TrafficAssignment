# TrafficAssignment
improved Traffic Assignment by Paired Alternative Segments (iTAPAS) algorithm for static multi-class traffic assignment problem with generalized link cost function.

**Generalized link cost function: c(v) = Σ ℿₚ * p(v) * t**
-   c(v): generalized link cost for link ij
-   t   : travel time on link ij
-   v   : travel speed on link ij
-   ℿₚ  : cost of parameter p
-   p(v): parameter p as a polynomial function of v

**Required properties of the generalized cost function**
-   Strictly positive
-   Monotonically non-decreasing
-   Continuously differentiable

**Arguments**
-   networkName : network from the repository https://github.com/anmol1104/TrafficAssignment/tree/main/src/Network
-   tol         : tolerance level for relative gap convergence
-   maxIters    : maximum number of iterations
-   maxRunTime  : maximum wall clock run time
-   log         : presents results for every iteration if log is on

**DataFiles (available at: https://github.com/anmol1104/TrafficAssignment/tree/main/src/Network)**
-   cost    : Enlists cost (ℿₚ) for all the parameters (p) of the generalized cost function
-   coef    : Enlists coefficients of p(v) for all the parameters (p) of the generalized cost function
-   class   : Enlists the relevant subset of parameters for the generalized cost function for each class
-   network : Details the topology of the network
-   demand  : Enlists OD pairs and corresponding deman


**IO Units**
    -   length  : miles
    -   time    : hour
    -   volume  : litre (of fuel)
    -   mass    : kg (of emissions)
    -   cost    : $

For more details on the iTAPAS algorithm refer to:
    - Bar-Gera, H., 2010. Traffic assignment by paired alternative segments. Transp. Res. Part B Methodological 44, 1022-1046.
    - Xie, J., Xie, C., 2016. New insights and improvements of using paired alternative segments for traffic assignment. Transp. Res. Part B Methodological 93, 406-424.
