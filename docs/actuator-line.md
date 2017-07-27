# Actuator line/sector model

The actuator line model (ALM) is implemented by applying body forces to
the flow field. These forces are calculated by sampling the velocity field along
the blades. At each blade point, the lift and drag forces are read from
precomputed tables and applied to the flow field using a Gaussian filter kernel
at a scale \\(\epsilon \\).

## Settings
The ALM in LESGO is implemented in a similar form as the [SOWFA](https://nwtc.nrel.gov/SOWFA) package from NREL, including the input
files. Input files for the NREL 5MW reference turbine are located in "inputATM."
The turbine array is specified in "inputATM/turbineArrayProperties." For example:


    ! Global Properties
    numberOfTurbines        1
    outputInterval          1
    updateInterval          1

    ! First turbine
    TURBINE_1 {
    turbineType             "NREL5MWRef"
    baseLocation            190.0, 315.0, 225.0
    numBladePoints          60
    epsilon                 12.5
    sampling                "atPoint"
    rotationDir             "cw"
    Azimuth                 232.0105
    RotSpeed                9.155
    Pitch                   0
    NacYaw                  0
    fluidDensity            1.23
    numAnnulusSections      1
    annulusSectionAngle     0.
    nacelleFlag             false
    nacelleCd               1.
    TSR                     6.
    tipALMCorrection        .true.
    optimalEpsilon          0.25
    }

Each turbine type is defined in a file that specifies basic details about the
gemoetry of the turbine, the turbine controller, and definitions of the blades.
The lift and drag lookup tables are located in the folder "inputATM/AeroData."

## Output
All output files relating to the turbines can be found in the "turbineOutput"
folder. The following The following quantities are written for each turbine in the files
"turbineOutput/TURBINE_#/...".
* Drag coefficient
* Drag force
* Lift coefficient
* Lift force
* Thrust force
* Tangential force
* Axial force
* Axial velocity
* Relative velocity
* Tangential velocity
* Yaw angle
* Pitch angle
* Rotational speed
* Power

# References
Martínez-Tossas LA, Churchfield M, Meneveau C. "[Large Eddy Simulation of wind turbine wakes: detailed comparisons of two codes focusing on effects of numerics and subgrid modeling](https://doi.org/10.1088/1742-6596/625/1/012024)." *Journal of Physics: Conference Series* **625** (2015). 012024.

Martínez-Tossas LA, Stevens RJAM, Meneveau C. "[Wind farm large-eddy simulations on very coarse grid resolutions using an actuator line model](https://doi.org/10.2514/6.2016-1261)." *34th Wind Energy Symposium* (2016). San Diego, California. AIAA Paper No. 2016-1261.
