# Wind turbine modeling
LESGO's flexibility to simulate atmospheric boundary layers, explore 
various subgrid scale models, and specify different inflow conditions, is 
complemented by the wind turbine modeling features. LESGO includes two wind 
turbine models:

* The **[actuator disk model](actuator-disk.html)** without rotation is well suited 
for large arrays with coarse grid resolution. It can be activated by setting 
the flag "USE_TURBINES" to true in "CMakeLists.txt".

*  The **[actuator line/sector model](actuator-line.html)** provides a more accurate 
turbine model. While the actuator line/sector model can also be used for large 
arrays, finer grid resolutions are needed compared to the actuator disk model. 
To activate the actuator line/sector model, set the flag "USE_ATM" to true in 
"CMakeLists.txt".
