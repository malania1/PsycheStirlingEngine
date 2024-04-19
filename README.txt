***************************************
* Beta Stirling Engine Simscape Model *
***************************************
Contributers: Marc Alania, Mohammed Almozahmi, Muhannad Alsulaiman, 
William “Chase” McPherson, Zoe Wilderspin

** Description **
This model may be used to perform basic simulation of a Beta Stirling
engine. Before the model can be run, the parameters must be initialized.
The parameters can be easily set and initialized using the two parameter
codes included in the file (Small_engine_params.m and Large_engine_params.m).
The parameters will be set to the values in whichever code was run last, so
any time a parameter is changed you must rerun the parameter code to implement
that change. The parameter code allows you to change the engine dimensions, hot 
and cold end temperatures, start up torque, alternator efficiency, and more! 
Feel free to modify the provided parameter codes or even create one of your own,
just make sure the names of each parameter remain the same or else the 
model won't recognize it. 

The model also includes two drive types, a rhombic drive and a flywheel,
and three possible working fluids, hydrogen, helium, and air. To swap out
these components, simply comment out the component in use and comment in the 
component you would like to use. 
