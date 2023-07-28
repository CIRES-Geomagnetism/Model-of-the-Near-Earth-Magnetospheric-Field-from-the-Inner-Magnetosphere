# Near-Earth-Magnetospheric-Field-Model-of-the-Inner-Magnetosphere
Contains scripts to compute the external and internal  potential fields using the model of hourly variations of the near-Earth magnetic field generated in the inner magnetosphere and its induced counterpart described in Fillion et al. (2023).  The provided scripts show an example on how the model should be used. 

Inputs
------
Coefficients_lmax_1_06-Feb-1997_04-Dec-2022.mat : Contains hourly time 
series of the model external and induced coefficients used to compute the 
predictions of the degree 1 model between February 6, 1997 and December 4, 
2022

Coefficients_lmax_2_06-Feb-1997_04-Dec-2022.mat : Contains hourly time 
series of the model external and induced coefficients used to compute the 
predictions of the degree 2 model between February 6, 1997 and December 4, 
2022

Coefficients_lmax_3_06-Feb-1997_04-Dec-2022.mat : Contains hourly time 
series of the model external and induced coefficients used to compute the 
predictions of the degree 3 model between February 6, 1997 and December 4, 
2022   

Script
------
Example.m : Show an example of how the model should be used. 
To provide a representation of the fields continous with time, the model 
coefficients are interpolated linearly.


Needs in path/or same directory:
design_SHA_Pol_e_static.m, design_SHA_Pol_i_static.m, gg2gm.m, jd2000.m, 
jd2date_v2.m, jd2000_to_serial.m, jd2000.m

Contact
-------
If you have questions, please contact Martin Fillion at 
martin.fillion@colorado.edu
