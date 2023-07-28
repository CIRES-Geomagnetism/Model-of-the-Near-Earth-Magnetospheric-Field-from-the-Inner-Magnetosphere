function [g10,g11,h11] = get_n1coef_igrf(time)

% Gets first three coefficients of IGRF model at specific epoch (adapted
% from https://github.com/wb-bgs/m_IGRF.
% 
% Usage: [g10,g11,h11] = IGRF(TIME, LATITUDE, LONGITUDE, ALTITUDE, COORD)
% 
% The inputs for
% the position can be scalars or vectors (in the latter case each should
% have the same number of elements or be a scalar), but TIME must be a
% scalar.
% 
% This function relies on having the file igrfcoefs.mat in the MATLAB
% path to function properly when a time is input. If this file cannot be
% found, this function will try to create it by calling GETIGRFCOEFS.
% 
% The IGRF is a spherical harmonic expansion of the Earth's internal
% magnetic field. Currently, the IGRF model is valid between the years 1900
% and 2020. See the health warning for the IGRF model here:
% http://www.ngdc.noaa.gov/IAGA/vmod/igrfhw.html
% 
% Reference:
% International Association of Geomagnetism and Aeronomy, Working Group 
% V-MOD (2010), International Geomagnetic Reference Field: the eleventh
% generation, _Geophys. J. Int._, _183_(3), 1216-1230, 
% doi:10.1111/j.1365-246X.2010.04804.x.
% 
% Inputs:
%   -TIME: Time to get the magnetic field values either in MATLAB serial
%   date number format or a string that can be converted into MATLAB serial
%   date number format using DATENUM with no format specified (see
%   documentation of DATENUM for more information).
% 
% Outputs:
%   - g10,g11,h11: first three gauss coefficients
% 
% See also: LOADIGRFCOEFS, GETIGRFCOEFS, IGRFLINE, DATENUM, IGRF11MAGM.

if nargout ~= 3
    error('Wrong number of output arguments')
end

%%% CHECK INPUT VALIDITY %%%
% Convert time to a datenumber if it is a string.
if ischar(time)
    time = datenum(time);
end

% Make sure time has only one element.
if numel(time) > 1
    error('igrf:timeInputInvalid', ['The input TIME can only have one ' ...
        'element.']);
end

%%% GET PROPER IGRF COEFFICIENTS %%%
[g, h] = loadigrfcoefs(time);

g10 = g(1,1);
g11 = g(1,2);
h11 = h(1,2);


