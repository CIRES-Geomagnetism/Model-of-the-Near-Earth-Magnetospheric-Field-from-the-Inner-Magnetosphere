function [g10,g11,h11] = get_n1coef_igrf_v2(time)

% Gets first three coefficients of IGRF model at specific epoch (adapted
% from https://github.com/wb-bgs/m_IGRF. v2 is updated to make calculation
% more efficient for large data volumes.
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
    error('Time must be numeric')
end

% get IGRF coefficient file
if ~exist('igrfcoefs.mat', 'file')
    getigrfcoefs;
end
load igrfcoefs.mat;

% Convert time to fractional years.
timevec = datevec(time);
for i=1:length(time)
    time(i) = timevec(i,1) + (time(i) - datenum([timevec(i,1) 1 1]))./(365 + double(...
        (~mod(timevec(i,1),4) & mod(timevec(i,1),100)) | (~mod(timevec(i,1),400))));
end

% Check validity on time.
yrs = cell2mat({coefs.year});
if ~isempty(find(time < yrs(1),1)) || ~isempty(find(time > yrs(end),1))
    error('igrf:timeOutOfRange', ['This IGRF is only valid between ' ...
        num2str(yrs(1)) ' and ' num2str(yrs(end))]);
end

% Get the nearest epoch that the current time is between.
for i=1:length(time)
    lastepoch(i) = find(yrs - time(i) < 0, 1, 'last');
    if isempty(lastepoch(i))
        lastepoch(i) = 1;
    end
    nextepoch(i) = lastepoch(i) + 1;
end

% Output either g and h matrices or gh vector depending on the number of
% outputs requested.
for i=1:length(time)

    % Get the coefficients based on the epoch.
    tempi_lastg = coefs(lastepoch(i)).g; tempi_lasth = coefs(lastepoch(i)).h;
    tempi_nextg = coefs(nextepoch(i)).g; tempi_nexth = coefs(nextepoch(i)).h;

    % If one of the coefficient matrices is smaller than the other, enlarge
    % the smaller one with 0's.
    if size(tempi_lastg, 1) > size(tempi_nextg, 1)
        tempi_smalln = size(tempi_nextg, 1);
        tempi_nextg = zeros(size(tempi_lastg));
        tempi_nextg(1:tempi_smalln, (0:tempi_smalln)+1) = coefs(nextepoch(i)).g;
        tempi_nexth = zeros(size(tempi_lasth));
        tempi_nexth(1:tempi_smalln, (0:tempi_smalln)+1) = coefs(nextepoch(i)).h;
    elseif size(tempi_lastg, 1) < size(tempi_nextg, 1)
        tempi_smalln = size(tempi_lastg, 1);
        tempi_lastg = zeros(size(tempi_nextg));
        tempi_lastg(1:tempi_smalln, (0:tempi_smalln)+1) = coefs(lastepoch(i)).g;
        tempi_lasth = zeros(size(tempi_nexth));
        tempi_lasth(1:tempi_smalln, (0:tempi_smalln)+1) = coefs(lastepoch(i)).h;
    end

    % Calculate g and h using a linear interpolation between the last and
    % next epoch.
    if coefs(nextepoch(i)).slope
        tempi_gslope_g10 = tempi_nextg(1,1);
        tempi_gslope_g11 = tempi_nextg(1,2);
        tempi_hslope_h11 = tempi_nexth(1,2);
    else
        tempi_gslope_g10 = (tempi_nextg(1,1) - tempi_lastg(1,1))/...
            diff(yrs([lastepoch(i) nextepoch(i)]));
        tempi_gslope_g11 = (tempi_nextg(1,2) - tempi_lastg(1,2))/...
            diff(yrs([lastepoch(i) nextepoch(i)]));
        tempi_hslope_h11 = (tempi_nexth(1,2) - tempi_lasth(1,2))/...
            diff(yrs([lastepoch(i) nextepoch(i)]));
    end
    g10(i) = tempi_lastg(1,1) + tempi_gslope_g10*(time(i) - yrs(lastepoch(i)));
    g11(i) = tempi_lastg(1,2) + tempi_gslope_g11*(time(i) - yrs(lastepoch(i)));
    h11(i) = tempi_lasth(1,2) + tempi_hslope_h11*(time(i) - yrs(lastepoch(i)));
    clear tempi*
end


end


