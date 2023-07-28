function Time_MSDN = jd2000_to_serial(Time_MJD2k)
% Time_serial = jd2000_to_serial(Time_MJD2k)
%
% Function that converts MJD2k vectors to matlab serial date number
%
% INPUTS:
%       - Time_MJD2k: time vector in Modified Julian day
% 
% OUTPUT:
%       - Time_MSDN: Matlab serial date number
%
% Martin Fillion, 24 Jan 2022

[year,month,day,UT] = ...
    jd2date_v2(Time_MJD2k);
hour = floor(UT);
minute = 60*(UT-hour);
Time_MSDN = datenum(year,month,day,...
    hour,minute,0);

end

