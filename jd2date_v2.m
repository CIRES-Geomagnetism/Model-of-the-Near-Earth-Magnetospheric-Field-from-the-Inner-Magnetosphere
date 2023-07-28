function varargout = jd2date_v2(mjd2000)
% [iy, im, id, ut] = jd2date(jd2000)
% Input: Modified Julian Day 2000 mjd2000(:)
% Output:  Year iy(:) [yyyy]
%         Month im(:) [1-12]
%           Day id(:) [1-31]
%               ut(:) [hours, 0-24] 

[iy, im,id, ih, min, sec] = datevec(mjd2000(:)+730486);
ut = ih+min/60+sec/3600;
if nargout == 0
    varargout{1} = datestr([iy, im,id, ih, min, sec], 0);
elseif nargout == 1
     varargout{1} = datestr([iy, im,id, ih, min, sec], 0);
else
    varargout{1} = iy;
    varargout{2} = im;
    varargout{3} = id;
    varargout{4} = ut;
end
