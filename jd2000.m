function jd2000 = jd2000(iy, im, id, varargin)
% jd2000 = jd2000(iy, im, id, ut)
% jd2000 = jd2000(iy, im, id)
% Input:  Year iy(:) [yyyy]
%        Month im(:) [1-12]
%          Day id(:) [1-31]
%              ut(:) [hours, 0-24] 
% Output: Modified Julian Day 2000 JD2000(:)
%         (0.0 = January 1, 2000, 00:00 UTC)

if nargin == 3; 
    ut = 0;
else
    ut = varargin{1}; 
end
% determine size of input arrays
max_size = max([size(iy); size(im); size(id); size(ut)]); 
max_length = max_size(1)*max_size(2);
% convert to matrix if input parameter is scalar
if length(iy)    == 1; iy = iy*ones(max_size); end;
if length(im)    == 1; im = im*ones(max_size); end;
if length(id)    == 1; id = id*ones(max_size); end;
if length(ut)    == 1; ut = ut*ones(max_size); end;
% check for equal length of all input parameters
if size(ut) ~= size(iy) | size(ut) ~= size(im) | size(ut) ~= size(id);
    error('Variables must be of equal size (or scalars)');
    return
end

jd2000 = datenum(iy, im, id+ut/24) - 730486;
