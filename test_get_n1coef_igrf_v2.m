clearvars
close all

%%

% time vector
year = [2000:2020]';
month = [1:12 1:9]';
day = [1:10 13 15 17 20 21:27]';
hour = [1:21]';

Time_MSDN = datenum(year,month,day,...
        hour,zeros(21,1),zeros(21,1));

% get dipole coefficients version v1
tic
for i=1:length(Time_MSDN)
    [g10_v1(i),g11_v1(i),h11_v1(i)] = get_n1coef_igrf(Time_MSDN(i));
end
toc

% get dipole coefficients version v2
tic
[g10_v2,g11_v2,h11_v2] = get_n1coef_igrf_v2(Time_MSDN);
toc

diff_g10 = g10_v1-g10_v2;
diff_g11 = g11_v1-g11_v2;
diff_h11 = h11_v1-h11_v2;