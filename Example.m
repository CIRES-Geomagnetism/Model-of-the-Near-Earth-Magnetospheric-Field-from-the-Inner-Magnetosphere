%-------------------------------------------------------------------------%
%
% Matlab script to compute predictions of the external, internal and total
% magnetic fields generated in the Earth's inner magnetosphere using 
% hourly time series of external and internal spherical harmonic 
% coefficients derived in Fillion et al. (2023).
% The model currently spans from February 2, 1997 to October 31, 2022 
%
%-------------------------------------------------------------------------%

clearvars 
close all 

%% inputs

% positions in geocentric coordinates and times at which the fields 
% are computed 
Radius = [6400; 6500; 6800]; % Distance from the Earth's center in km
Colatitude = [50; 90; 115]; % in degree
Longitude = [0;100;250]; %in degree
Time_MJD2k = [jd2000(2016,3,12,3.2);...
    jd2000(2017,7,16,8.6);...
    jd2000(2018,11,21,14.8)]; % in Modified Julian Day 2000 (MJD2k). Use jd2000 
% to convert from date to MJD2k    

%% Convert geocentric coordinates to CD coordinates (coefficients are given in CD coordinates)

for i=1:length(Time_MJD2k)
    % compute ith matlab serial date number
    time_MSDN = jd2000_to_serial(Time_MJD2k(i));
    % get dipole coefficients for ith Time_MJD2k
    [g10,g11,h11] = get_n1coef_igrf(time_MSDN);
    [Colatitude_gm(i,1), Longitude_gm(i,1)] ...
        = gg2gm(Colatitude(i,1), Longitude(i,1), 1, [],[],[g10 g11 h11]);
    clear tempi*
end

%% Static parameters (must not be changed)

a = 6371.2; % Earth Radius

lmax = 2; % degree 2 model

if lmax ~= 1 && lmax ~= 2 && lmax ~= 3
    error('Truncation degree must be either 1, 2 or 3')
end

%% Load coefficient time series and coefficient time vector

% Coefficient time series are given in nT

% load coefficients
struct_coef = load(['Coefficients_lmax_' num2str(lmax) ...
    '_06-Feb-1997_04-Dec-2022.mat']);
    
Coef_ext = struct_coef.ext_coef;
Coef_int = struct_coef.int_coef;
Time_coef = struct_coef.Time_MJD2k;

%% Check if query points are within validity period of the model

if ~isempty(find(Time_MJD2k < Time_coef(1),1)) || ...
    ~isempty(find(Time_MJD2k > Time_coef(end),1))
    error('Some query points are outside the validity period of the model')
end

%% Interpolate coefficients time series to query points

% default is linear interpolation
coef_ext_inter = interp1(Time_coef,Coef_ext',Time_MJD2k); 
coef_ext_inter = coef_ext_inter';
coef_int_inter = interp1(Time_coef,Coef_int',Time_MJD2k); 
coef_int_inter = coef_int_inter';

%% Compute external, internal and total magnetic fields (in nT)

% design forward problem matrices for the external field
[A_r_ext, A_theta_ext, A_phi_ext] = ...
    design_SHA_Pol_e_static(Radius./a,...
    Colatitude_gm,Longitude_gm,...
    lmax,lmax);
% design forward problem matrices for the internal field
[A_r_int, A_theta_int, A_phi_int] = ...
    design_SHA_Pol_i_static(Radius./a,...
    Colatitude_gm,Longitude_gm,...
    lmax,lmax);

% compute external field
B_ext_r = sum(A_r_ext.*coef_ext_inter',2);
B_ext_theta = sum(A_theta_ext.*coef_ext_inter',2);
B_ext_phi = sum(A_phi_ext.*coef_ext_inter',2);

% compute internal field
B_int_r = sum(A_r_int.*coef_int_inter',2);
B_int_theta = sum(A_theta_int.*coef_int_inter',2);
B_int_phi = sum(A_phi_int.*coef_int_inter',2);

% compute total field
B_tot_r = B_ext_r + B_int_r;
B_tot_theta = B_ext_theta + B_int_theta;
B_tot_phi = B_ext_phi + B_int_phi;



