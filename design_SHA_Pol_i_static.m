function [A_r, A_theta, A_phi] = design_SHA_Pol_i_static(rho,theta,phi,...
    nmax,mmax)

% [A_r, A_theta, A_phi] = design_SHA_Pol_i_static(rho,theta,phi,...
%    nmax,mmax)
%
% Calculate design matrices A_i that connects the vector of REAL,
% Schmidt-normalized, spherical harmonic expansion coefficients,
%        x = g_{n}^{m} and h_{n}^{m} 
% of a poloïdal internal potential field
% and the magnetic component B_r, B_theta, B_phi:
%        B_r     = A_r*x
%        B_theta = A_theta*x
%        B_phi   = A_phi*x
%
% Inputs:   rho(:)                  radius [units of reference radius]
%           theta(:), phi(:)        co-latitude, longitude [deg]
%           nmax, mmax              maximum degree and order
%
% (Optimized version)
%
% Martin Fillion, 26/07/2019 (from Arnaud Chulliat 2011 script)

rad = pi/180;

N_data = length(theta);
if (length(phi) ~= N_data)
    error('design_SHA_Pol_i_static: theta and phi have different dimensions')
end

[~,dim2] = size(theta);
if dim2 > 1 
    theta = theta'; 
    phi = phi'; 
end;

cos_theta = cosd(theta);
sin_theta = sind(theta);

% convert to row vector if input parameter is scalar

if length(rho) == 1; rho = rho*ones(N_data,1); end;

% number of parameters

N_coeff = mmax*(mmax+2) + (nmax-mmax)*(2*mmax+1);

A_r     = zeros(N_data,N_coeff);
A_theta = zeros(N_data,N_coeff);
A_phi   = zeros(N_data,N_coeff);

k = 0;

for n = 1:nmax

    rn1 = rho.^(-(n+2));

    Pnm = legendre(n, cos_theta, 'sch')';

    dPnm = zeros(size(Pnm,1),n+1);
    dPnm(:,n+1) = (sqrt(n/2).*Pnm(:,n));        % m=n
    dPnm(:,1) = -sqrt(n*(n+1)/2.).*Pnm(:,2);    % m=0
    if n > 1                                    % m=1
        dPnm(:,2) = (sqrt(2*(n+1)*n).*Pnm(:,1) - ...
            sqrt((n+2)*(n-1)).*Pnm(:,3))/2;
    end           
    for m = 2:n-1                               % m=2...n-1
        dPnm(:,m+1) = (sqrt((n+m)*(n-m+1)).*Pnm(:,m) - ...
            sqrt((n+m+1)*(n-m)).*Pnm(:,m+2))/2;
    end
    if n == 1
        dPnm(:,2) = sqrt(2)*dPnm(:,2);
    end   

    for m = 0:min(n,mmax)

        alpha = m*phi*rad;
        cos_phi = cos(alpha);
        sin_phi = sin(alpha);
        
        if m==0
            k = k+1;
            A_r(:,k)     = (n+1).*rn1(:).*Pnm(:,1);
            A_theta(:,k) = -rn1(:).*dPnm(:,1);
            A_phi(:,k)   = -rn1(:)*0;  
        else
            k = k+1;
            A_r(:,k)     = (n+1).*rn1(:).*Pnm(:,m+1).*cos_phi(:);
            A_theta(:,k) = -rn1(:).*dPnm(:,m+1).*cos_phi(:);
            A_phi(:,k)   = ...
                rn1(:).*Pnm(:,m+1)./sin_theta(:).*(m).*sin_phi;
        
            k = k+1;
            A_r(:,k)     = (n+1).*rn1(:).*Pnm(:,m+1).*sin_phi(:);
            A_theta(:,k) = -rn1(:).*dPnm(:,m+1).*sin_phi(:);
            A_phi(:,k)   = ...
                -rn1(:).*Pnm(:,m+1)./sin_theta(:).*m.*cos_phi(:);
        end

    end
end  










