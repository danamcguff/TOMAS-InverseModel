function Dp = DMPS_getDp(Z_i,charge,T,p)
% This function gives the diameter correspoonding to e-mobility (C kg-1 s)
% with electron number (charge) at specific T (K) and p (mmHg or torr).
% -------------------------------------------------------------------------
% created: 2019/04/3, yhuang@caltech.edu
% -------------------------------------------------------------------------
% This function can handle matrix, vector, or scalar input of Z_i. If input
% is a matrix, should be careful about the size of matrix (m*n), since this
% function transforms matrix to vector (m*n,1) and MATLAB creates a huge
% Jacobian matrix (m*n,m*n), which may eat all the memory of your computer.
% For example: a 201*201 matrix input of Z_i, will create a 201*201*201*201
% Jacobian matrix, which occupies about 12.2 GB memory.
% comment & test: 2019/04/05, YH
% update the guess FUN of Dp to improve fslove output, 2019/04/07, YH
%% calculate viscosity and mean free path of air
mu0 = 1.8324e-5; % reference viscosity, kg m-1 s-1
T0 = 296.15; % reference temperature, K
S = 110.4; % Sutherland constand for dry air, K
mu = mu0*(T/T0)^(3/2)*(T0+S)/(T+S);
% viscosity of air, Schlichting 1968, kg m-1 s-1

fp0 = 6.73e-8; % reference mean free path, m
p0 = 760; % reference pressure, mmHg
fp = fp0*(T/T0)*(p0/p)*(1+S/T0)/(1+S/T);
% mean free path of air molecule, Willeke 1976, m

%%
eC = 1.60217646e-19; % C, electron charge in coulomb
c0 = charge*eC/3/pi/mu;
% constants for calculating the slip correction factor
% numbers from:  from Hutchins et al 1995 AS&T
a = 1.2310; % +/- 0.0022 S.E.
b = 0.4695; % +/- 0.0037 S.E.
c = 1.1783; % +/- 0.0091 S.E.

Zi = Z_i(:);
options = optimset('display','off');
f_Zi = @(x) c0./x.*(1+2*fp./x.*(a+b*exp(-c*x/2/fp)));
Dp = fsolve(@(x) f_Zi(x)-Zi, Dp_guess(Zi), options); % m
if size(Z_i,2) > 1
Dp = reshape(Dp,size(Z_i,2));
end

% ---------------- guess function for Dp ---------------- %
    function Dp_g = Dp_guess(mobility)        
        x = logspace(-9,-5,1e3);
        y = f_Zi(x);
        Dp_g = spline(y,x,mobility);
    end

end

