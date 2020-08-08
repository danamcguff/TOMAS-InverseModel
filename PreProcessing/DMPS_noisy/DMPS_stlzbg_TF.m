function [sTF,Stlzbg] = DMPS_stlzbg_TF(DMA,Vs,ne)
% return theoretical calculation of transfer function and diameter etc. at
% a specific voltage (Vs) and charge number (ne) for one type of DMA
% -------------------------------------------------------------------------
% sTF is a matrix: [Z_tilde,TF,Z,Dp]
% Stlzbg is a struct: beta,delta,G,sigma,Zstar,Dstar,Vs,Dpstar
% -------------------------------------------------------------------------
% created 2019/04/07, yhuang@caltech.edu
%% characterize parameters-------------------%
Qa = DMA.flow(1); % aerosol inlet flow, m3 s-1
Qc = DMA.flow(2); % classified outlet flow, m3 s-1
Qsh = DMA.flow(3); % sheath flow, m3 s-1
Qex = DMA.flow(4); % excess flow, m3 s-1
beta = (Qa+Qc)/(Qsh+Qex);
delta = (Qc-Qa)/(Qc+Qa);

L = DMA.cnfg(1); % m, length of column
r2 = DMA.cnfg(2); % m, outer radius
r1 = DMA.cnfg(3); % m, inner radius, from B&W pg 559
TK = DMA.TK; % K, temperature
Ptorr = DMA.Ptorr; % torr, pressure
gamma = (r1/r2)^2;
kappa = L*r2/(r2^2-r1^2);
Agamma = (-1/2*(1+gamma)*log(gamma)-(1-gamma))^-1;

eC = 1.60217646e-19; % C, electron charge in coulomb
kB = 1.38064852e-23; % m2 kg s-2 K-1, Boltzmann constant

% reference e-mobility
Zstar = (Qsh+Qex)/4/pi/L/Vs*log(r2/r1); % m2 V-1 s-1 (C kg-1 s)

%% calculate sigma
% Non-radial-flow-corrected G, from Stolzenburg 88 Eq2.64 & Apndx B
% These equations assume Poiseuille Flow (not Plug flow).
Ir_w = @(w) Agamma^2 * (1-gamma)^-1 * (...
    - 1/2*w^2*((1-gamma)*log(w)-(1-w)*log(gamma))^2 ...
    + (1/2*w^2*(1-gamma)+1/3*w^3*log(gamma))*...
      ((1-gamma)*log(w)-(1-w)*log(gamma)) ...
    + 1/4 *(1-w^2)*(1-gamma)^2 ...
    + 5/18*(1-w^3)*(1-gamma)*log(gamma) ...
    + 1/12*(1-w^4)*(log(gamma))^2);
G = 4*(1+beta)^2/(1-gamma) * (Ir_w(gamma) + (2*(1+beta)*kappa)^-2);

Dstar = kB*TK/ne/eC*Zstar; % reference diffusivity, m2 s-1
sigma = sqrt(G*Dstar*4*pi*L/(Qsh+Qex));
% Stolzenburg form of sigma^2 = G*D*4*pi*L/(Qsh+Qex)
Stlzbg.beta = beta;
Stlzbg.delta = delta;
Stlzbg.G = G;
Stlzbg.sigma = sigma;
Stlzbg.Zstar = Zstar;
Stlzbg.Dstar = Dstar;
Stlzbg.Vs = Vs;
Stlzbg.Dpstar = DMPS_getDp(Zstar,ne,TK,Ptorr); % m, particle diameter

%% Stolzenburg transfer function
z_tilde = (0.1:0.01:1.9)';
sTF(:,1) = z_tilde; % dimensionless mobility

Int_err = @(x) x.*erf(x)+1/sqrt(pi)*exp(-x.^2);
Omega = @(x) sigma/sqrt(2)/beta/(1-delta) ... % diffusive form
    .*(Int_err((x-(1+beta))/sqrt(2)./sigma) ...
     + Int_err((x-(1-beta))/sqrt(2)./sigma) ...
     - Int_err((x-(1+beta*delta))/sqrt(2)./sigma) ...
     - Int_err((x-(1-beta*delta))/sqrt(2)./sigma));
sTF(:,2) = Omega(z_tilde);
sTF(:,3) = sTF(:,1)*Stlzbg.Zstar;
sTF(:,4) = DMPS_getDp(sTF(:,3),ne,TK,Ptorr); % convert Z to Dp (m)

end
