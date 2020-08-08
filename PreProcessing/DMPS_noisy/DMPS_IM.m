function IM = DMPS_IM(Dpi,Dpb,DMA,CPC,polarity,Lsmpl)
% Dpi: m, the peak position, same length of counting vector
% Dpb: m, particle size in each bin
% DMA: structure, L, r1, r2, Qsh, Qex, Qa, Qm, T, Ptorr
% CPC: parameters for detection efficiency
% polarity: +1/-1, particle charge
% Use Wiedensohler 1988 charging probability distribution
% Lsmpl: m, length of sampling tube
% ----------------------------------------------------------------------- %
% IM: Inversion Matrix
% ----------------------------------------------------------------------- %
% last modified by 2019/04/07,yhuang@caltech.edu
%% constants and parameters
eC = 1.60217646e-19; % C, electron charge in coulomb
kB = 1.38064852e-23; % J/K, Boltzmann constant

%%
TK = DMA.TK; % K, temperature
Ptorr = DMA.Ptorr; % torr, pressure
mu0 = 1.8324e-5; % reference viscosity, kg m-1 s-1
T0 = 296.15; % reference temperature, K
S = 110.4; % Sutherland constand for dry air, K
viscosity = mu0*(TK/T0)^(3/2)*(T0+S)/(TK+S);
% viscosity of air, Schlichting 1968, kg m-1 s-1

fp0 = 6.73e-8; % reference mean free path, m
p0 = 760; % reference pressure, mmHg
mfp = fp0*(TK/T0)*(p0/Ptorr)*(1+S/T0)/(1+S/TK);
% mean free path of air molecule, Willeke 1976, m

% constants for calculating the slip correction factor
% numbers from:  from Hutchins et al 1995 AS&T
gam1 = 1.2310; % +/- 0.0022 S.E.
gam2 = 0.4695; % +/- 0.0037 S.E.
gam3 = 1.1783; % +/- 0.0091 S.E.

chgs = 8; % highest number of charges considered
Nbin = length(Dpb); % number of particle size bin
P_chg = Calc_rough_charge_prob(Dpb,chgs);
% P_chg is the charge probability matrix, Nbin*chg
Zpi = zeros(chgs,Nbin);
for iie = 1:chgs
    Zpi(iie,:) = Calc_Zp_fromDp(Dpb,iie);
end % for iie...

%% for calculation of G and omega
Qa = DMA.flow(1); % aerosol inlet flow, m3 s-1
Qc = DMA.flow(2); % classified outlet flow, m3 s-1
Qsh = DMA.flow(3); % sheath flow, m3 s-1
Qex = DMA.flow(4); % excess flow, m3 s-1
beta = (Qa+Qc)/(Qsh+Qex);  %Baron & Willeke  18-31
delta = (Qc-Qa)/(Qc+Qa);  %Baron & Willeke  18-34

L = DMA.cnfg(1); % m, length of column
r2 = DMA.cnfg(2); % m, outer radius
r1 = DMA.cnfg(3); % m, inner radius, from B&W pg 559
gamma = (r1/r2)^2;
kappa = L*r2/(r2^2-r1^2);
Agamma = (-1/2*(1+gamma)*log(gamma)-(1-gamma))^-1;

G = Calc_G();

%% calculate inversion matrix
Ncnt = length(Dpi); % number of CPC counting bins
Zp_star = Calc_Zp_fromDp(Dpi,1);
IM = zeros(Ncnt,Nbin);
for ic = 1:Ncnt
    IM(ic,:) = Calc_omega(Zp_star(ic));
end % for ic...

%% correct for sampling tube loss
Di = kB*TK*(1+(2*mfp./Dpb).*(gam1+gam2*exp(-gam3*Dpb./(2*mfp))))...
    ./(3*pi*viscosity*Dpb);
eta_t = Tube_loss(Lsmpl,Di,Qa);
IM = IM.*repmat(eta_t,Ncnt,1);

%% correct for CPC3010/CPC3025 counting efficiency
% Calculating the CPC detection efficiency for particles of diameter Dp
% based on the heterogeneous nucleation theory
if CPC.flag == 3010 || CPC.flag == 3025
    a = CPC.a;
    D1 = CPC.D1;
    D2 = CPC.D2;
    f0 = @(x) (x>=D2*log(a-1)+D1).*(1 - a./(1+exp((x-D1)/D2)));
    % x is the diameter in nm, EQ.1 from Mertes et al. AS&T 1995
    % https://doi.org/10.1080/02786829508965310
elseif CPC.flag == 3772 || CPC.flag == 3776
    a = CPC.a;
    b = CPC.b;
    d1 = CPC.d1;
    d2 = CPC.d2;
    f0 = @(x) (x>=d2*log(b/a-1)+d1).*(a - b./(1+exp((x-d1)/d2)));
end
eta_c = f0(Dpb*1e9);

IM = IM.*repmat(eta_c,Ncnt,1);

IM(IM<1e-10) = 0;
%% ************************************************************************
    function G = Calc_G()
        % Non-radial-flow-corrected G, from Stolzenburg 88 Eq2.64 & Apndx B
        % These equations assume Poiseuille Flow (not Plug flow).
        
        G = 4*(1+beta)^2/(1-gamma) * (Calc_Igamma(gamma) + ...
            (2*(1+beta)*kappa)^-2);
        
        function Igamma = Calc_Igamma(w)
            Igamma = Agamma^2 * (1-gamma)^-1 * ...
                (-1/2*w^2*((1-gamma)*log(w)-(1-w)*log(gamma))^2 + ...
                (1/2*w^2*(1-gamma)+1/3*w^3*log(gamma))* ...
                ((1-gamma)*log(w)-(1-w)*log(gamma)) + ...
                1/4*(1-w^2)*(1-gamma)^2 + ...
                5/18*(1-w^3)*(1-gamma)*log(gamma) + ...
                1/12*(1-w^4)*(log(gamma))^2);
        end % function Igamma
        
    end % function Calc_G

%% ************************************************************************
    function omega = Calc_omega(Zpstar)
        % Diffusional transfer function:
        % Omega(Z~) B&W 18-32 / Stolzenburg's thesis 1988
        % Dpm: row vector of diameter in meter,
        % Zp_star: si unit, scalar
        % chg: charge number
        omega = zeros(1,Nbin);
        
        for ie = 1:chgs
            sigma = sqrt(G*4*pi*L*kB*TK*Zpstar/((Qsh+Qex)*eC*ie));
            % calculate sigma: Baron & Willeke 18-35
            
            Zp_tilda = Zpi(ie,:)/Zpstar;
            % Baron & Willeke 18-32
            x1 = (Zp_tilda-(1+beta))./((sqrt(2))*sigma);
            x2 = (Zp_tilda-(1-beta))./((sqrt(2))*sigma);
            x3 = (Zp_tilda-(1+beta*delta))./((sqrt(2))*sigma);
            x4 = (Zp_tilda-(1-beta*delta))./((sqrt(2))*sigma);
            Int_err = @(x) x.*erf(x)+exp(-x.^2)/sqrt(pi);
            omega = omega + sigma/(sqrt(2)*beta*(1-delta)).*...
                (Int_err(x1)+Int_err(x2)-Int_err(x3)-Int_err(x4)).*...
                P_chg(:,ie)';
            % omega has been weighted by charge probability
        end % for ie...
        
    end % function Calc_omega

%% ************************************************************************
    function Zp = Calc_Zp_fromDp(Dpm,chg)
        % Calculate Zp from particle diameter
        % Zp: si unit
        % Dpm: meter, can be a vector
        % chg: unit charge carried by particle, scalar
        Cc = 1 + (2*mfp./Dpm).*(gam1 + gam2*exp(-gam3*Dpm./(2*mfp)));
        Zp = chg*eC*Cc./(3*pi*viscosity*Dpm);
    end % function Calc_Zp_fromDp

%% ************************************************************************
%     function Dpm = Calc_Dp_fromZp(Zp,chg)
%         % Calculate Dp from electrical mobility
%         % Zp: si unit, scalar
%         % Dpm: meter, scalar
%         % chg: unit charge carried by particle, scalar
%         f_Dp = @(x) x - chg*e*...
%             (1+(2*mfp*1e9./x).*(gam1+gam2*exp(-gam3./(2*mfp*1e9./x))))...
%             /(3*pi*viscosity*Zp)*1e9; % 1e9 converts m to nm
%         Dpm = fzero(f_Dp,[1 1e4])/1e9;
%     end % function Calc_Dp_fromZp

%% ************************************************************************
    function charge_prob = Calc_rough_charge_prob(Dp_m,chg)
        % Dp: Particle size in meter
        % chgs:  maximum charges calculated
        %
        % Return a matrix of length(Dp)*chgs
        
        charge_prob=zeros(length(Dp_m),chg); % initialize variable
        % parameter calcs
        Dp_nm = Dp_m*1e9;
        
        ep0 = 8.854e-12; % C^2/(N m^2)
        c1 = eC./(sqrt(4.0*pi*pi*ep0*kB*TK*Dp_m));
        c2 = 2.0*pi*ep0*kB*TK*Dp_m./(eC*eC);
        
        % constants for charging probability calculation
        ratio_of_p_to_n_ion_conc = 1.0;
        ratio_of_p_to_n_ion_mobl = 0.875;
        c3 = log(ratio_of_p_to_n_ion_conc*ratio_of_p_to_n_ion_mobl);
        
        p = [-2.3197,  0.6175,   0.6201, -0.1105, -0.1260, 0.0297; % -1
            -26.3328, 35.9044, -21.4608,  7.0867, -1.3088, 0.1051; % -2
             -2.3484,  0.6044,   0.4800,  0.0013, -0.1553, 0.0320; % +1
            -44.4756, 79.3772, -62.8900, 26.4492, -5.7480, 0.5059];% +2
        % parameters are from Wiedensohler, J Aerosol Sci, 1988
        % -0.1553 (-0.1544 in Wiedensohler 1988) from Baron & Willeke
        o = (0:1:5)';
        f = @(x) p(x,:)*(log10(Dp_nm)).^o;
        power1 = f(2+polarity); % polarity == -1 (-), == 1 (+)
        power2 = f(3+polarity);
        
        % charging calculation
        charge_prob(:,1)=10.0.^power1;
        idp = Dp_nm>=20;
        charge_prob(idp,2) = 10.0.^power2(idp);
        idp = Dp_nm>=70;
        for charge = 3:chg
            charge_prob(idp,charge) = c1(idp).* ...
                exp(-((charge+polarity*c2(idp)*c3).^2)./(2*c2(idp)));
        end % for charge...
    end % function calc_rough_charge_prob

%% ************************************************************************
    function eta = Tube_loss(Leff,Di,Q)
        % calculate the penetration efficiency in a cylindrical tubing
        % only account for diffusion loss (for small particles)
        % Gormley and Kennedy model (Baron and Willeke, EQ 19-20, 19-21).
        
        % Di: diffusion coefficient of particle, m2 s-1
        % Leff: effective length, m
        % Q: Volumetric flow rate, m3 s-1
        % eta: Penetration efficiency
        
        Q_si = Q; % Volumetric flow rate in m3 s-1
        miu = pi*Di*Leff/Q_si; % Dimensionless parameter
        eta = heaviside(miu-0.02).*...
            (0.81905*exp(-3.6568*miu)+0.09753*exp(-22.305*miu)+...
            0.03250*exp(-56.961*miu)+0.01544*exp(-107.62*miu))...
            +(heaviside(miu+1e-3)-heaviside(miu-0.02)).*...
            (1-2.5638*miu.^(2/3)+1.2*miu+0.1767*miu.^(4/3));
        % since heaviside(0) = 0.5, set miu+1e-3 just in case miu = 0
    end % function Tube_loss

end % function Matrix_Generation
