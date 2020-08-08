function PNSD = DMPS_inv(Dpbx,DMA,CPC,init,flag)
% dir: directory stores the code
% DMA: basic information about DMA and raw data, struct
% CPC: parameters for detection efficiency
% init: initial guess, can be [] or a PNSD struct
% flag: 1 - run inversion, return PNSD & iIM, 0 - no inversion, only iIM
% PNSD: particle number size distribution, struct
% ----------------------------------------------------------------------- %
% Last modify: 04/08/2019, yhuang@caltech.edu

%% Initialize the inputs
Tscan = DMA.t; % scan time
Nscn = length(Tscan); % total scan number
[Conc,Dpi] = DMA_vDp(DMA); % ascending Dpi (m, vector) and Conc (matrix)

%% calculate intergrated inversion matrix (IM)
polarity = -1;
Lsmpl = 8; % adjusted to make number conc consistent
if isempty(init)
    Dpb = Dpbx;
else
    Dpb = Dpbx(1:54); % this is for DMA1, since IM is small for Dpb > 31 nm
end
% Dpb = Dpbx;
Nbin = length(Dpb); % bin number
IM = DMPS_IM(Dpi,Dpb,DMA,CPC,polarity,Lsmpl);

% use geometric mean to find the new bin boundaries
% and calculate the width of each bin
b0 = sqrt(Dpb(2:end).*Dpb(1:end-1));
b_ = [Dpb(1)^2/b0(1) b0]; % left boundary
b1 = [b0 Dpb(end)^2/b0(end)]; % right boundary
dlgDp = log10(b1./b_); % calculate the delta\log10(Dp)
iIM = IM.*repmat(dlgDp,size(IM,1),1); % integrated form of inversion matrix
iIM(iIM<1e-10) = 0;

PNSD.iIM = iIM;
PNSD.Dpb = Dpb;
PNSD.dlgDp = dlgDp;
PNSD.Tscan = Tscan;

dNdlgDp = zeros(Nbin,Nscn); % to store particle number size distribution
info_e = zeros(2,Nscn); % store inversion information
%% inversion starts from here
if flag
    for i = 1:Nscn
        B_sm = Conc(:,i);
        if isempty(init)
            Guess = 10*interp1(Dpi,B_sm,Dpb,'linear','extrap')';
        else
            B_sm = smoothdata(B_sm,'gaussian','SamplePoints',log(Dpi));
            B_sm(B_sm < 0) = 0;
            Guess = init.dNdlgDp(i,1:Nbin)';
            idx = find(Dpb>=Dpi(1)&Dpb<=Dpi(end));
            Guess1 = 10*interp1(Dpi,B_sm,Dpb(idx),'linear')';
            Guess(idx) = Guess(idx) + Guess1;
        end
        Guess(Guess<0) = 0;
        W_sm = 1./sqrt(B_sm);
        W_sm(isinf(W_sm)) = 0;
        [Sln,n_iter,err] = DMPS_TNNLS(iIM, B_sm, 0.5, 10, W_sm, Guess);
        % Use Total Non-Negative Least Square method (TNNLS)
        Sln(isnan(Sln)) = 0;
        info_e(:,i) = [n_iter; err];
        Sln = smooth(log(Dpb),Sln,'moving');
        Sln(Sln<0) = 0;
        dNdlgDp(:,i) = Sln;
        
        figure(11)
        loglog(Dpb,dNdlgDp(:,i),Dpb,Guess,Dpi,B_sm);
        title(num2str(i));
        axis([1e-9 1e-6 1e-3 1e4])
        legend('Inv','Gus','Raw')
        pause(0.1) % This is used to show the result
    end % for i...
    % hold off
    
    TNum = dlgDp*dNdlgDp; % total num conc, # cm-3
    TSur = Dpb.^2.*dlgDp*dNdlgDp*pi*1e12; % total surf conc, um2 cm-3
    TVol = Dpb.^3.*dlgDp*dNdlgDp*pi/6*1e18; % total volume conc, um3 cm-3
    PNSD.dNdlgDp = dNdlgDp';
    PNSD.TNum = TNum';
    PNSD.TSur = TSur';
    PNSD.TVol = TVol';
    PNSD.info_e = info_e';
end % if flag...

end % function DPDMA_inv
