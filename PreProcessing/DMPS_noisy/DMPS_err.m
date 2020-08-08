function [Nerr,Err,info, PNSD] = DMPS_err(dNdlgDp,Dpb,DMA,CPC,range)
% dNdlgDp: particle number size distribution, time x bin
% Dpb: size bin of dNdlgDp, row vector, m
% DMA, CPC: information of DMA conditions and CPC detection efficiency
% range: size range to calculate number uncertainty, m
% ----------------------------------------------------------------------- %
% Nerr: uncertainty of size distribution
% Err: uncertainty of number and volume concentration in 'range'
% info: information about the calculation process
% ----------------------------------------------------------------------- %
% created 04/12/2019, yhuang@caltech.edu

%%
[~,Dpi] = DMA_vDp(DMA);
PNSD = DMPS_inv(Dpi',DMA,CPC,[],0);
iIM = PNSD.iIM;
dlgDp = PNSD.dlgDp;

if ~isempty(range)
    idx = find(Dpi>=range(1) & Dpi<=range(2));
else
    idx = [];
end

n_dlgDp = pchip(Dpb,dNdlgDp,Dpi); % interpolated PNSD at Dpi
Nscn = size(n_dlgDp,1);
Nerr = zeros(Nscn,length(Dpi)); % error of dNdlgDp
info = zeros(Nscn,2); % number of iteration, residual error
Err = zeros(Nscn,2); % error of number & volume in specified range
for i = 1:Nscn
    F = n_dlgDp(i,:)';
    Ce_ = 1./sqrt(round(iIM*F*DMA.flow(1)*60000*0.5,0)); % relative error
    Ce_(isinf(Ce_)) = 0;
    Ce = iIM*F.*Ce_; % uncertainty of cpc
    We = 1./sqrt(Ce);
    We(isinf(We)) = 0;
    Guess = 10*Ce;
    [Fe,info(i,1),info(i,2)] = DMPS_TNNLS(iIM, Ce, 0.5, 10, We, Guess);
    Nerr(i,:) = Fe';
    if ~isempty(idx)
        Err(i,1) = sqrt((dlgDp(idx).^2)*(Fe(idx).^2)); % # cm-3
        Err(i,2) = 1e18*pi/6*...
            sqrt((Dpi(idx)'.^6.*dlgDp(idx).^2)*(Fe(idx).^2)); % um3 cm-3
    end
end

end
