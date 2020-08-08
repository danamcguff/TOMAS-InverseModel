function [Conc,Dpi] = DMA_vDp(DMA)
% return diameter for voltage at each step in ascending order, same as Conc
% applied in DMPS_inv.m
% created 04/08/2019, yhuang@caltech.edu

[Volt,idx] = sort(DMA.v); % stepping voltage
Conc = DMA.c(idx,:); % particle number concentration in measured bin
Nstp = length(Volt);
Dpi = zeros(Nstp,1);
for i = 1:Nstp
    [~,Stlzbg] = DMPS_stlzbg_TF(DMA,Volt(i),1);
    Dpi(i) = Stlzbg.Dpstar; % from Volt to Dpi (m)
end % for i...

end % function DMA_vDp...
