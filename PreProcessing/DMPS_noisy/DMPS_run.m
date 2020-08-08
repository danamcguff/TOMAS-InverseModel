% This script is written to run size inversion from raw data 110618
clear; clc
% dir_code = '/home/yuanlong/Documents/MATLAB/DPDMA_inv/DMPS/';
dir_code = '';
% !! dir_code stores the directory of the code, note the '/' at the end
DMPS_prep
% DMPS_prep prepares DMA and CPC working conditions for inversion
% Hauke-type DMA-1 (11 cm) and DMA-2 (28 cm), CPC-3025 and CPC-3010

%% start inversion here
load([dir_code 'Dpb.mat']);
% file 'Dpb.mat' stores the diameter of PNSD, row vector
% unit of meter, bin for inverted particle size, same as Jorma's data

PNSD2 = DMPS_inv(Dpb,DMA2,CPC3010,[],1);
% run DMA2 inversion first
PNSD1 = DMPS_inv(Dpb,DMA1,CPC3025,PNSD2,1);
% use PNSD2 as the initial guess for PNSD1 inversion

%% combine PNSD1 and PNSD2
% DMA-1 gives 2.7-32 nm, DMA-2 gives 15-600 nm
% use 2.7-20 nm from DMA-1 and 32-600 nm from DMA-2, interpolate 20-32 nm
Dpb = PNSD2.Dpb;
Dpb_ = [Dpb(1:44) Dpb(55:end)]; % 2.7 - 20 nm & 32 - 600 nm
dNdlgDp_ = [PNSD1.dNdlgDp(:,1:44), PNSD2.dNdlgDp(:,55:end)];
dNdlgDp = spline(Dpb_,dNdlgDp_,Dpb);
dNdlgDp(dNdlgDp<0) = 0;
% use spline to interpolate values betwwen 20-32 nm
dlgDp = PNSD2.dlgDp;
Tscan = PNSD2.Tscan;

PNSD.dNdlgDp = dNdlgDp;
PNSD.Dpb = Dpb;
PNSD.dlgDp = dlgDp;
PNSD.Tscan = Tscan;
% recalculate TNum, TSur, TVol
PNSD.TNum = dlgDp*dNdlgDp'; % total num conc, # cm-3
PNSD.TSur = Dpb.^2.*dlgDp*dNdlgDp'*pi*1e12; % total surf conc, um2 cm-3
PNSD.TVol = Dpb.^3.*dlgDp*dNdlgDp'*pi/6*1e18; % total volume conc, um3 cm-3

clear dNdlgDp dNdlgDp_ Dpb Dpb_ dlgDp Tscan

%% Comparison I: temporal profile of total number concentration
load([dir_code 'JJdata.mat']);
Dpb = JJinv(1,3:end);
dNdlgDp = JJinv(2:end,3:end);
TNum = JJinv(2:end,2);

figure(1)
set(gcf,'position',[150 250 1000 420]);
subplot(1,2,1)
plot(PNSD.Tscan, PNSD.TNum, 'r-', PNSD.Tscan, TNum,'bo')
ylim([0 1e4])
legend('YH','JJ')
ylabel('Number Concentration (# cm^{-3})');
datetick
subplot(1,2,2)
plot(TNum,PNSD.TNum, 'ro')
hold on
plot(1:max(TNum),1:max(TNum),'k-');
hold off
xlabel('JJ'); ylabel('YH')
axis([2e3 1e4 2e3 1e4])
suptitle('Comparison between JJ and YH Inversion')

%% Comparison II: particle number size distribution
figure(2)
Nscn = size(PNSD.dNdlgDp,1);
for i = 1:Nscn
    semilogx(PNSD.Dpb,PNSD.dNdlgDp(i,:),Dpb,dNdlgDp(i,:))
    legend('YH','JJ')
    title(['Scan ' num2str(i)]);
    xlabel('Diameter (nm)')
    ylabel('dN/dlgD_p (# cm^{-3})')
    pause(0.2)
end
clear Nscn i





