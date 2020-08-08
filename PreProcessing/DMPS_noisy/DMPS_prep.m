%% specify DMA parameters here
% DMA-1
Qa = 1.5; % lpm
Qc = 1.5; % lpm
% Qa is aerosol inlet flow, Qc is classified outlet flow
Qex = 9.75; % lpm
Qsh = Qex + Qa - Qc; % lpm
% Qex is excess flow, Qsh is sheath flow
L = 0.108; % m, length of column
r1 = 0.025; % m, outer radius
r2 = 0.0335; % m, inner radius, from Birmili et al. AS&T, 1997
DMA1.flow = [Qa; Qc; Qsh; Qex]*1e-3/60; % flow condition, m3 s-1
DMA1.cnfg = [L; r2; r1]; % configuration, m

% DMA-2
Qa = 1.0; % lpm
Qc = 1.0; % lpm
% Qa is aerosol inlet flow, Qc is classified outlet flow
Qex = 6.575; % lpm, a little offset (6.60) to be consistent
Qsh = Qex + Qa - Qc; % lpm
% Qex is excess flow, Qsh is sheath flow
L = 0.28; % m, length of column
r1 = 0.025; % m, outer radius
r2 = 0.0335; % m, inner radius, from Birmili et al. AS&T, 1997
DMA2.flow = [Qa; Qc; Qsh; Qex]*1e-3/60; % flow condition, m3 s-1
DMA2.cnfg = [L; r2; r1]; % configuration, m

DMA1.TK = 297; % K
DMA1.Ptorr = 760; % torr
DMA2.TK = 297; % K
DMA2.Ptorr = 760; % torr

clear Qa Qc Qsh Qex L r2 r1

%% load raw data from .DAT file
fname = 'DM110618.DAT';
fileID = fopen([dir_code fname],'r');
data = textscan(fileID,'%s%s%[^\n\r]');
fclose(fileID);
x = data{1,1}; y = data{1,2};
t1 = zeros(length(x)/42,1);
t2 = zeros(size(t1));
d1 = zeros(14,length(x)/42);
d2 = zeros(22,length(x)/42);
v1 = cellfun(@str2num,x(3:16));
v2 = cellfun(@str2num,x(21:42));
for i = 1:42:length(x)
    t1(floor(i/42)+1) = datenum([x{i}(2:end) y{i}(1:end-1)],...
        'mm-dd-yyyyHH:MM:SS');
    t2(floor(i/42)+1) = datenum([x{i+17}(2:end) y{i+17}(1:end-1)],...
        'mm-dd-yyyyHH:MM:SS');
    d1(:,floor(i/42)+1) = cellfun(@str2num,y(i+2:i+15));
    d2(:,floor(i/42)+1) = cellfun(@str2num,y(i+20:i+41));
end

DMA1.t = t1;
DMA1.v = v1;
DMA1.c = d1;
DMA2.t = t2;
DMA2.v = v2;
DMA2.c = d2;

clear fname fileID data x y t1 t2 d1 d2 v1 v2 i

%% test if the diameter at one voltage is consistent with the setpoint
% % this section can be commented out
% DMA = DMA2;
% Vs = 3700.33;
% ne = 1;
% [sTF,Stlzbg] = DMPS_stlzbg_TF(DMA,Vs,ne);
% figure(1)
% semilogx(sTF(:,4)*1e9,sTF(:,2),'linewidth',2);
% set(gca,'ytick',0:0.2:1)
% xlabel('Diameter nm');
% ylabel('\Omega');
% title(['TF of Dp^* = ' num2str(round(Stlzbg.Dpstar*1e9,2)) ' nm']);
% 
% clear DMA Vs ne

%% CPC detection efficiency
% DMA-1 uses CPC3025 for ultrafine particle counting
CPC3025_curve = ...
    [1.919,  1;
     1.934,  1.6;
     1.977,  2.9;
     2.076,  5.7;
     2.261,  11;
     2.605,  21;
     3.151,  35;
     3.969,  52;
     5.268,  71;
     7.231,  87;
     10.457, 98;
     14.247, 100;
     19.123, 100;
     25.572, 100];
% curve raw data from Fig.7 in Mordas et al. AS&T 2008
% https://doi.org/10.1080/02786820701846252
f0 = @(a,D1,D2,x) (x>=D2*log(a-1)+D1).*(1 - a./(1+exp((x-D1)/D2)));
% EQ.1 from Mertes et al. AS&T 1995
% https://doi.org/10.1080/02786829508965310
g = fittype(f0,'independent','x');
f = fit(CPC3025_curve(:,1),CPC3025_curve(:,2)/100,...
    g,'StartPoint',[1.7,4.3,1.5]);
% fx = logspace(0,2,101);
% fy = f0(f.a,f.D1,f.D2,fx);
% figure(4)
% semilogx(fx,fy,'-',CPC3025_curve(:,1),CPC3025_curve(:,2)/100,'o')
% legend('Fit','Data')
% xlabel('Diameter (nm)'); ylabel('Detection Efficiency')
% title('CPC 3025')
% % plot to check the fitting results, can be commented
CPC3025.a = f.a;
CPC3025.D1 = f.D1;
CPC3025.D2 = f.D2;
CPC3025.flag = 3025;

% DMA-2 uses CPC3010 for normal size counting
CPC3010.a = 1.7;
CPC3010.D1 = 4.3;
CPC3010.D2 = 1.5;
CPC3010.flag = 3010;
% Parameters from Mertes et al. AS&T 1995
% https://doi.org/10.1080/02786829508965310

CPC3772.a = 0.893;
CPC3772.b = 1.480;
CPC3772.d1 = 12.554;
CPC3772.d2 = 2.554;
CPC3772.flag = 3772;

CPC3776.a = 0.958;
CPC3776.b = 10.002;
CPC3776.d1 = 0;
CPC3776.d2 = 1.355;
CPC3776.flag = 3776;
% parameters for CPC 3772 and 3776 from Hermann et al. JAS 2007
% https://doi.org/10.1016/j.jaerosci.2007.05.001

clear CPC3025_curve f0 g f fx fy
