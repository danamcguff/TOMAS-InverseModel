%% This script returns the error of particle number in a specific range

clear; clc
close all

scan_min = 10; % minutes it takes to measure full size distribution
S = rng;
N= 27; % number of scenarios to loop over

for i_scenario=1:N
    
% load the particle number size distribution data
dir_code = '';
% remember to change this directory
load( [dir_code sprintf('DLMdata%02i.mat', i_scenario)] ); % this file contains the inverted PNSD
Tscan = DLMinv(3:end,1);
Dpb = DLMinv(1,3:end-1); %meters
dNdlgDp = DLMinv(3:end,3:end-1); % # cm-3
TNum = DLMinv(3:end,2); % total number conc, # cm-3

%model and measurement time series
t_inv = Inventories.time_sec./(3600*24); %time points for inventory vals (days)
deltat = t_inv(2)-t_inv(1); %Model time step
tq = Tscan(1):deltat:Tscan(end);
dt = (Tscan(2)-Tscan(1)); %Assume we have measurement at each hour
Nscan = ( dt*60*24 )/scan_min;

% use geometric mean to find the new bin boundaries
% and calculate the width of each bin
b0 = sqrt(Dpb(2:end).*Dpb(1:end-1));
b_ = [Dpb(1)^2/b0(1) b0]; % left boundary
b1 = [b0 Dpb(end)^2/b0(end)]; % right boundary
dlgDp = log10(b1./b_); % calculate the delta\log10(Dp)
TVol = (Dpb.^3.*dlgDp*dNdlgDp'*pi/6*1e18)'; % total volume conc, um3 cm-3
dN = dNdlgDp*diag(dlgDp); 
Ntot = sum( dN, 2 );

%% run this section to retrieve uncertainty
DMPS_prep % run this first to prepare DMA and CPC info

[E_dist36,E_N36,Info36, PNSD36] = DMPS_err(dNdlgDp,Dpb,DMA1,CPC3025,...
    [3 6]*1e-9);
E_N36( isnan(E_N36) ) = 0;
% check DMPS_err.m for more information about inputs and outputs
% E_N36(:,1) is the number uncertainty between 3-6 nm
% E_N36(:,2) is the volume uncertainty between 3-6 nm

[E_dist10,E_N10,Info10, PNSD10] = DMPS_err(dNdlgDp,Dpb,DMA2,CPC3010,...
    [1e-8 Dpb(end)]);
E_N10( isnan(E_N10) ) = 0;

% E_N10(:,1) is the number uncertainty >10 nm
% E_N10(:,2) is the volume uncertainty >10 nm, this can be used as the
% volume uncertainty over the entire range, since the volume of particles <
% 10 nm contributes < 1e-6

Dp_combined = [PNSD36.Dpb PNSD10.Dpb];
E_combined = [E_dist36 E_dist10]./sqrt(Nscan);
E_PNSD = interp1(Dp_combined,E_combined',Dpb)';

w_pnsd = randn( size(E_PNSD) ).*E_PNSD;
PNSD_noisy = dNdlgDp + w_pnsd;
dN_noisy = PNSD_noisy*diag(dlgDp); 
TVol_noisy = (Dpb.^3.*dlgDp*PNSD_noisy'*pi/6*1e18)'; % total volume conc, um3 cm-3
w_vol = 1e18*pi/6*sqrt((Dpb.^6.*dlgDp.^2)*(w_pnsd'.^2))';
idx = find(b_>=10e-9 );
w_n10 = sqrt((dlgDp(idx).^2)*(w_pnsd(:, idx)'.^2))';
TN10_noisy = sum( dN_noisy(:,idx), 2 );
idx = find(b_>=3e-9 & b1<=6e-9);
TN36_noisy = sum( dN_noisy(:,idx), 2 );
w_n36 = sqrt((dlgDp(idx).^2)*(w_pnsd(:, idx)'.^2))';


% kstart = find( Tscan == 0.5 );
% fig = figure;
% Nc = length(Tscan(kstart:end) ); C = jet(Nc);
% axes1 = axes('Parent',fig,'YMinorTick','on','XMinorTick','on',...
%     'XScale','log', 'NextPlot', 'replacechildren','ColorOrder',C);
% hold('on');
% for i = kstart:length( Tscan )
%     x = [ Dpb, fliplr(Dpb) ];
%     y = [dNdlgDp(i,:)+E_PNSD(i,:), fliplr( dNdlgDp(i,:)-E_PNSD(i,:) ) ];
%     patch = fill( x, y, C(i-kstart+1,:) );
%     set(patch, 'edgecolor', 'none');
%     set(patch, 'FaceAlpha', 0.3);
% end
% colorbar;
% colormap('jet');
% caxis([Tscan(kstart), Tscan(end)]);

%% true inventory variables
idx = find(b_>=10e-9 );
TN10 = sum( dN(:,idx), 2 );

idx = 2:4; % bins for 3 nm and 6 nm diameter
TN36 = sum( dN(:,idx), 2 );


%% Noisy measurements

Z = interp1(t_inv, [Inventories.N36.val, Inventories.N10.val, Inventories.Vol.val], ...
     Tscan );
 rng(S);

%Add white noise to true inventory variable
 % st. dev. of noise based on uncertainty calculated above for each scan
 % divided by sqrt( N ) where N= number of scans averaged for hourly
 % measurements
w(:,1) = randn( length(E_N36(:,1)),1).*real( E_N36(:,1)./sqrt(Nscan) )+Z(:,1);
w(:,2) = randn( length(E_N10(:,1)),1).*real( E_N10(:,1)./sqrt(Nscan) )+Z(:,2);
w(:,3) = randn( length(E_N10(:,1)),1 ).*real( E_N10(:,2)./sqrt(Nscan) )+Z(:,3);

Z_noisyMeas = w;
Z_noisyMeas( Z_noisyMeas<0) = 1e-4; %do not allow noisy value to go negative
k2 = find( Tscan==7 ); k1 = find( Tscan==4 );
w_normz(i_scenario,:) = mean( abs(w(k1:k2,:)./Z(k1:k2,:) ) );
Z_std_normz(i_scenario,:) = std( Z_noisyMeas(k1:k2,:) )./mean(Z_noisyMeas(k1:k2,:) );
Z_std_normz_true(i_scenario,:) = std( Z_noisyMeas(k1:k2,:) )./mean(Z(k1:k2,:) );

%% Smooth noisy measurements to mimic the procedure we would take with real data
dt_hrs = dt*24;
Nbins = length( Dpb );
d = 1; m = round(11/dt_hrs); % 11-hr window
if rem(m,2) == 0
    m = m+1; % horizon must be odd
end
[b,g]= sgolay(d,m); clear setpt
Ts = length(Tscan);
ind_range = round(m*0.75):round(Ts - m*0.75); % data leftover does not include ends
dx = zeros(Ts,2, 3); setpt = Z_noisyMeas;
dx_sep = zeros(Ts,2,3); setpt_sep = [TN36_noisy, TN10_noisy, TVol_noisy];
for j = 1:3
    for p = 0:d
        dx(:,p+1,j) = conv(setpt(:,j), factorial(p)/(-dt_hrs*3600)^p * g(:,p+1), 'same');
        dx_sep(:,p+1,j) = conv(setpt_sep(:,j), factorial(p)/(-dt_hrs*3600)^p * g(:,p+1), 'same');
    end
    setpt_smooth(:,j) = dx(ind_range,1,j);
    setpt_smooth_sep(:,j) = dx_sep(ind_range,1,j);
    rate_smooth(:,j) = dx(ind_range,2,j);
    setpt_raw(:,j) = setpt(ind_range,j);
end
Tscan_s = Tscan(ind_range);
tq = Tscan_s(1):deltat:Tscan_s(end);
Z_noisyInterp_nofilt = interp1(Tscan_s, setpt_raw, tq );
Z_noisyInterp = interp1(Tscan_s, setpt_smooth, tq );
Z_noisyInterp_sep = interp1(Tscan_s, setpt_smooth_sep, tq );

k1_s = find( tq==4); k2_s = find( tq==7); 
Z_stdFilt_normz_true(i_scenario,:) = std( Z_noisyInterp(k1_s:k2_s,:) )./mean(Z(k1:k2,:) );
Z_stdFilt_normz(i_scenario,:) = std( Z_noisyInterp(k1_s:k2_s,:) )./mean(Z_noisyInterp(k1_s:k2_s,:) );
Z_avg(i_scenario,:) = mean( Z_noisyInterp(k1_s:k2_s,:) );
Z_stdFilt(i_scenario,:) = std( Z_noisyInterp(k1_s:k2_s,:) );
 N36_std(i_scenario) =  mean(E_N36(k1:k2,1) );
N36_var(i_scenario) =  mean(E_N36(k1:k2,1).^2 ); N36_true(i_scenario) = mean(Z(k1:k2,1) );

%plot results from inventory variables
figure; yyaxis right
hold on
plot( tq, Z_noisyInterp_nofilt(:,3), '-.' );
plot( tq, Z_noisyInterp(:,3), ':','LineWidth', 1.5);
% plot( Tscan, Z(:,3)+w_vol, ':','LineWidth', 1.5);
plot( t_inv, Inventories.Vol.val );
ylabel('Volume Concentration \mum^3 cm^{-3}')
 
yyaxis left
hold on
C = [0.49, 0.18, 0.56]; %RGB for purple
plot( tq, Z_noisyInterp_nofilt(:,1), '-.' );
plot( tq, Z_noisyInterp(:,1),':','LineWidth', 1.5 );
% plot( Tscan, Z(:,1)+w_n36,':','LineWidth', 1.5 );
plot( t_inv, Inventories.N36.val);
plot( tq, Z_noisyInterp_nofilt(:,2), '-.', 'Color', C );
% plot( Tscan, Z(:,2)+w_n10, '-.', 'Color', C );
plot( tq, Z_noisyInterp(:,2), ':','Color',C, 'LineWidth', 1.5);
plot( t_inv, Inventories.N10.val, '-','Color',C ); 
% set( gca, 'YScale', 'log');
ylabel('Number Concentration cm^{-3}')
% legend('noisy N_{3-6}','N_{3-6}',  'noisy N_{10}','N_{10}', ...
%     'noisy Vol', 'Vol');
legend('noisy N_{3-6}','Filtered noisy N_{3-6}','N_{3-6}',  'noisy N_{10}',...
    'Filtered noisy N_{10}','N_{10}', 'noisy Vol', ...
    'Filtered noisy Vol', 'Vol');
xlabel('Days')
title('Inventory time-series, from full simulated size distribution')


T = length(tq);
time_smooth = tq.*(24*3600); %seconds

rate_smooth =diff( [setpt_smooth(1,:)', setpt_smooth']' )./(dt_hrs*3600);
rate_smoothInterp = interp1( Tscan_s, rate_smooth, tq ,'previous');

%% Create text-file in FORTRAN77 format for time, observation, and rate of each measurement
%%time is in SECONDS

%%% W R I T E %%% ---> Meas1.txt, Meas2.txt, Meas3.txt
x = 0;
% % Measurement 1  : setpt_nuc
name1 = sprintf('Meas%03i_inventory01.txt', i_scenario);
fid = fopen(name1, 'w');
%First, let's write the number of observations for allocating arrays
fprintf(fid, '%12i\n', T );
% fprintf(fid, '%12i\n', N);

while x < T
    x = x+1;
    fprintf(fid, '%6i ', round( time_smooth(x) ) );
%     fprintf(fid, '%14.12f ', number(x,:) );
    if rate_smoothInterp(x,1) < 0
         fprintf(fid, '%14.12E %14.12E\n', Z_noisyInterp(x,1), rate_smoothInterp(x,1) );
    else
         fprintf(fid, '%14.12E  %14.12E\n', Z_noisyInterp(x,1), rate_smoothInterp(x,1) );
    end
end
fclose(fid);

x = 0;
% % Measurement 2 : setpt_emis
name1 = sprintf('Meas%03i_inventory02.txt', i_scenario);
fid = fopen(name1, 'w');
%First, let's write the number of observations for allocating arrays
fprintf(fid, '%12i\n', T );
% fprintf(fid, '%12i\n', N );

while x < T
    x = x+1;
    fprintf(fid, '%6i ', round( time_smooth(x) ) );
%     fprintf(fid, '%14.12f ', number(x,:) );
    if rate_smoothInterp(x,1) < 0
         fprintf(fid, '%14.12E %14.12E\n', Z_noisyInterp(x,2), rate_smoothInterp(x,2) );
    else
         fprintf(fid, '%14.12E  %14.12E\n', Z_noisyInterp(x,2), rate_smoothInterp(x,2) );
    end
end
fclose(fid);

x = 0;
% % Measurement 3 : setpt_soa
name1 = sprintf('Meas%03i_inventory03.txt', i_scenario);
fid = fopen(name1, 'w');
%First, let's write the number of observations for allocating arrays
fprintf(fid, '%12i\n', T );
% fprintf(fid, '%12i\n', N );

while x < T
    x = x+1;
    fprintf(fid, '%6i ', round( time_smooth(x) ) );
%     fprintf(fid, '%14.12f ', number(x,:) );
    if rate_smoothInterp(x,1) < 0
        fprintf(fid, '%14.12E %14.12E\n', Z_noisyInterp(x,3), rate_smoothInterp(x,3) );
    else
         fprintf(fid, '%14.12E  %14.12E\n', Z_noisyInterp(x,3), rate_smoothInterp(x,3) );
    end
end
fclose(fid);


end
