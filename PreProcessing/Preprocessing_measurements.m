% This is a function to calculate the size distribution parameters used in 
% inverse method with TOMAS code. These parameters are called inventory
% variables here.
% 1) Calculate several values at each time: several size cuts for each (0th - 3rd) moment
% 2) Smooth inventories and full size distribution with a 3-day moving average 
% 3) Interpolate hourly measurements to the model time step ( 5 min )
% 4) Write results to text files set up so FORTRAN code can easily read

function Preprocessing_measurements(Site)

load(sprintf('Measured_%s_20062007.mat',Site) );
% "Measured_%Site%_20062007.mat" : MATLAB workspace containing the
% following data retrieved from http://ebas.nilu.no ...
%%dp_nm : mean diameter of each bin [nm]
%%number : Number concentration in each bin [1/cm3] -- DMPS
%%ndistbn : Number concentration distribution [1/cm3]/log10(dp) -- DMPS
%%time_dmps : datetime array of number measurements
  tnum=time_dmps; 
  
Nz = 16; %number of inventories
Nmoments = 4; %number of moments
Nsizes = Nz/Nmoments; %number of size ranges for each moment
%size ranges for inventory variables
bdy_low = [3, 10, 6, 3];
bdy_high = [6,dp_nm(end), dp_nm(end), dp_nm(end)];

%Define size-cuts with inequality or strict inequality?
Cond_lb = [ {'Ineq'}, {'Ineq'}, {'Strict'}, {'Ineq'}];
Cond_ub = [ {'Ineq'}, {'Ineq'}, {'Ineq'}, {'Ineq'}];

%name of each inventory
inventory_name = [{'N_{3-6} cm^{-3}'},{'Dp_{3-6} \mu m cm^{-3}'}...
    , {'SA_{3-6} \mu m^2 cm^{-3}'}, {'Mdry_{3-6} \mu g m^{-3}'} ...
    , {'N_{10} cm^{-3}'}, {'Dp_{10} \mu m cm^{-3}'}, {'SA_{10} \mu m^2 cm^{-3}'} ...
    , {'Mdry_{10} \mu g m^{-3}'}, {'N_{6} cm^{-3}'}, {'Dp_{6} \mu m cm^{-3}'} ...
    ,  {'SA_{6} \mu m^2 cm^{-3}'}, {'Mdry_{6} \mu g m^{-3}'} ...
    , {'N_{3} cm^{-3}'}, {'Dp_{3} \mu m cm^{-3}'} ...
    ,  {'SA_{3} \mu m^2 cm^{-3}'}, {'Mdry_{3} \mu g m^{-3}'}  ];

dt = 3600; %measurement frequency is 1 hr
deltat = 5*60; %model time-step is 5 min
m = 3; %window of X hours
d = 1; % degree of polynomial to fit

  start_meas = datetime(2006,5,29,0,0,0);
NumDays = 368;

T = length(tnum);
vol_box = 7.5e19; %[cm3/box]
rho = 1310; %Assumed from model [kg/m3]
N_nuc_all = zeros(T,1); N_tot_all = zeros(T,1); SA = zeros(T,1); V_all = zeros(T,1);
N80_all = zeros(T,1); Z = zeros(T,Nz);
Nbins = length( dp_nm );
for t = 1:T
    for k = 1:Nbins
        SA(t) = SA(t) + number(t,k)*100^3*(10^-9*dp_nm(k))^2*pi; % 1/m
        V_all(t) = V_all(t) + number(t,k)*(10^-9*dp_nm(k))^3*pi/6*(1e6)^3; %um3/cm3
        if dp_nm(k) <= 6 && dp_nm(k) >= 3
            N_nuc_all(t) = N_nuc_all(t) + number(t,k);
        end
        if dp_nm(k) >= 10
            N_tot_all(t) = N_tot_all(t) + number(t,k);
        end
        if dp_nm(k) > 80
            N80_all(t) = N80_all(t) + number(t,k);
        end
        for n = 1:Nsizes
            for j = 1:Nmoments
                indx = (n-1)*Nmoments+j;
                if strcmp( Cond_lb{n}, 'Ineq' )
                    if strcmp( Cond_ub{n}, 'Ineq' )
                        if dp_nm(k) >= bdy_low(n) && dp_nm(k) <= bdy_high(n)
                            Z(t, indx) =  Z(t, indx) + ...
                                (dp_nm(k)*1e-3)^(j-1)*number(t,k);
                        end
                    else
                        if dp_nm(k) >= bdy_low(n) && dp_nm(k) < bdy_high(n)
                            Z(t, indx) =  Z(t, indx) + ...
                                (dp_nm(k)*1e-3)^(j-1)*number(t,k);
                        end
                    end
                else
                    if strcmp( Cond_ub{n}, 'Ineq' )
                        if dp_nm(k) > bdy_low(n) && dp_nm(k) <= bdy_high(n)
                            Z(t, indx) =  Z(t, indx) + ...
                                (dp_nm(k)*1e-3)^(j-1)*number(t,k);
                        end
                    else
                        if dp_nm(k) > bdy_low(n) && dp_nm(k) < bdy_high(n)
                            Z(t, indx) =  Z(t, indx) + ...
                                (dp_nm(k)*1e-3)^(j-1)*number(t,k);
                        end
                    end
                end
                
            end
        end
    end
end
%convert 3rd moment to Aerosol Volume
for n=1:Nsizes
    j = Nmoments;
    indx = (n-1)*Nmoments+j;
    Z(:,indx) = pi/6.*Z(:,indx);
end


N_nuc = N_nuc_all;
N_10 = N_tot_all;
Vdry = V_all;
CCN = N80_all;
t_meas = tnum;

save(sprintf('%s_inventories_and_CCN.mat',Site), 'N_nuc', 'N_10', 'Vdry', 'CCN', 'Z', 't_meas');

%% Smooth with Savitzky-Golay filter
T = length( t_meas );
if rem(m,2) == 0
    m = m+1; % horizon must be odd
end
ind_range = round(m*0.75):round(T - m*0.75); % data leftover does not infclude ends
[b,g] = sgolay(d,m); % Coefficients for polynomial
% initialize matrices
rate = zeros(T,Nz);
dx = zeros(T,2, Nz);

setpt_numb = number;
setpt = Z;
%%%%%Get smoothed measurements and their derivatives
for j = 1:Nz
    for p = 0:d
        dx(:,p+1,j) = conv(setpt(:,j), factorial(p)/(-dt)^p * g(:,p+1), 'same');
    end
    setpt_smooth(:,j) = dx(ind_range,1,j);
    rate_smooth(:,j) = dx(ind_range,2,j);
    setpt_raw(:,j) = setpt(ind_range,j);
end

for j = 1:Nbins
    for p = 0:d
        ndist_dx(:,p+1,j) = conv(setpt_numb(:,j), factorial(p)/(-dt)^p * g(:,p+1), 'same');
    end
    ndist_setpt_smooth(:,j) = ndist_dx(ind_range,1,j);
    ndsit_rate_smooth(:,j) = ndist_dx(ind_range,2,j);
    ndist_setpt_raw(:,j) = setpt_numb(ind_range,j);
end



for t = 1:T
    if t == 1
        rate(t,:) = zeros(1,Nz);
    else
        rate(t,:) =  ( setpt(t,:) - setpt(t-1,:) )./dt;
    end
end
rate_raw(:,:) = rate(ind_range,:);
time_smooth = t_meas(ind_range);
Tsm = length(time_smooth);

Z_filt = zeros(Tsm,Nz);
for t = 1:length( ndist_setpt_smooth(:,1) )
    for k = 1:Nbins
        for n = 1:Nsizes
            for j = 1:Nmoments
                indx = (n-1)*Nmoments+j;
                if dp_nm(k) >= bdy_low(n) && dp_nm(k) <= bdy_high(n)
                    Z_filt(t, indx) =  Z_filt(t, indx) + ...
                        (dp_nm(k)*1e-3)^(j-1)*ndist_setpt_smooth(t,k);
                end
                
            end
        end
    end
    if t == Tsm
        rate_fromSmooth(t,:) = (Z_filt(t,:)- Z_filt(t-1,:) )./dt;
    else
        rate_fromSmooth(t,:) = (Z_filt(t+1,:)- Z_filt(t,:) )./dt;
    end
end
%convert 3rd moment to Aerosol Volume
for n=1:Nsizes
    j = Nmoments;
    indx = (n-1)*Nmoments+j;
    Z_filt(:,indx) = pi/6.*Z_filt(:,indx);
    rate_fromSmooth(:,indx) = pi/6.*rate_fromSmooth(:,indx);
end

k1 = datefind(start_meas, time_smooth);
k2 = datefind(start_meas+calendarDuration(0,0,NumDays), time_smooth);


%linearly Interpolate for 5-minutely data
tspan = etime( datevec( time_smooth(k2) ), datevec( time_smooth(k1) ) )/60; %minutes
hrSpan = tspan/60; %hours
show_meas = time_smooth(k1) + calendarDuration(0,0,0,0,0:5:tspan,0);

Ndt = 3600/deltat;
Z_interp = zeros( hrSpan*Ndt+1, Nz ); rate_interp = zeros( hrSpan*Ndt+1, Nz );
Z_interp_filtDerivative = zeros( hrSpan*Ndt+1, Nz );
rate_fromSmooth = zeros( hrSpan, Nz );
for i = 0:hrSpan-1
    Hrs = i*Ndt;
    for j = 1:Ndt
        if j == 1
            Z_interp(i*Ndt+j,:) = Z_filt(k1+i,:);
            Z_interp_filtDerivative(i*Ndt+j,:) = Z_filt(k1+i,:);
        else
            Z_interp(i*Ndt+j,:) = Z_filt(k1+i,:) +( Z_filt(k1+i+1,:)- ...
                Z_filt(k1+i,:) )./Ndt*(j-1);
            Z_interp_filtDerivative(i*Ndt+j,:) = Z_filt(k1+i,:) +rate_smooth(k1+i,:).*((j-1)*deltat);
        end
        rate_interp(i*Ndt+j,:) = ( Z_filt(k1+i+1,:)- Z_filt(k1+i,:) )./dt;
    end
    rate_fromSmooth(i+1,:) = ( Z_filt(k1+i+1,:)- Z_filt(k1+i,:) )./dt;
end
% rate_fromSmooth(end+1,:) = (Z_filt(k2+1,:)- Z_filt(k2,:) )./dt;
Z_interp(hrSpan*Ndt+1,:) = Z_filt(k2,:); 
Z_interp_filtDerivative(hrSpan*Ndt+1,:) = Z_filt(k2,:); 
rate_interp(hrSpan*Ndt+1,:) = rate_smooth(k2,:);
rate_fromSmooth(hrSpan+1,:) = ( Z_filt(k1+i+2,:)- Z_filt(k1+i+1,:) )./dt;
        
nan_indx = find( isnan( sum(Z_filt, 2) )==1 );
%Remove first data point after period of no data
for j = 1:length( nan_indx )
    if j == length( nan_indx)
        if nan_indx(j) < Tsm
            Z_filt( nan_indx(j)+1, :) = NaN;
        end
    elseif j == 1
        if nan_indx(j) ~= 1
            Z_filt( nan_indx(j)-1,:) = NaN;
        end
    elseif nan_indx(j+1)-nan_indx(j) ~= 1
        Z_filt( nan_indx(j)+1, : ) = NaN;
    elseif nan_indx(j)-nan_indx(j-1) ~= 1
        Z_filt( nan_indx(j)-1, : ) = NaN;
    end
end
% figure; plot( show_meas, Z_interp(:,1), '-s', 'Color', 'r', 'MarkerSize',2); hold on;
% % plot( time_smooth(k1:k2), Z_filt(k1:k2,1), 'c*' )
% plot( time_smooth(k1:k2), setpt_raw(k1:k2,1),'*', 'Color', 'k', 'MarkerSize',2 ); 

% k1 = 1; k2 = length( time_smooth );
%%%%%%%%%%%%%%%% P L O T S %%%%%%%%%%%%%%%%%%%%%%%%%%
for j = [1,5,16] %1:Nz
        figure; hold on;
        plot( time_smooth, setpt_raw(:,j),'r*', 'MarkerSize', 2); 
%         plot( show_meas, Z_interp(:,j),'b', 'LineWidth', 1); 
%         plot( show_meas, Z_interp_filtDerivative(:,j),'k--', 'LineWidth', 1); 
        plot( time_smooth, Z_filt(:,j), 'b');
        legend( 'Raw','Smoothed'); 
        ylabel(inventory_name{j} ); hold off
%           figure; plot( setpt_smooth(:,j), Z_filt(:,j), '*', 'MarkerSize',2);
%           xlabel('Filtered Inventory'); ylabel('Inventory from filtered Size Distribution');
%           title( inventory_name{j});
          Error(:,j) = ( abs( setpt_smooth(:,j)-Z_filt(:,j) ))./nanmean( setpt_smooth(:,j))*100;
%         figure; plot( time_smooth, rate_smooth(:,j));
%         ylabel(['derivative of ',inventory_name{j}] );
end

for j = 1:Nbins
    dndlogdp(:,j) = ndist_setpt_smooth(:,j)./log10( dk_nm(j+1)/dk_nm(j) );
end
k = find( time_smooth==datetime(2006,7,1) );
kf = find( time_smooth==datetime(2006,7,4) );
% data = log10(dndlogdp(k:k+240,:) );
data = log10( dndlogdp );
j = find(data==-Inf);
ncol=1;
for i = 1:length(j)
    if j(i) <= length( data(:,1) )
        data( j(i),1 ) = -5;
    else
        m = ceil( j(i)/length(data(:,1) ) );
        if m> ncol
            ncol=m;
        end
        data(j(i)-(ncol-1)*length(data(:,1)),ncol) = -5;
    end
end
figure;
% pcolor(datenum(time_smooth ), log10(dp_nm), data' );
% datetick('x', 'yyyy-mm-dd HH:MM:SS',  'keepticks');
pcolor([1:kf-k+1]./24, log10(dp_nm), log10( dndlogdp(k:kf,:) )' );
xlabel('Days')
lim = ceil( max(max( data )) );
clim = linspace(-lim, lim, lim+1 );
caxis( [clim(1), clim(end)] );
colorbar; clabel = 10.^clim;
colorbar('YTickLabel',clabel)
shading flat; h = gca;
set(h, 'YTick', [0, log10(3), 1,2,3] ); set(h, 'YTickLabel',{'1','3','10','100','1,000'})
ylabel('d_p [nm]'); 
% 
N = 25; %number of lines to plot
   X = linspace(0,pi*3,1000);
   Y = bsxfun(@(x,n)sin(x+2*n*pi/N),X.',1:N);
   C = linspecer(N);
fig = figure; 
axes1 = axes('Parent',fig,'YMinorTick','on','XMinorTick','on',...
    'XScale','log', 'NextPlot', 'replacechildren','ColorOrder',C);
hold on;
for j = 1:25
    plot(dp_nm, 10.^data(k+j-1,:), 'DisplayName', string(time_smooth(k+j-1) ) );
end
legend('show');
xlabel('D_p [nm]');
ylabel('dN/dlog_{10}d_p [cm^{-3}]');

%% Write values to input file!
%%% W R I T E %%% ---> Meas1.txt, Meas2.txt, Meas3.txt
parfor j = 1:Nz
    name1 = sprintf('Meas_inventory%02i.txt', j);
fid = fopen(name1, 'w');

%First, let's write the number of observations for allocating arrays
fprintf(fid, '%12i\n', length( show_meas ) );
fprintf(fid, '%12i\n', Nbins );

X = ndist_setpt_smooth(k1:k2,:);
Y = Z_filt(k1:k2,j); Y_raw = setpt_raw(k1:k2,j);
Ysum = sum( Z_filt(k1:k2,:), 2 );
dYdt = rate_fromSmooth(:,j); %dYdt = rate_interp(:,j);% dYdt = rate_smooth(:,j); %
time_upd = 0;
for x = 1:(k2-k1+1)

    fprintf(fid, '%6i ',  time_upd );
    xk = datefind( time_smooth(x), t_meas );
    if isnan( Ysum(x) )
        fprintf(fid, '%14.12f ', -50000*ones(Nbins,1) );
        Y(x) = -50000;
        Y_raw(x) = -50000;
        dYdt(x) = -50000;
    else
        fprintf(fid, '%14.12f ', X(x,:) );
    end
    if dYdt(x) < 0
         fprintf(fid, '%14.12E %14.12E\n', Y(x), dYdt(x) );
    else
         fprintf(fid, '%14.12E  %14.12E\n', Y(x), dYdt(x) );
    end
    time_upd = time_upd + dt;
end
fclose(fid);

end

fid = fopen('MeasInfo.txt', 'w');
fprintf(fid, '%i ', datevec( show_meas(1) )); 
fprintf(fid, '\n%i\n', Nbins);
fprintf(fid, '%4.3f ', dp_nm);
fprintf(fid, '\n%s\n', Site);
fclose(fid);

fclose('all');
