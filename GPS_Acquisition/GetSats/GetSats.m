clear
clc
close all
addpath('functions')

c = 299792458;

getNewTLE = true; % fetches most recent TLE from celelstrack if true
plotDoppler = true; % plots doppler for each el mask if true
constellations = 5; %if Constellation is 5 go "deep space" if statement of sgp4init to false 
duration = 1*60*60; % seconds
TLEO = 60; % seconds
El_mask_vec = [0, 10, 20, 30]; % degrees
El_mask_vec = 0;
bin_size = 5*60; % seconds

% % start time in UTC
year = 2016;
mon = 6;
day = 22;%18;
hr = 9;%10;
minute = 33;
sec = 06;%25;

rxPos_lla = [33.952790; -117.319165; 415];
Fc = 137.5e6;   % Hz
% bad_IDs = [23545, 23546, 25158, 25413, 25476, 25478, 25479, 25480,...
%     25482, 25983, 25986, 33060, 33061, 33062, 33063, 33064, 33065,...
%     40089, 41181, 41180, 40088, 40090, 21576];
bad_IDs = [];

%%%%%%%%%%%%%%%%%% code begins %%%%%%%%%%%%%%%%%%%%%%%%%
num_bins = ceil(duration/bin_size);
duration = num_bins*bin_size;
npts = ceil(duration/TLEO);
[offset, ~] = mdh2days(year, mon, day, hr, minute, sec);
satrec = GenerateLEOSats(constellations, TLEO, offset, npts, getNewTLE, true);
rxPos = lla2ecef(rxPos_lla);

bin_vec = bin_size/2:bin_size:duration;
tVecLEO = (0:npts-1)*TLEO;
num_ticks = 12;
skip_size = floor(length(tVecLEO)/num_ticks);
xvec = datenum(year,mon,day,hr,minute,sec) + tVecLEO/(24*3600);
xvec_bin = datenum(year,mon,day,hr,minute,sec) + bin_vec/(24*3600);

El_mask_rad_vec = El_mask_vec/180*pi;
num_avail_sats = zeros(length(El_mask_vec), npts);
for ll = 1:length(El_mask_vec)
    El_mask_rad = El_mask_rad_vec(ll);
    legends = {};
    IDs = [];
    for ii = 1:length(satrec)
        rs = satrec(ii).r_s_ecef;
        vs = satrec(ii).v_s_ecef;
        satnum = satrec(ii).satnum;
        [~, El] = GetAzEl(rs(1,:),rs(2,:),rs(3,:),rxPos);
        ind_av = (El>El_mask_rad);
        num_avail_sats(ll, :) = num_avail_sats(ll, :) + ind_av;
        if plotDoppler
            if any(ind_av) && ~any(bad_IDs==satnum)
                rs_av = rs(:, ind_av);
                vs_av = vs(:, ind_av);
                delr = rxPos*ones(1,size(rs_av,2)) - rs_av;
                delv = -vs_av;
                d = sqrt(sum((delr).^2));
                fDhat = NaN(1, npts);
                fDhat(ind_av) = -(sum(delv.*delr)./d)/c*Fc;
                figure(ll+1);
                hold on;
                plot(xvec, fDhat, 'linewidth', 2);
                IDs = [IDs, ii];
                legends{end + 1} = ['SV ID = ', num2str(ii)];
            end
        end
    end
    if plotDoppler
        figure(ll+1);
        title(['Doppler frequency (Hz), Elevation mask = ', num2str(El_mask_vec(ll)),' deg'])
        legend(legends)
        xticks(xvec(1:skip_size:end))
        datetick('x','HH:MMPM','keepticks')
    end
end
num_avail_sats_bin = zeros(length(El_mask_vec), length(bin_vec));
for ll = 1:length(El_mask_vec)
    num_avail_sats_ll = reshape(num_avail_sats(ll, :), size(num_avail_sats,2)/num_bins, num_bins);
    num_avail_sats_bin(ll,:) = mean(num_avail_sats_ll);
end
figure(1)
bar(xvec_bin'*ones(1,length(El_mask_vec)), num_avail_sats_bin','stacked')
xvec_bin_ticks = [2*xvec_bin(1)- 1*xvec_bin(2), xvec_bin, 2*xvec_bin(end)-1*xvec_bin(end-1)];
xticks(xvec_bin_ticks)
datetick('x','HH:MMPM','keepticks')