function satrec = GenerateLEOSats(constellations, T_LEO, offset, npts, use_TLE, flag)
% inputs:
% constellations:   array containing the index of constellations
%                   1 Global Star
%                   2 Orbcomm
%                   3 Iridium
%                   4 Iridium Next
%                   5 GPS Operational
% T_LEO:            Sampling time (seconds)
% offset:           Offset from epoch (seconds)
% npts:             Number of points to simulate

% dtr = pi/180;s
% min_per_day = 60*24;

%***************Get TLE Orbital Elements**********************
satrec0 = ParseTLEFile(constellations, use_TLE);
% L = length(satrec0);
% for ii = 1:length(satrec0)
%     satrec0(ii).inclo = (pi/2)*(ii-1)/L;
% end
%satrec0 = satrec0([15,21,54]);
% %****************Get Station Location**************************
% stationfiles=dir([path2data,'\*.mat']);
% if isempty(stationfiles)
%    error('No station files file')
% else
%    load([path2data,'\',stationfiles(1).name]);
%    origin_llh=[lct.station.rx_latitude*dtr;...
%                lct.station.rx_longitude*dtr;...
%                lct.station.rx_altitude];
%    rx_name=lct.station.rx_name;
% end
% 
% %********************Initialize spg4****************************
% %set in sgp4init
% %global tumin mu radiusearthkm xke j2 j3 j4   
% %global opsmode
% rx_name=sscanf(rx_name,'%s');
% fprintf('\n')
% fprintf('Satellite ID %5i \n',satrec.satnum)
% fprintf('Station %s: Lon=%6.4f Lat %6.4f Alt=%6.2f m \n',...
%             rx_name,origin_llh(2)/dtr,origin_llh(1)/dtr,origin_llh(3))

%********************Run SPG4 for Multiple Orbits****************            
%NOTE:  dt, npts, and tsince offset tailored to object 32765 parameters
for ii = 1:length(satrec0)
    epochdaysVec(ii) = satrec0(ii).epochdays;
end
max_epochs = max(epochdaysVec);
for ii = 1:length(satrec0)
    if flag
        tsince_offset(ii) = (offset - epochdaysVec(ii))*1440; %minutes
    else
        tsince_offset(ii) = offset/60 + (max_epochs - epochdaysVec(ii))*1440; %minutes
    end
end

dt = T_LEO/60;  %10 sec
tVec = (0:npts-1)*T_LEO;

for nn = 1:length(satrec0)
%   Expand dt, npts, and tsince_offset for broader selection capabilities
    tsince = tsince_offset(nn) + [0:npts-1]*dt;
    
%     if (satrec0(nn).epochyr < 57)
%         Eyear= satrec0(nn).epochyr + 2000;
%     else
%         Eyear= satrec0(nn).epochyr + 1900;
%     end
%     [Emon, Eday, Ehr, Emin, Esec] = days2mdh(Eyear, 0);
%     UTsec = Ehr*3600 + Emin*60 + Esec;

    xsat_ecf = zeros(3,npts);
    vsat_ecf = zeros(3,npts);
    for k = 1:npts
       [satrec_nn, xsat_ecf(:,k), vsat_ecf(:,k)] = spg4_ecf(satrec0(nn), tsince(k));
    end
    %Scale state vectors to mks units
    xsat_ecf = xsat_ecf*1000;  %m
    vsat_ecf = vsat_ecf*1000;  %mps
    satrec_nn.r_s_ecef = xsat_ecf;
    satrec_nn.v_s_ecef = vsat_ecf;
    satrec_nn.time = tVec;
    satrec(nn) = satrec_nn;
end
end