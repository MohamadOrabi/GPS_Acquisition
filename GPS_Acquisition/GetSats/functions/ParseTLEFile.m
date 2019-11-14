function satrec = ParseTLEFile(constellations, url_flag)
% current_dir = '../OrbitCode/TLERetriever/Current/';
% archived_dir = '../OrbitCode/TLERetriever/Archived/';
current_dir = 'data/Current/';
archived_dir = 'data/Archived/';
count = 1;
time_stamp = datestr(datetime);
time_stamp = time_stamp(1:end-3);
time_stamp = strrep(time_stamp,' ','_');
time_stamp = strrep(time_stamp,'-','_');
time_stamp = strrep(time_stamp,':','_');
constellation_names = {'GlobalStar', 'Orbcomm', 'Iridium', 'IridiumNext', 'GPS'};
% 1 Global Star
% 2 Orbcomm
% 3 Iridium
% 4 Iridium Next
% 5 GPS
% 6 Stored TLE File

for nn = 1:length(constellations)
    indx = constellations(nn);
    const_name = constellation_names{indx};
    switch indx
        case 1
            tle_url = 'http://celestrak.com/NORAD/elements/globalstar.txt';
        case 2
            tle_url = 'http://celestrak.com/NORAD/elements/orbcomm.txt';
        case 3
            tle_url = 'http://celestrak.com/NORAD/elements/iridium.txt';
        case 4
            tle_url = 'http://celestrak.com/NORAD/elements/iridium-NEXT.txt';
        case 5
            tle_url = 'https://celestrak.com/NORAD/elements/gps-ops.txt';
        otherwise
            error('Unknown constellation');
    end
    if url_flag
        file_path = [current_dir,const_name,'_',time_stamp];
        websave(file_path, tle_url);
        file_path = [file_path, '.txt'];
    else
        path_listings = dir([archived_dir, const_name,'_*']);
        if isempty(path_listings)
            error([const_name, ' TLE file not available in Archived folder']);
        else
            file_path = [archived_dir, path_listings(1).name];
        end
    end
%     file_path='data/Archived/Orbcomm1_5_Sep_2018_18_10.txt';
    file_path='GPS-22-6-2016.txt';

    fid = fopen(file_path, 'rb');
    numLines = 1;
    while ~feof(fid)
        fgetl(fid);
        numLines = numLines + 1;
    end
    fclose(fid);
    NumSats = floor(numLines/3);
    fid = fopen(file_path, 'rb');
    for ii = 1:NumSats
        L0 = fgetl(fid);
        L1 = fgetl(fid);
        L2 = fgetl(fid);
        if true
            satrecii = twoline2rvMOD(L1, L2);
            satrecii.constellation = indx;
            satrec(count) = satrecii;
            count = count + 1;
        end
    end
    fclose(fid);
end