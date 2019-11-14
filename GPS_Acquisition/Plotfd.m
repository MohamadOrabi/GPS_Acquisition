clc; clear all;
%close all;
sats = [14];
sats = [1,14,22,31];
sats = [3,10,11,22,31];


for sat = sats
    f_filename = "fd_"+sat+".csv";
    if exist(f_filename, 'file')

        fd = csvread("fd_"+sat+".csv");
        figure
        hold on

        plot(fd(1:end-1));

        legend("fd");
        title(['PRN# ',num2str(sat)]);
    end
end
sum(fd(47850:end)/1575.42e6*3e8*0.001)