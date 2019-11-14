clc; clear all;
%close all;
sats = [1,3,10,11,14,22,31];
c = 299792458;

for sat = sats
    f_filename = "debug_vec_"+sat+".csv";
    if exist(f_filename, 'file')    

        debug = csvread("debug_vec_"+sat+".csv");
        
        a_f0 = debug(1);
        a_f1 = debug(2);
        a_f2 = debug(3);
        
        figure
        hold on
        t = linspace(0,10,2e5);
        deltat_sv = (a_f0 + a_f1 * t + a_f2 * t.^2);
        plot(t(1:end),(deltat_sv - 1*deltat_sv(1))*c);
        %debug = debug(23870:end);
        %plot((0:length(debug(1:end-1))-1)*0.001, (debug(1:end-1)-debug(1)));

        legend("Deltat_sv");
        title(['PRN# ',num2str(sat)]);
    end
end
