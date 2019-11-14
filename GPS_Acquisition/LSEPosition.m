clc; clear all;

sats = [3];

for sat = sats 
    if exist("transmission_time_" + sat + ".csv",'file')
       
        t = csvread("transmission_time_"+sat+".csv");
        x = csvread("x_kvec_"+sat+".csv");
        y = csvread("y_kvec_"+sat+".csv");
        z = csvread("z_kvec_"+sat+".csv");
        
        r = [x,y,z];
        figure
        hold on
        plot(x);
        
        figure
        hold on
        
        plot((0:length(t)-1)*0.001, t);

        legend("Transmission Time");
        title(['PRN# ',num2str(sat)]);
    end

end