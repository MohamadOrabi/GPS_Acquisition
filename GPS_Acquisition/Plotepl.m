clc; clear all;
%close all;
sats = [22];

for sat = sats
    ee_filename = "ee_"+sat+".csv";
    pp_filename = "ee_"+sat+".csv";
    ll_filename = "ee_"+sat+".csv";

    if exist(ee_filename, 'file')

        ee = csvread("ee_"+sat+".csv");
        pp = csvread("pp_"+sat+".csv");
        ll = csvread("ll_"+sat+".csv");
        
        figure
        hold on

        plot((1:length(pp)-1)/1e3,ee(1:end-1),'x');
        plot((1:length(pp)-1)/1e3,pp(1:end-1),'o');
        plot((1:length(pp)-1)/1e3,ll(1:end-1),'+');

        legend('early','prompt','late');
        title(['PRN# ',num2str(sat)]);
    end
end

figure
for i = 1:10:length(pp)
    plot([-1 0 1],[ee(i) pp(i) ll(i)])
    ylim([0,3.0028e+08])
    hold off
    title(num2str(i/1e3))
    pause(1e-9)
end
