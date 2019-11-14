clc; clear all;
%close all;
%sats = [3];
sats = [3,10,11,22,31];
sats = 26
%sats = 1:32;
%sats = 31;


for sat = sats
    phase_filename = "phase_error_"+sat+".csv";
    if exist(phase_filename, 'file')

        phase_error = csvread("phase_error_"+sat+".csv");
        figure
        hold on
        t = (1:length(phase_error(1:end-1)))/1000-1;
        %t = (1:length(phase_error(1:end-1)))-1;
        plot(t,phase_error(1:end-1));%*180/pi);
        %plotting 90 and -90 degrees bounds
        %plot(t,ones(1,length(phase_error)-1)*90,'r')
        %plot(t,ones(1,length(phase_error)-1)*-90,'r')

        title(['PRN# ',num2str(sat)]);
        legend("Phase Error");
        title(['PRN# ',num2str(sat)]);
%         stdii = [];
%         for ii = 1:length(phase_error)-1000
%             stdii(ii) = norm(phase_error(ii:ii+1000))/sqrt(1000);
%         end
%         figure; plot((1:length(stdii))/1000, stdii, (1:length(stdii))/1000, 0*stdii + 0.35,'--')
    end
end