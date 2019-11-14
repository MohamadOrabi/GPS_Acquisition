%sats 10 and 4 are bad (We lose 4 at some point)
clc; clear all;
%close all;
sats = [31];
start = 1234
IQbits = 0
j=1;



for sat = sats
    I_filename = "I_"+sat+".csv";
    Q_filename = "Q_"+sat+".csv";
    if exist(I_filename, 'file')
        I = csvread("I_"+sat+".csv");
        Q = csvread("Q_"+sat+".csv");
        
%         %Decoding
%         for i = start:20:length(I)-20
%            IQbits(j) = abs(sum(I(i:i+19)))/sum(I(i:i+19)); 
%            j=j+1;
%         end
%         IQbits = (-IQbits+1)/2;
%         preamble = [1 0 0 0 1 0 1 1]*-1 +1;
%         index = strfind(IQbits,preamble)
%         diff(index)
%      
%         figure
%         stairs(IQbits)
%         ylim([0,1.5]); grid;
        
        figure
        hold on
        t = (1:length(I)-1)/1e3;
        plot(t,I(1:end-1));
        plot(t,Q(1:end-1));

        legend("I","Q");
        title(['PRN# ',num2str(sat)]);
    end
end