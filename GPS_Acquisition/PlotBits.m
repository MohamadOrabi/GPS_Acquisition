clc; clear all;
%close all;
sats = [10];

for sat = sats
    bits_filename = "decoded_bits_"+sat+".csv";
    if exist(bits_filename, 'file')

        bits = csvread("decoded_bits_"+sat+".csv");
        figure
        hold on
        
        prn_in_bits = [1 0 0 0 0 0]
        %preamble = [1 0 0 0 1 0 1 1];
        preamble = [1 0 0 0 1 0 1 1]*-1 +1;
        sat
        pindex = strfind(bits',preamble);
        pindex = unique([pindex , strfind(bits',1-preamble)])


        stairs(bits(1:end-1));

        legend("Bits");
        title(['PRN# ',num2str(sat)]);
    end
end

start = 495;
ids = 0;
j = 1;
for i = start:300:start+300*50
    
TLM_HOW = bits(i:i+30+21)'*-1+1;
frameID = TLM_HOW(end-2:end);
p = TLM_HOW(1:8);
if p*-1+1 == preamble
    TLM_HOW = TLM_HOW*-1+1;
end

%TLM_HOW = bits(i-30-21:i)'*-1+1;
%frameID = TLM_HOW(1:3);
%ids((i-start)/300+1) = binaryVectorToDecimal(frameID);
ids(j) = binaryVectorToDecimal(frameID);
j = j+1;
end
ids

