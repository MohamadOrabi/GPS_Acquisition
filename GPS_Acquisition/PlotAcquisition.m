sat = [1,3,10,11,14,22,31];

for i = 1:length(sat)
figure
corrs = csvread("corrs_" + sat(i) + ".csv");
surf(corrs)
title(['PRN#', num2str(sat(i))]);
end