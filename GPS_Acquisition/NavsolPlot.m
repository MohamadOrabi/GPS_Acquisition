%plotting Navsol and Channel 
%2439

clc;

load('channel.mat')
load('navsol.mat')

i = find(channel(14,2439:end) == 3) + 2438;

for j = [7 8]%1 : size(channel,1)-1
    figure
    %plot(1:length(i),channel(j,i))
    plot(channel(2,i),channel(j,i)- channel(j,i(1))+957.7) 
    title(['Row: ', num2str(j)])
end

