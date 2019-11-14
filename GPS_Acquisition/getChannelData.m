%Col 1: k (k is incremented every code length (1ms))
%Col 2: Transmission Time (from tracking the satellite)
%Col 3,4,5: Deocded x, y, and z respectively
%Week Number is 144, TOW is 120822 s

clc; clear all;%close all
sats = [1,3,10,11,14,22,31] %[3,14,22,31];
%sats = 26
%sats = [22,31];
c = 299792458;
time = 292;     %in seconds
count_max = 20;
threshold = 0.001;

t = zeros(time*1E3,32);
range = zeros(time*1E3,32);
xs = zeros(time*1E3,32);
ys = zeros(time*1E3,32);
zs = zeros(time*1E3,32);

load('channel.mat')
load('navsol.mat')

index = 1;
max_start = 0;
%Fill ALLData cell

shifts = zeros(1,max(sats)); %Stores the initial shift in the pseudo range for each satellite

for sat = sats 
    if exist("Channel_Data_" + sat + ".csv",'file')
        try
          
            AllData{sat} = csvread("Channel_Data_"+sat+".csv");
            
            Data = AllData{sat};
            start = Data(1,1);            
            t(start:end,sat) = Data(:,1)/1e3;
            range(start:end,sat) = Data(:,2)*c;
            xs(start:end,sat) = Data(:,3);
            ys(start:end,sat) = Data(:,4);
            zs(start:end,sat) = Data(:,5);
            
            %Getting max start
            if start>max_start
                max_start = start;
            end

            %Plot PseudoRange
            figure
            plot(t(start:end,sat),range(start:end,sat));
            title(['PRN# ', num2str(sat)])
            range(start,sat)

            %NavsolPlot
            i = find(channel(14,2439:end) == sat) + 2438;
            tmin(sat) = max([t(start,sat) , channel(2,i(1))]);
            interpolated_t = linspace(max([t(start,sat) , channel(2,i(1))]), min([t(end,sat) , channel(2,i(end))]),300);
            r1 = interp1(t(start:end,sat),range(start:end,sat),interpolated_t);
            r2 = interp1(channel(2,i),channel(8,i),interpolated_t);
            shifts(sat) = r1(1) - r2(1);
            figure
            hold on
            plot(interpolated_t,r1 - r1(1),'r');
            plot(interpolated_t,r2-r2(1));
            title(['PRN# ', num2str(sat), ' ~ With NavSol'])

            figure
            error = (r1-r1(1)) - (r2-r2(1));
            plot(interpolated_t,error);
            title(['PRN# ', num2str(sat), ' ~ Error'])
            %plot(channel(2,i),channel(j,i)- channel(j,i(1))+957.7) 

            %plot(t(start:end,sat),xs(start:end,sat));

%             if Data(1,2)<Data(end,2)
%                 ylim([Data(1,2)-0.1,Data(end,2)])
%             else 
%                 ylim([Data(end,2)-0.1,Data(1,2)])
%             end
            
            
            %sat = sat + 1;

        catch
            warning("Sat " + sat + " did not plot")
        end
    end
end
shifts = [sats; tmin(sats); shifts(sats); shifts(sats)/c];
tmin = [sats;tmin(sats)];

sats_to_use = [1,3,11,14,22,26,31]; %You can add 4
sats_to_use =[1,3,10,11,14,22,31]


%initial guess
ex = 0; %e:estimate
ey = 0;
ez = 0;
cdt = 0;
r = [ex;ey;ez;cdt];

H_k = zeros(length(sats_to_use),4);
h_k = zeros(length(sats_to_use),1);
range_k = zeros(length(sats_to_use),1);
deltax = zeros(4,1);

%output = zeros(time*1e3 - max_start + 1, 3);
output = zeros(10000,3);
%range(:,31) = range(:,31) - 1e-3*c;

if (false)  %Instead of commenting and uncommenting the entire code
close all
for k = max_start:max_start + length(output)%time*1e3
%for k = 240e3-1000:240e3
    
    %Gauss-Newton Iterations
    err = Inf;
    count = 1;
    if k == max_start
        xs_vec_all = [];
        range_all = [];
        for i = 1:length(sats_to_use)
          xs_vec_all =  [xs_vec_all, [xs(k,sats_to_use(i));ys(k,sats_to_use(i));zs(k,sats_to_use(i))]];
          range_all(i, 1) = range(k,sats_to_use(i));
        end
        r = mean(xs_vec_all, 2);
        r = 6371e3*r/norm(r);
        d_all = sqrt(sum((xs_vec_all - r*ones(1, size(xs_vec_all, 2))).^2))';
        dt_all = mean(range_all - d_all);
        r = [r; dt_all];
    end  
    while (err>threshold && count<count_max)
        %Get Jacobian at ex
       
        ex_vec = r(1:3);    
        cdt = r(4);
        for i = 1:length(sats_to_use)
            xs_vec = [xs(k,sats_to_use(i));ys(k,sats_to_use(i));zs(k,sats_to_use(i))];        
            H_k(i,:) = [(ex_vec - xs_vec)'./norm(ex_vec-xs_vec) , 1 ];
            h_k(i,:) = norm(ex_vec - xs_vec) + cdt;
            range_k(i) = range(k,sats_to_use(i));
        end

        deltax = (transpose(H_k)*H_k)\transpose(H_k)*(range_k-h_k);
        r = r + deltax;
        err = norm(deltax);
        count = count + 1;
    end
    clc; count
    output(k-max_start+1,1:3) = r(1:3)';
end
output_lla = ecef2lla(output');
end
%Get Satellite Position LLA
sattt = 22;
satt_pos = [xs(max_start:end,sattt)';ys(max_start:end,sattt)';zs(max_start:end,sattt)'];   %Y has a - sign !!!!
xs_lla = ecef2lla(satt_pos);






