clc; clear all;%close all
sats = [1 2 3 4 5];

%Rows: Time Cols: Sat#

R_sat = 24e6; %Km
R_earth = 6e6; %Km

bias = 5e6;


points = 100;
xs = zeros(points,length(sats));
ys = zeros(points,length(sats));
zs = zeros(points,length(sats));
bs = zeros(points,length(sats));
rs = zeros(points,length(sats));

degTorad = pi/180;
[xr,yr,zr] = sph2cart(5*degTorad,5*degTorad,R_earth);
receiver_position = [xr,yr,zr];

initial_azimuth = [-20,20,30,-20,5];
initial_elevation = [20,20,-20,-30,5];

initial_azimuth = initial_azimuth.*degTorad;
initial_elevation = initial_elevation.*degTorad;


%initialize sat positions
figure
hold on
[x y z] = sphere;
surf(x*R_earth,y*R_earth,z*R_earth)
scatter3(xr,yr,zr,'x')
for i = 1:length(sats) 
    [x,y,z] = sph2cart(initial_azimuth(i), initial_elevation(i),R_sat);
    xs(1,i) = x;
    ys(1,i) = y;
    zs(1,i) = z;
    rs(1,i) = norm([x-xr,y-yr,z-zr]);
    
    vel = rand(1,3) * 20e3;
    scatter3(x,y,z,'*')
    for t = 2:points
        xs(t,i) = xs(t-1,i) + vel(1);
        ys(t,i) = ys(t-1,i) + vel(2);
        zs(t,i) = zs(t-1,i) + vel(3);
        %scatter3(xs(t,i),ys(t,i),zs(t,i),'*')
        rs(t,i) = 3 * randn + norm([xs(t,i)-xr,ys(t,i)-yr,zs(t,i)-zr]) + 1*bias;
    end
    plot3(xs(:,i),ys(:,i),zs(:,i))
end

threshold = 1e-4;
count_max = 100;




for k = 1:points   
    %Gauss-Newton Iterations
    if k == 1
        xs_vec_all = [];
        range_all = [];
        for i = 1:length(sats)
          xs_vec_all =  [xs_vec_all, [xs(k,i);ys(k,i);zs(k,i)]];
          range_all(i, 1) = rs(k,i);
        end
        r = mean(xs_vec_all, 2);
        r = R_earth*r/norm(r);
        d_all = sqrt(sum((xs_vec_all - r*ones(1, size(xs_vec_all, 2))).^2))';
        dt_all = mean(range_all - d_all);
        r = [r; dt_all];
    end  
    err = Inf;
    count = 1;
    while (err>threshold && count<count_max)
        %Get Jacobian at ex
        ex_vec = r(1:3);    
        cdt = r(4);
        H_k = zeros(length(sats),4);
        h_k = zeros(length(sats),1);
        range_k = zeros(length(sats),1);
        for i = 1:length(sats)
            xs_vec = [xs(k,i);ys(k,i);zs(k,i)];        
            H_k(i,:) = [(ex_vec - xs_vec)'./norm(ex_vec-xs_vec) , 1 ];
            h_k(i,:) = norm(ex_vec - xs_vec) + cdt;
            range_k(i) = rs(k,i);
        end

        deltax = (H_k'*H_k)\H_k'*(range_k-h_k);
        r = r + deltax;
        err = norm(deltax);
        count = count + 1;
    end
    count
    output(k,1:3) = r(1:3)';
end

Error = output - receiver_position;
figure
scatter3(Error(:,1), Error(:,2), Error(:,3),'.')
xlabel('z');
ylabel('y');
zlabel('x');
axis equal
%plot(Error)
%ylim([-1,1])
