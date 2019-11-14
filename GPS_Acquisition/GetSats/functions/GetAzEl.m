function [Az,El]=GetAzEl(Xsv,Ysv,Zsv,Rx)
Az=zeros(size(Xsv));
El=Az;
[lat,lon]=ECEF_to_LLA(Rx);
for i=1:size(Xsv,1)
    for j=1:size(Xsv,2)
        r_Rx_Svij=[Xsv(i,j);Ysv(i,j);Zsv(i,j)]-Rx;
        [Elij,Azij]=LLA_to_ENU(lat,lon,r_Rx_Svij);
        Az(i,j)=Azij;
        El(i,j)=Elij;
    end
end
end
    
function [lat,lon]=ECEF_to_LLA(r)
% This function converts ECEF coordinates to WGS 84 Longitude Latitude
% Source for algorithm: B. Bowring (1976). "Transformation from Spatial to Geographical Coordinates". Surv. Rev. 23 (181): 323–327.
% Source for semi-major axis (a) and first eccentricity (e) parameters: http://www.arsitech.com/mapping/geodetic_datum/
a=6378137.0;
e=0.08181919092890624;
x=r(1);
y=r(2);
z=r(3);
p=sqrt(x^2+y^2);
% Initializes Newton-Raphson
k_new=1/(1-e^2);
err_threshold=0.0001;
err=1;
% Iterates Newton-Raphson
while err>err_threshold
    k_old=k_new;
    ci=(p^2+(1-e^2)*z^2*k_old^2)^(3/2)/(a*e^2);
    k_new=1+(p^2+(1-e^2)*z^2*k_old^3)/(ci-p^2);
    err=abs(k_new-k_old);
end
k=k_new;

lon=atan2(y,x); % Calculates longitude
lat=atan2(z*k,p);   % Claculates latitude
if lon>pi
    lon=lon-2*pi;   % Ensures longitude is [-pi,pi)
end
end

function [El,Az]=LLA_to_ENU(lat,lon,r_Rx_Sv)
% This function calculates the Elevation and Azimuth angles of the SV
% Source for algorithm: https://en.m.wikipedia.org/wiki/Molodensky_transformations#Molodensky_transformation
R=[-sin(lon),cos(lon),0;-sin(lat)*cos(lon),-sin(lat)*sin(lon),cos(lat);cos(lat)*cos(lon),cos(lat)*sin(lon),sin(lat)]; % Build the ECEF to ENU rotation matrix
r_ENU=R*r_Rx_Sv;    % Express the vector from the receiver to the SV in ENU coordinates    
El=asin(r_ENU(3)/norm(r_ENU));  % Calculates the elevation angle
Az=atan2(r_ENU(1),r_ENU(2));    % Calculates the azimuth angle
end