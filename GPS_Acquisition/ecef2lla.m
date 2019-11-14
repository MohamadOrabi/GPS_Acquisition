function [R_lla] = ecef2lla(R)
% This function converts ECEF coordinates to WGS 84 Longitude Latitude
% Source for algorithm: B. Bowring (1976). "Transformation from Spatial to Geographical Coordinates". Surv. Rev. 23 (181): 323–327.
% Source for semi-major axis (a) and first eccentricity (e) parameters: http://www.arsitech.com/mapping/geodetic_datum/
a = 6378137.0;
e = 0.08181919092890624;
N = size(R,2);
for ii=1:N
    r = R(:,ii);
    x = r(1);
    y = r(2);
    z = r(3);
    p = sqrt(x^2+y^2);
    k_new = 1/(1-e^2);
    err_threshold = 0.0001;
    err = 1;

    while err>err_threshold
        k_old = k_new;
        ci = (p^2+(1-e^2)*z^2*k_old^2)^(3/2)/(a*e^2);
        k_new = 1+(p^2+(1-e^2)*z^2*k_old^3)/(ci-p^2);
        err = abs(k_new-k_old);
    end
    k = k_new;

    lon(ii) = atan2(y,x);
    lat(ii) = atan2(z*k,p);
    R_N = a/(sqrt(1 - e^2*sin(lat(ii))^2));
    if lon(ii) > pi
        lon(ii) = lon(ii) - 2*pi;
    end
    alt(ii) = p/cos(lat(ii)) - R_N;
end
lat = lat/pi*180;
lon = lon/pi*180;
R_lla = [lat; lon; alt];
end