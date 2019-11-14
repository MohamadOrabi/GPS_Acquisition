clc;
ephemeris = read_rinex_nav('jplm1930.16n');

%time = 120786; % time vector for computing position
time = 120822;
PRN  = [27]; % PRN of satellite

% initialize constants and variables
svid = ephemeris(:,1);
m0   = ephemeris(:,2);
dn   = ephemeris(:,3);
e    = ephemeris(:,4);
a    = (ephemeris(:,5)).^2;
omg0 = ephemeris(:,6);
i0   = ephemeris(:,7);
w    = ephemeris(:,8);
odot = ephemeris(:,9);
idot = ephemeris(:,10);
cuc  = ephemeris(:,11);
cus  = ephemeris(:,12);
crc  = ephemeris(:,13);
crs  = ephemeris(:,14);
cic  = ephemeris(:,15);
cis  = ephemeris(:,16);
toe  = ephemeris(:,17);
iode = ephemeris(:,18);
GPS_week = ephemeris(:,19);
toc=ephemeris(:,20);
af0= ephemeris(:,21);
af1= ephemeris(:,22);
af2= ephemeris(:,23);
TGD=ephemeris(:,24);

meu = 3.986005e14;         % earth's universal gravitational [m^3/s^2]
odote = 7.2921151467e-5;   % earth's rotation rate (rad/sec)
lightspeed = 2.99792458e8; % speed of light (m/s)

F = -4.442807633e-10; % Constant, [sec/(meter)^(1/2)]

indx1 = find(svid == PRN(1));

% compute positions for single time
tsat = time(1);
for j = 1:4
    [blah indx2] = min(abs(tsat-toe(indx1)));
    if tsat-toe(indx1(indx2)) < 0 
        if indx2 == 1
        else
            indx2 = indx2 - 1;
        end
    end
    index = indx1(indx2);
    n0 = sqrt(meu/a(index)^3);
    t = tsat-toe(index);
    n = n0 + dn(index);
    m = m0(index) + n*t;
    
    m_dot=n;                                  % Calculate Velocity
    
    E = kepOrb2E(m,e(index));
    
    %Compute relativistic correction term
    dtr = F * e(index) * sqrt(a(index)) * sin(E);
    
    % Compute satellite clock correction    
    clkCorr= (af2(index) * (tsat-toc(index)) + af1(index)) * (tsat-toc(index)) + ...
        af0(index);
    
    t = t - clkCorr;
    
    E_dot=m_dot/(1-e(index)*cos(E));          % Calculate Velocity    
    
    v = atan2(sqrt(1-e(index)^2)*sin(E),cos(E)-e(index));
     
    v_dot=sin(E)*E_dot*...                    % Calculate Velocity
        (1+e(index)*cos(v))/(sin(v)*(1-e(index)*cos(E)));  
    
    phi = v + w(index);
    
    phi_dot=v_dot;                            % Calculate Velocity
    
    du = cus(index)*sin(2*phi) + cuc(index)*cos(2*phi);
    dr = crs(index)*sin(2*phi) + crc(index)*cos(2*phi);
    di = cis(index)*sin(2*phi) + cic(index)*cos(2*phi);
    
    du_dot=2*(cus(index)*cos(2*phi)-cuc(index)*sin(2*phi))*phi_dot; % Calculate Velocity
    dr_dot=2*(crs(index)*cos(2*phi)-crc(index)*sin(2*phi))*phi_dot; % Calculate Velocity
    di_dot=2*(cis(index)*cos(2*phi)-cic(index)*sin(2*phi))*phi_dot; % Calculate Velocity
        
    u = phi + du;
    r = a(index)*(1-e(index)*cos(E)) + dr;
    i = i0(index) + di + idot(index)*t;
    
    u_dot=phi_dot+du_dot;                         % Calculate Velocity
    r_dot=a(index)*e(index)*sin(E)*E_dot+dr_dot;  % Calculate Velocity
    i_dot=idot(index)+di_dot;                     % Calculate Velocity

    xp = r*cos(u);
    yp = r*sin(u);
    
    xp_dot=r_dot*cos(u)-r*sin(u)*u_dot;           % Calculate Velocity
    yp_dot=r_dot*sin(u)+r*cos(u)*u_dot;           % Calculate Velocity

    omg = omg0(index) + (odot(index) - odote)*t - odote*toe(index);
    
    omg_dot=odot(index) - odote;                  % Calculate Velocity
    
    XS = xp*cos(omg) - yp*cos(i)*sin(omg);
    YS = xp*sin(omg) + yp*cos(i)*cos(omg);
    ZS = yp*sin(i);
    
    spos_lla = ecef2lla([XS YS ZS]')
end

% ------------------------------------------------
function E = kepOrb2E(M,e)
% Inputs:  - mean anomaly in radians
%          - eccentricity
% Output: Eccentric anomaly

if (-pi < M < 0) | (M > pi)
    E = M - e;
else
    E = M + e;
end

check = 1;

while check > 10e-10
    E_new = (E + (M - E + e * sin(E))/(1 - e * cos(E)));
    check = abs(E_new - E);
    E = E_new;
end
end

