function [ G, H, E, sigma, B, eta ] = MaxEntPro( qs, Ts, Rn, theta, Thmax, Thmin, Ids, rho_s, par, OPT1 )

% disp('Function in m')
% ----------------------------------------------------------------------- %

% ##INPUT
%qs     --> Specific humidity   [kg/kg]
%Ts     --> Surface temperature [K]
%Rn     --> Net radiation       [W/m^2]
%OPT1   --> OPT1 = 1: assume canopy
%           OPT1 = 2: assume bare soil
%           OPT1 = 3: assume water surface
%           OPT1 = 4: assume snow surface

% ##OUTPUT
%G      --> Ground Heat Flux    [W/m^2]
%H      --> Sensible Heat Flux  [W/m^2]
%E      --> Latent Heat Flux    [W/m^2]
%sigma  --> parameter; refer to [Wang et al 2011]
%B      ->> Bowen ratio: [H/E]

% ##CONSTANTS
Rv=461.5;              %Gas constant for water vapor[J/kg/K]
Lv=2.5E+06;            %latent heat of vaporization [J/kg]
rho_w = 1000;          %density of water            [kg/m^3]
rho_a = 1.225;         %density of air              [kg/m^3]
Cp=1005;               %Heat capacity               [J/kg/K]
g = 9.81;              %gravity accelaration
Av = Lv*Lv/Rv/Cp;      %

% ----------------------------------------------------------------------- %

% ##
idx = 1/6;             %index
T0 = 300;              %Reference temperature       [K]
rhoCp = Cp*rho_a;      % 
% albedo = 0.05;       %albedo of liquid water
% alfa = 0.75;         %1 or 0.75
% beta = 4.7;
% gam1 = 15; 
% gam2 = 9;
C01 = 1.7;             %unstable C0
C02 = 1.2;             %stable C0
ka = 0.4;              %kappa
if OPT1 == 1
    z = 10;                %at                          [m]
elseif OPT1 == 2 
    z = 0.2;
elseif OPT1 == 3 || OPT1 == 4
    z = 2;
end

if isempty(par)
    alpha = 1; eta0  = 10/3; beta  = 2;
else
    alpha = par.alpha; eta0  = par.eta0; beta  = par.beta;
end

[m,n] = size(Rn);
eta = NaN(m,n);
if isempty(theta)
    eta = ones(m,n);
else
    if isempty(Thmin) 
        Thmin = prctile(theta,1,2); %min(theta,[],2);
    end
    if isempty(Thmax)
        Thmax = prctile(theta,99,2);%max(theta,[],2);
    end
    %Ths = max(theta,[],2);
    Ths = Thmax;
    Ths = repmat(Ths,1,n);
    Thmax = repmat(Thmax,1,n);
    Thmin = repmat(Thmin,1,n);
    if OPT1 ~= 3 
        if isempty(theta)
            eta = ones(m,n);
        else
            if isempty(Thmin) 
                Thmin = prctile(theta,1,2);
                Thmin = repmat(Thmin,1,n);
            end
            if isempty(Thmax)
                Thmax = prctile(theta,99,2);
                Thmax = repmat(Thmax,1,n);
            end
    
            for i = 1 : m
                for j = 1 : n
                    dumeta = eta0*(theta(i,j) - Thmin(i,j))/(Thmax(i,j) - Thmin(i,j)).*alpha;
                    if dumeta > 1
                        eta(i,j) = 1;
                    elseif dumeta <= 0.25
                        eta(i,j) = 0.25;
                    else
                        eta(i,j) = dumeta;
                    end
                end
            end
        end
    end
end
% ## Thermal inertial
rho_s(rho_s > rho_w ) = rho_w ;
Iw    = sqrt(rho_w*4.18*10^3*0.58); %thermal inertia of water
Cp_i  = 2000; % heat capacity of ice
ks_s  = 0.3;  % thermal conductivity of snow
ks_i  = 2.0;  % thermal conductivity of ice
if isempty(rho_s)
    rho_s = 500; % Default snow density
end
Isnow = sqrt(rho_s.*Cp_i.*ks_s);     % thermal inertia of snow, rho_s: snow density

Ih1 = rhoCp*C01*sqrt(ka*z)*(ka*z*g/rhoCp/T0)^idx;    % unstable I0
Ih2 = rhoCp*C02*sqrt(ka*z)*(ka*z*g/rhoCp/T0)^idx;    % stable I0
if OPT1 == 2
    if isempty(Ids)
        Ids = 800; %Default thermal inertia of dry soil [600-1000] (Farouki 1982; Wang et al. 2010) Huang et al 2016 use 800
    end
    Is  = sqrt(Ids^2 + theta.*Iw^2);
    %Is  = sqrt(Ids.^2 + theta.*Iw.^2);
    %Is  = sqrt(Ids.^2 + Thmax.*Iw.^2);
    Ih1 = rhoCp*C01*sqrt(ka*z)*(ka*z*g/rhoCp/T0)^idx;    % unstable I0
    Ih2 = rhoCp*C02*sqrt(ka*z)*(ka*z*g/rhoCp/T0)^idx;    % stable I0
    a1  = Is./Ih1;            %unstable soil
    a2  = Is./Ih2;            %stable soil
end

if OPT1 == 3
    Iwsi = Iw;
elseif OPT1 == 4
    Iwsi = Isnow;
end

Rn = double(Rn);
sigma = NaN(m,n); 
B = NaN(m,n); 
H = NaN(m,n); 
E = NaN(m,n);
G = NaN(m,n);
if OPT1 == 1
    qsat = Calc_qsat( Ts, -50 );
    %qs   = qsat;
    %qsat = qs;
    for i = 1 : m
        for j = 1 : n
            sigma(i,j) = eta(i,j)*Lv^2.*qsat(i,j)./(Cp.*Rv.*Ts(i,j).^2);
            B(i,j)     = 6.*((1+11./36.*sigma(i,j)).^0.5-1);
            E(i,j)     = Rn(i,j)./(1+1./B(i,j));
            H(i,j)     = Rn(i,j)./(1+B(i,j));
            G(i,j)     = 0;
        end
    end
elseif OPT1 == 2
    Ts
     esat = Calc_esat( Ts );
     
     qsat = e2q(esat,1e5);
     qss  = (theta./Ths).^beta.*qsat;
    for i = 1 : m
        for j = 1 : n
            if isnan(Ts(i,j)) || isnan(Rn(i,j)) || isnan(qss(i,j)) || isnan(Is(i,j)) %|| isnan(eta(i,j))

            else
                % Av = Lv*Lv/Rv/Cp; 
                x1         = Av/(Ts(i,j))^2;
                sigma(i,j) = x1*qss(i,j);
                B(i,j)     = 6*((1+11/36*sigma(i,j))^0.5-1);
                bdsg       = B(i,j)/sigma(i,j);
                x0         = Rn(i,j)/2;

                if Rn(i,j) >= 0 % unstable condition
                   H(i,j) = fzero(@(y) ((1+B(i,j))*y-Rn(i,j))*abs(y)^idx+bdsg*a1(i,j)*y, x0);
                   E(i,j) = B(i,j)*H(i,j); 
                   G(i,j) = Rn(i,j) - H(i,j) - E(i,j);
                else            % stable condition
                   H(i,j) = fzero(@(y) ((1+B(i,j))*y-Rn(i,j))*abs(y)^idx+bdsg*a2(i,j)*y, x0);
                   E(i,j) = B(i,j)*H(i,j); 
                   G(i,j) = Rn(i,j) - H(i,j) - E(i,j);
                end
            end
        end
    end 

elseif OPT1 == 3 || OPT1 == 4
    a1  = Iwsi./Ih1;         %unstable soil
    a2  = Iwsi./Ih2;         %stable soil

    esat = Calc_esat( Ts );
    qss  = e2q(esat,1e5);

    for i = 1 : m
        for j = 1 : n
            if isnan(Ts(i,j)) || isnan(Rn(i,j)) %|| isnan(qs(i,j))

            else
                x1         = Av/(Ts(i,j))^2;
                sigma(i,j) = x1*qss(i,j);
%                 sigma(i,j) = Lv * qss(i,j)/Rv/(Ts(i,j)^2);
                B(i,j)     = 6*((1+11/36*sigma(i,j))^0.5-1);
                bdsg       = B(i,j)/sigma(i,j);
                x0         = Rn(i,j)/2;
                if Rn(i,j) >= 0 % unstable condition
                   H(i,j) = fzero(@(y) ((1+B(i,j))*y-Rn(i,j))*abs(y)^idx+bdsg*a1*y, x0);
                   E(i,j) = B(i,j)*H(i,j); 
                   G(i,j) = Rn(i,j) - H(i,j) - E(i,j);
                else            % stable condition
                   H(i,j) = fzero(@(y) ((1+B(i,j))*y-Rn(i,j))*abs(y)^idx+bdsg*a2*y, x0);
                   E(i,j) = B(i,j)*H(i,j); 
                   G(i,j) = Rn(i,j) - H(i,j) - E(i,j);
                end
            end
        end
    end 
end

eta0
