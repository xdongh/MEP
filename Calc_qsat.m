function qs = Calc_qsat(tair,phi_s)
    Rv=461.5;              %Gas constant for water vapor[J/kg/K]
    Lv=2.5E+06;            %latent heat of vaporization [J/kg]
    T0=273.15;
    g  = 9.81;

    esat = Calc_esat( T0 );
    
    qs = 0.622.*esat./1e5.*exp(Lv/Rv.*(1/T0 - 1./tair)).*exp(-g*phi_s./Rv./tair);
end