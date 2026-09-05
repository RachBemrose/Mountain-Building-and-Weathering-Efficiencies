function dy = compute(t,y)
load setuprun.mat
dy = zeros(7,1); 

%%% ----- define globals 

global F_sub_0
global F_sio3w_0
global F_co3w_0
global F_pel_0
global F_dep_0
global F_ridge_0
global F_meta_0
global F_hydro_0
global F_wedge_0
global Cdeep_old
global OA_0
global CS_0
global OC_0
global MAN_0
global temp_0
global phi_0
global trackermatrix
global trimmed_tracker
global ticker

%% Stellar Evolution Model

if strcmp(solar_evolution,'on')
    if strcmp(temp_function, 'jacob_clim')
        S = (LumiFit_0(t,star_size)) .* (1/(semi_major_axis^2) .* ((1 - e^2)^(1/2))); % Corrected solar luminosity fits from Rushby et al. (2013) - between 0.45 M and 1 M and distance scaling for jacob's interpolator
        solar_L = 3.846e26;
        solar_flux = S .* ((solar_L) / ((4*pi) * (semi_major_axis* 1.496e11)^2) * ((1 - e^2)^(1/2)));
    else 
    
        S = LumiFit_0(t,star_size); % Corrected solar luminosity fits from Rushby et al. (2013) - between 0.45 M and 1 M.
        solar_L = 3.846e26;
        solar_flux = S .* ((solar_L) / ((4*pi) * (semi_major_axis* 1.496e11)^2) * ((1 - e^2)^(1/2)));
    end
else
    S = linspace(0.5,0.75,1e4);
%     S = 1;
    solar_L = 3.846e26;
    solar_flux = S .* ((solar_L) / ((4*pi) * (semi_major_axis* 1.496e11)^2) * ((1 - e^2)^(1/2)));
end


%% Planetary Albedo

if strcmp(variable_albedo,'on')
    albedo = ((1.4891 - (0.0065979 .* (y(5)))) + ((8.567e-6) .* ((y(5)) .^2))); % Caldiera and Kasting 1992
elseif strcmp(variable_albedo,'icetransition')
    albedo_0 = 0.3;
    albedo_i = 0.6;
    T_0 = 278;
    T_I = 273;
        if y(5) <= T_I
            albedo = albedo_i;
        elseif y(5) > T_I && y(5) < T_0 
            albedo = albedo_0 + ( albedo_i - albedo_0) .* ( (y(5) - T_0).^2 / (T_I - T_0).^2 );
        elseif y(5) >= T_0
            albedo = albedo_0; %above 278, use albedo_0
        end
elseif strcmp(variable_albedo,'icetransition_variable')
    albedo_0 = 0.3;
    albedo_i = 0.6;
    T_0 = 278;
    T_I = 273;
        if y(5) <= T_I
            albedo = albedo_i;
        elseif y(5) > T_I && y(5) < T_0 
            albedo = albedo_0 + ( albedo_i - albedo_0) .* ( (y(5) - T_0).^2 / (T_I - T_0).^2 );
        elseif y(5) >= T_0
            albedo = ((1.4891 - (0.0065979 .* (y(5)))) + ((8.567e-6) .* ((y(5)) .^2))); % Caldiera and Kasting 1992; above 278, use variable albedo, definitely works more effectively.
        end
elseif strcmp(variable_albedo,'jacob_alb')
     albedo = jacob_alb(y(7),S); % 
elseif strcmp(variable_albedo,'jacob_alb_low')
    albedo = jacob_alb_lowpressure(y(7),S); %
elseif strcmp(variable_albedo,'ocean_continent_split')
    albedo = ocean_frac_energy(ocean_fraction, albedo_0);
else
    albedo = albedo_0;
end
%% Effective Temperature
% 
 eff_temp = (  ((S .* solar_flux) .* (1 - albedo)) / (4.* 5.67e-8)    ).^(0.25);

%% H2O Greenhouse

if strcmp(water_vapour_greenhouse, 'on')
    H2O_greenhouse = 29.7 + 0.5 .* (y(5) - 285); %WHAK 1982 T0 = 285
else
    H2O_greenhouse = 0;
end

%% Surface Pressure

surface_pressure = 0.937.*(Planet_Radius)^2.4; %(for less than 5 Earth masses)

%% CH4 Greenhouse

if strcmp(methane_greenhouse, 'on')   
    % - Temperature dependent methane emissions from the biosphere - %
        if strcmp(biosphere,'on')
            u_max = 0.75; %0.75 (doublings per day) or 0.38
            T_opt = 28; % 28 or 9
            delta_T = 18; %28 or 18
            T_ref = 1400;
            cyano = u_max - u_max .* ( ( (T_opt - (y(5)-273)) / delta_T ).^2 ) .* exp( - ( ((y(5)-273) - T_opt) / T_ref ).^3 ); %temperature dependent growth curve
            psi = 611; % mol yr-1. 1273 for 'Lovelockian' or 611 to ensure 400e12 mol yr-1 at 18 C.
            flux = (psi * cyano) * 1e12; % mols of CH4 emitted per year by the biosphere. 
        else
            flux = 400e12;    
        end

    escape = 3.7e-5 .* flux;
    fch4 = flux - escape;
    M_0 = 700; % (http://cdiac.ornl.gov/pns/current_ghg.html)
    M = fch4/1.76e10;
    N_0 = 1;
    N = 1;
    f_0 = 0.47 * log(1+2.01e-5 .* (M_0.*N_0).^0.75 + (5.31e-15).*M_0.*(M_0.*N_0).^1.52);
    f = 0.47 * log(1+2.01e-5 .* (M.*N).^0.75 + (5.31e-15.*M) .*(M.*N).^1.52);
    CH4_forcing = 0.036*(sqrt(M) - sqrt(M_0)) - ((f - f_0)); %IPCC (http://www.grida.no/publications/other/ipcc_tar/?src=/climate/ipcc_tar/wg1/222.htm) http://www.climatechange2013.org/images/report/WG1AR5_Ch08SM_FINAL.pdf
    CH4_greenhouse = CH4_forcing * Climate_sensitivity;
else
    CH4_greenhouse = 0;
end
    

%% CO2 Greenhouse Parameterisations
%The IPCC formulations are actually no longer used, but were incorporated
%in order to test the radiative-convective climate model. 

if strcmp(Greenhouse_Parameterisation,'CK')
    CO2_greenhouse = 815.17  + (4.895*10^7)*(y(5)^-2) -  (3.9787*10^5)*(y(5)^-1)  -6.7084*((log10(  y(7)  ))^-2) + 73.221*((log10( y(7)  ))^-1) -30882*(y(5)^-1)*((log10(  y(7) ))^-1);  
elseif strcmp(Greenhouse_Parameterisation,'IPCC1')
    C = y(7) * 1e6; %into ppm
    CO2_IPCC = 5.35 .* log(C/280); % (http://www.grida.no/publications/other/ipcc_tar/?src=/climate/ipcc_tar/wg1/222.htm)
    CO2_greenhouse = Climate_sensitivity .* CO2_IPCC; %http://www.grida.no/climate/ipcc_tar/wg1/pdf/tar-06.pdf (6 of 68)
elseif strcmp(Greenhouse_Parameterisation,'IPCC2')
    C = y(7) * 1e6; %into ppm
    CO2_IPCC = 4.841 .* log(C/280) + (0.0906 .* (sqrt(C) - sqrt(280)));
    CO2_greenhouse = Climate_sensitivity .* CO2_IPCC; %http://www.grida.no/climate/ipcc_tar/wg1/pdf/tar-06.pdf (6 of 68)
elseif strcmp(Greenhouse_Parameterisation,'IPCC3')
    C = y(7) * 1e6; %into ppm
    g_0 = log (1+1.2 .* 280 + 0.005 .* 280 .^(2) + 1.4e-6 .* 280 .^ (3));
    g= log (1+1.2 .* C + 0.005 .* C .^(2) + 1.4e-6 .* C .^(3));
    CO2_IPCC = 3.35 .* (g - g_0); 
    CO2_greenhouse = Climate_sensitivity .* CO2_IPCC; %http://www.grida.no/climate/ipcc_tar/wg1/pdf/tar-06.pdf (6 of 68)
elseif strcmp(Greenhouse_Parameterisation,'WHAK')
    CO2_greenhouse = 4.32 .* log(y(7)/280e-6); % a doubling of CO2 would result in a ~3K increase in temperature, as per the IPCC Fifth Assessment Report 
end

CO2_greenhouse(CO2_greenhouse<0) = 0; %prevents the CO2 greenhouse parameterisation from going negative.



%% Global Average Temperature

if strcmp(temp_function,'opacity_based')
    temp =  eff_temp .* (1 + 0.75 .* ( 0.0297 .* exp ( (0.0111 .* y(5)) + (0.2190 .* log(y(7))) ) ) ); % Ben and Colin's opacity based parameterisation
elseif strcmp(temp_function,'jacob_clim')
%%%%% altered here BM
%     temp = (merged_clim_alb(y(7),S)) ;
    %%%%% altered here BM
    temp = (jacob_clim_alb_fixed(y(7),S)) ;
%     temp = (jacob_clim_alb(y(7),S)) ;
    column_integrated_H2O = jacob_H2O(y(7),S);
    OLR = jacob_olr(y(7),S);
elseif strcmp(temp_function,'jacob_clim_low')
     temp = (jacob_clim_lowpressure(y(7),S)) + 273; 
     column_integrated_H2O = jacob_H2O_lowpressure(y(7),S);
     OLR = jacob_olr_lowpressure(y(7),S);
else
    temp = (  eff_temp + H2O_greenhouse + CO2_greenhouse + CH4_greenhouse) ; 
    column_integrated_H2O = jacob_H2O_lowpressure(y(7),S);
     OLR = jacob_olr_lowpressure(y(7),S);
end



%% pCO2 %%

if strcmp (dynamic_phi, 'on')
    temp_present = 288;
    kc2 = -0.0448;
   %phi = phi_0 * exp( kc2 *( temp_present - temp ) ) * (  y(1) * ((1-y(6))^2) / ( ((1-phi_0)^2)*OA_0 ) .* (y(1)/OA_0) );
    phi = phi_0 * exp( kc2 *( temp_present - temp ) ) ;
    long_term_stable_CO2 = (280e-6); %1PAL
    CO2atm = ( (y(1) * phi) / (phi_0 * OA_0) ) * (long_term_stable_CO2);   %(pCO2) From Bergman et al. 2004 (Noam's thesis, equation 5.14 pg.138)
    CO2relative = CO2atm / (280e-6);
    phi = phi_0;
elseif strcmp (dynamic_phi, 'Kump')   
    phi = phi_0;
    long_term_stable_CO2 = (280e-6); %1PAL
    CO2_carbon_dependence = ((y(1)/OA_0).^2 .* long_term_stable_CO2); % Kump & Arthur (1999) carbon partitioning reservoir size dependence 
        TB = 27; %Omta et al. (2011) carbon partitioning temperature dependence
        delta_T = y(5) - 261; %this 288 is wrong! 261 K gives 1 PAL at 288 K. 288 here gives 0 * PAL at 288, thereby underestimating pCO2 significantly.
    CO2_temp_dependence = (((CO2_carbon_dependence/1e6) .* delta_T) / TB);
    CO2atm =  (CO2_carbon_dependence + CO2_temp_dependence);% / Planet_Radius.^2;
    CO2relative = CO2atm / (280e-6);  
 else
   phi = phi_0;
end


           
%% C deep
% if t > 3e9
%   Cdeep = ( -0.4*t / (3e9) ) + 0.75 ; % Proportion of carbon that reaches the deep mantle. Sleep and Zahnle 2001
% else
%   Cdeep = Cdeep_old;
% end

Cdeep = 0.75;

% if t > 3e9
%    DS = ((3.35e-10) .* (t)) + 1; % Source depth. Sleep and Zahnle 2001
% else
%    DS = 2;
% end

DS = 2;

DS_relative = DS/2;

%% heat flow

if strcmp(heat_flow_evolution, 'on')
    Q = geothermal_heat_flux(t,heat_flow,Planet_Mass); % final function varies radiogenic heating with planet mass (Seagar et al, pg 385). seafloor creation rates  go as the heat flow squared
else
    Q = inline('1 + 0*x*y*z');
end

%% Biotic Weathering Enhancement

CO2_C3_plant_limit = 0.000150;
CO2_C4_plant_limit = 0.000010;

bio_enhancement = K;

% if strcmp(time_dependent_biotic_enhancement, 'on')
%     if le(t, plant_colonisation_time)
%        bio_enhancement = K;
%     elseif ge(t, plant_colonisation_time) && ge(y(7), CO2_C3_plant_limit)
%             u_max = 0.75; %0.75 (doublings per day) or 0.38
%             T_opt = 273+28; % 28 or 9
%             del_T = 273+28; % 28 or 18
%             T_ref = 1400;
%             C3 = 389.3; %Still et al. (2003)
%             C4 = 18.6;
%             C_ratio = C4/C3;
%       growth = u_max - u_max .* ( ( (T_opt - (y(5)-273)) / del_T ).^2 ) .* exp( - ( ((y(5)-273) - T_opt) / T_ref ).^3 ); %temperature dependent growth curve
%      if ge(y(7), CO2_C3_plant_limit)
%         bio_enhancement = (K + ((1 - ocean_fraction) * growth));
%      else
%         bio_enhancement = (K + (1 - ocean_fraction * growth)* (C_ratio));
%      end
%     else
%      bio_enhancement = K;
%     end
% else
% bio_enhancement = K;   
% end 

% temp_C = 0:1:50;
B_t_e_23 = 0.8352;
B_t_e_25 = 0.3442;
temp_dependent_productivity_23 = B_t_e_23^-1 .* (1 - (((y(5)-273) - 25)/25).^2);
temp_dependent_productivity_25 = B_t_e_25^-1 .* (((1 - (y(5)-273)/50).^(1/2)) - ((1 - (y(5)-273)/50).^2));


% CO2 dependence
pCO2_e = 280;
alpha_min = 10;
alpha_1_2_24 = 181; %ppm
alpha_1_2_25 = 124;
CO2_dependent_productivity_24 = ((y(7) - alpha_min) .* (pCO2_e - alpha_min + alpha_1_2_24)) / ((y(7) - alpha_min + alpha_1_2_24) .* (pCO2_e - alpha_min));
CO2_dependent_productivity_25 = ((y(7) - alpha_min) .* (pCO2_e - alpha_min + alpha_1_2_25)) / ((y(7) - alpha_min + alpha_1_2_25) .* (pCO2_e - alpha_min));

% NPP (eq. 22)

epsilon = 0.25;
B = epsilon + (1-epsilon) .* temp_dependent_productivity_23 .* CO2_dependent_productivity_24;
B_25 = epsilon + (1-epsilon) .* temp_dependent_productivity_25 .* CO2_dependent_productivity_25; %probably use this one going forward, see Honing page 7 for discussion
% figure;plot(temp_C,B)
% hold on
% plot(temp_C,B_25)


%% weathering

if strcmp(silicate_weathering,'SZ')
    fco2 = inline ( ' v .* (w ^ x) .* exp((y - 288) / z )'); % Equation 8 in Zahnle and Sleep 2002
    
    C = 0.176;
    evap = (Planet_Radius.^2) .* ((1 - C) .* ocean_fraction); %write comments! WHat is this?
    g = (6.67e-11 .* (Planet_Mass * 5.98e24)) ./ (Planet_Radius * 6371000).^2;
    %rachel_weathering_function = 32.82.*Planet_Radius.^3 - 88.38.*Planet_Radius.^2 + 76.35.*Planet_Radius - 20.25;
    WR = Weathering_Efficiency(Planet_Radius);

%     F_sio3w = (fco2(F_sio3w_0,CO2relative,B,temp,B_kinetic)) * bio_enhancement;
 %   F_sio3w = (fco2(F_sio3w_0,CO2relative,B,temp,B_kinetic)) .* bio_enhancement .* erosion_rate * evap;
  %  F_sio3w = (fco2(F_sio3w_0,CO2relative,B,temp,B_kinetic)) .* rachel_weathering_function;
%     F_sio3w = (fco2(F_sio3w_0,CO2relative,B,temp,B_kinetic)) .* bio_enhancement;
% updated F_sio3w = (fco2(F_sio3w_0,CO2relative,B,temp,B_kinetic)) .* Weathering_Efficiency(Planet_Radius)
WR = Weathering_Efficiency(Planet_Radius) ./ Weathering_Efficiency(1);
F_sio3w = 0.5.*(fco2(F_sio3w_0,CO2relative,B,temp,B_kinetic)) .* WR;

elseif strcmp(silicate_weathering,'COPSE')
    ks1 = 7537.69;
    ks2 = 0.03;
    fco2 = inline('u .* (exp(v .* ((w - 288)/(w .* 288))) .* (exp(x .* (w - 288))).^0.65 .* sqrt(y) )'); %relative concentration i.e. PAL
    F_sio3w = (fco2(F_sio3w_0,ks1,temp,ks2,CO2relative)) * bio_enhancement;
    
end

gco2=inline( ' exp(0.05*(x-288)).*((y).^(z))'); %Dependence of carbonate weathering on temperature



%% Flux calculations

F_sub = F_sub_0 .* (Q.^2) .* (y(3)/OC_0); %
F_sub_arc = (1-Cdeep) * F_sub; % Proportion of subducted C is expelled through arc volcanoes
F_sub_man = Cdeep * F_sub; %Proportion of subducted C goes to Mantle
% 
% F_co3w = F_co3w_0 .* (gco2(temp,CO2relative,G));
F_co3w = F_co3w_0 .* (gco2(temp,CO2relative,G)) .* rachel_weathering_function .* evap .* (y(2)/CS_0);
% F_co3w = F_co3w_0 .* (gco2(temp,CO2relative,G)) .* (y(2)/CS_0) ;

fraction_pel = 0.542; %dep and pel represent 24e12 present day, dep = 11e12 and pel = 13e12.
F_pel = ( (F_co3w + F_sio3w) * fraction_pel); %Pelagic carbon (from OA to OC reservoir)
F_dep = ( (F_co3w + F_sio3w) * (1 - fraction_pel) ); %Carbon deposition (from OA reservoir to CON reservoir)


F_ridge = F_ridge_0 .* (y(4)/MAN_0) * (Q.^2) .* DS_relative;
F_meta = (F_meta_0 .* (Q.^2) .* (y(2)/CS_0)); %"The metamorphic flux Fmeta is assumed proportional to the product of the heat flow and the continental carbonate reservoir, QRcon.

F_hydro = F_hydro_0 .* CO2relative^A .* (Q.^2);
F_wedge = F_wedge_0 .* (y(3)/OC_0); % .* (Q(t,heat_flow,Planet_Mass).^2);

%% ----- Reservoir Differentials----- %%

    dy(1) =  F_ridge + F_co3w + F_meta + F_sub_arc - F_pel - F_dep - F_hydro; % ocean/atmosphere
% dy(1) = 0;

    dy(2) = F_dep + F_wedge - F_meta - F_co3w; % R_con continental sediments
% dy(2) = 0;

    dy(3) =  F_hydro + F_pel - F_sub - F_wedge; % oceanic crust
% dy(3) = 0;

    dy(4) = F_sub_man - F_ridge; % mantle 
% dy(4) = 0;
    

%% ----- Iterative Parameterisations ----- %%   

     dy(5) = temp - y(5) ; %%% global avg. surf. temp.
%     
     dy(6) = phi - y(6);
    
     dy(7) = CO2atm - y(7);
 
    %  %% - trackers
    % ticker = ticker + 1;
    % trackermatrix(ticker,1) = t;
    % trackermatrix(ticker,2) = temp;
    % trackermatrix(ticker,3) = phi;
    % trackermatrix(ticker,4) = albedo;
    % trackermatrix(ticker,5) = Q;
    % trackermatrix(ticker,6) = S;
    % trackermatrix(ticker,7) = CO2atm;
    % trackermatrix(ticker,8) = F_ridge;
    % trackermatrix(ticker,9) = F_sub;
    % trackermatrix(ticker,10) = F_co3w;
    % trackermatrix(ticker,11) =  F_pel;
    % trackermatrix(ticker,12) =  F_dep;
    % trackermatrix(ticker,13) = F_meta;
    % trackermatrix(ticker,14) = F_hydro;
    % trackermatrix(ticker,15) = F_wedge;
    % trackermatrix(ticker,16) = F_sio3w;
    % trackermatrix(ticker,17) = CO2_greenhouse;
    % trackermatrix(ticker,18) = CH4_greenhouse;
    % trackermatrix(ticker,19) = eff_temp;
    % trackermatrix(ticker,20) = H2O_greenhouse;
    % trackermatrix(ticker,21) = column_integrated_H2O;
    % trackermatrix(ticker,22) = OLR;
    % trackermatrix(ticker,23) = bio_enhancement;
    % trackermatrix(ticker,24) = B_25;
    % 
    % 
    % 
   


    format short G
    disp('Time                  Temp(K)     CO2(bars)')
    disp([t,temp,CO2atm,Planet_Radius])
end
    