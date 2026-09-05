       %%% --- Optimised S and Z Model --- %%%
%%% --- Andrew Rushby, Ben Mills, Martin Johnson and Mark Claire (2013)
%%% ---%%%
global trimmed_tracker

% model_version = runscript_version;

counter = 1;
for x = 0.5 % <--- this will be the variable that changes for each run; can be single value or array
% for x = 0.5 
%% PLANET & STAR OPTIONS     

%-planet-%
Planet_Property = 'Radius'; %choose Mass or Radius, which will scale the other. See code for relationships
  Planet_Mass = [3.93] ; %Earth Masses
    %or%
  Planet_Radius = x;
ocean_fraction = 0.7; %<--- scales the OC and CS reservoirs, not very robust (yet))
variable_albedo = 'jacob_alb'; % jacob_alb, or 'ocean_continent_split', 'icetransition_variable'
    albedo_0 = 0.31; %initial albedo
heat_flow_evolution = 'on'; %age-dependant geothermal heat flux
dynamic_phi = 'Kump'; %Ocean/atmosphere CO2 partitioning

%% PLANET SCALING %%

if strcmp(Planet_Property, 'Radius')
    if le(Planet_Radius, 1)
    Planet_Mass = Planet_Radius .^ 3.268; %These are from Barnes et al. 2015 - HITE
    elseif gt(Planet_Radius, 1)
    Planet_Mass = Planet_Radius .^ 3.65;
    end 
elseif strcmp(Planet_Property, 'Mass')
    if Planet_Mass > 0.01 < 1
    a = 1;
    b_ = 0.306;

    elseif Planet_Mass > 1 < 10
    a = 1;
    b_ = 0.274;                 % Seager et al. 2007, originally Sotin et al. 2007.Valencia et al. 2006 has b between 0.267 and 0.272. Grasset has b between 0.252 and 0.285

    elseif Planet_Mass > 10 < 100
    a = 1.01;
    b_ = 0.2685;                %Grasset in Seager 2007. Mean value.
    
    Planet_Density = (1*(Planet_Mass.^(1-3*(b_)))); 
    Planet_Gravity = (1*(Planet_Mass.^(1-2*(b_)))); 
    Planet_Mantle = (1*(Planet_Mass.^(0.3094))); 
    end
end


%-stellar parameters-%

semi_major_axis = 0.552; %in AU
e = 0.18; %Eccentricity. Earth = 0.0167 
solar_evolution = 'on';
star_size = 0.78; %in solar masses
star_age = Agefit(star_size)*1e9; %main sequence lifetime is calculated from fits to stellar evolution model.

%% PLANETARY GREENHOUSE

Climate_sensitivity = 0.8; % min = 0.41 and max = 1.2
water_vapour_greenhouse = 'off'; % <--- needs updating.
Greenhouse_Parameterisation = 'IPCC2'; %IPCC2 probably best of the IPCC bunch.
methane_greenhouse = 'on'; %<--- constant 3K wamring. Will be updated later.
temp_function = 'jacob_clim'; %default(parameterised) or jacob_clim (for radiative-convective fits).
biosphere = 'off'; % <--- Only considered when methane_greenhouse is on. Not quite working yet, leave off.

%% WEATHERING OPTIONS

heat_flow = 0.2; % 0.7 for high geothermal heat flow, 0.2 for low heat flow.
A = 0.23; %(0.23, 0.4, 0.7, 1) dependence of carbonatization on amount of dissolved carbonate in sea water. Quantitive depletion corresponds to 1.
silicate_weathering = 'SZ'; %SZ or COPSE
B = 0.3; %silicate weathering sensitivity S&V: between 0.3 and 0.4. Default = 0.4. Berner's version = 0.5
B_kinetic = 13.7; %silicate weathering calibration constant for runoff kinetics. Default = 13.7. Berner's version = 11.1
G = 0.3; % carbonate weathering sensitivity
time_dependent_biotic_enhancement = 'off';
    plant_colonisation_time = 4.2e9;
    K = 0.5; %pre-plant dampening factor for terrestrial weathering. Default = 0.5;
    


%% RUN OPTIONS    

whenstart = -5e9; %start run
whenend = 23e9; %end run
%disp(star_age)

%% FIGURE OPTIONS
   working_variable_{counter} = x;
   cyclecolour(1) = 'r';
   cyclecolour(2) = 'g';
   cyclecolour(3) = 'b';
   cyclecolour(4) = 'c';
   cyclecolour(5) = 'm';
   cyclecolour(6) = 'y';
   cyclecolour(7) = 'r';
   cyclecolour(8) = 'r';
   cyclecolour(9) = 'g';
   cyclecolour(10) = 'b';
   cyclecolour(11) = 'k';
   cyclecolour(12) = 'c';
   cyclecolour(13) = 'm';
   cyclecolour(14) = 'y';
   cyclecolour(15) = 'g';
   cyclecolour(16) = 'g';
   cyclecolour(17) = 'b';
   cyclecolour(18) = 'c';
   cyclecolour(19) = 'm';
   cyclecolour(20) = 'y';
   figurecolour = cyclecolour(counter); 
 
   linechoice = '-';
   namechoice = num2str(working_variable_{counter});
   
 

    
    save setuprun.mat
    dataOut{counter}=trackerout{counter};
    trackerout{counter} = trimmed_tracker;
    counter = counter+1;

    plot(trimmed_tracker(:,1)/1e9, trimmed_tracker(:,2),'LineWidth',2)
xlabel('Time (Gyr)')
ylabel('Temperature (K)')
grid on
  
% 
% plot(1, 'ko', 'MarkerFaceColor','g', 'MarkerSize',10)
% 
% hold on
% 
% %top line is Kepler, bottom line is Earth    
end

% disp('wrapper finished')
% whos

