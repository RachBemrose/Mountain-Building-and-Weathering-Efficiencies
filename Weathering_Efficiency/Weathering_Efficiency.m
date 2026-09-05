function Bemrose_Weathering_Rates = Weathering_Efficiency(~)

Radius = 1.02;

R_p = [0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8];

Computed_Weathering_Rates = [0.211,0.2337,0.29247,0.40274,0.67657,1,1.2737,1.6057,2.2082,3.4651,5.38,9.0118,14.954,22.45];

Bemrose_Weathering_Rates = interp1(R_p,Computed_Weathering_Rates,Radius,'spline');

end

% Tau Ceti E = 6.1467 when radius is 1.527
% Kepler 442b = 2.6385 when radius is 1.34
%HD 20794d = 40.422 when radius is 2.04
% Trappist 1E = 0.74252 when radius is 0.92
% Kepler 452b = 10.568 when radius is 1.63
% Proxima Centuri B = 1.0574 when radius is 1.02
% Wolf1069b = 1.2189 when radius is 1.08