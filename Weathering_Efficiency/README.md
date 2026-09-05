This script was developed as a function derived from the original computed values of W(R) in the methodology and results. The script calls the relevant function in 'updated_compute.m', whilst using a spline interpolation to compute a weathering efficiency value for any given radius, with no other data required.

## Interpretations
R_p represents planetary radius.
Computed_weathering_rates represents the original values calculated with the derived W(R) equation as well as the values predicted from the look up table and correlational fit.
Bemrose_weathering_rates takes the above values and generates a spline interpolation that is then subsequently used as the Weathering Efficiency function.
