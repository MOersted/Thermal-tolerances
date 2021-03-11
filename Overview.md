# Guide to scripts
Supporting information to the publication:
"A unifying model to estimate thermal tolerance limits across static, dynamic and fluctuating exposures to thermal stress in ectotherms"
Authors: Lisa Bjerregaard Jørgensen, Hans Malte, Michael Ørsted, Nikolaj Andreasen Klahn & Johannes Overgaard

______
As a practical application of the mathematical framework presented in the paper, we here provide R-scripts to derive parameters of the Thermal Tolerance Landscape (TTL) and use these to assess thermal tolerance limits. Here, we present two scripts and a guide to these (READ_ME.pdf). Which script you should use depend on the type of input data which depends on the type of experiment conducted:

1)	“TTL_from_Static.R”. This script derives TTL parameters from static experiments, where time to failure (tcoma) is measured at one or more constant temperatures. An input data template is provided (static_input.csv)

2)	“TTL_from_Dynamic.R”. This script derives TTL parameters from dynamic experiments, where the maximal temperatures tolerated (dynamic CTmax, dCTmax) are measured using one or more ramping rates. An input data template is provided (dynamic_input.csv). 

Both scripts use the derived TTL parameters to convert between and within static and dynamic measurements and input data of natural temperature fluctuations can be added to assess when failure occurs based on the derived TTL parameters.
