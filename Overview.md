# Guide to scripts
Supporting information to Jørgensen et al. (2021) "A unifying model to estimate thermal tolerance limits across static, dynamic and fluctuating exposures to thermal stress in ectotherms". https://www.researchsquare.com/article/rs-347544/v1. Please cite appropriately.
Authors: Lisa Bjerregaard Jørgensen, Hans Malte, Michael Ørsted, Nikolaj Andreasen Klahn & Johannes Overgaard

______

As a practical application of the mathematical framework presented in the paper, we here provide R-scripts to derive parameters of the Thermal Death Time (TDT) curve and use these to assess thermal tolerance limits. Here, we present two scripts and a guide to these (READ_ME.pdf). Please make sure to read through the READ_ME file before using the scripts. Which script you should use depend on the type of input data which depends on the type of experiment conducted:

1)	“TDT_from_Static.R”. This script derives TDT parameters from static experiments, where time to failure (tcoma) is measured at one or more constant temperatures. An input data template is provided (static_input.csv)

2)	“TDT_from_Dynamic.R”. This script derives TDT parameters from dynamic experiments, where the maximal temperatures tolerated (dynamic CTmax, dCTmax) are measured using one or more ramping rates. An input data template is provided (dynamic_input.csv). 

Both scripts use the derived TDT parameters to convert between and within static and dynamic measurements and input data of natural temperature fluctuations can be added to assess when failure occurs based on the derived TTL parameters. We supply templates for input data, with headers that correspond to the respective scripts, please see the READ_ME guide for appropriate use of these.
