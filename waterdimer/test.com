%nproc=1
%mem=4GB
%Chk=gradients-qmjob.chk
#P force hf/sto-3g  Charge NoSymm Prop=(Read, Field) punch=derivatives Guess=Read

Generated by Zonkey

0 1
O 1.39518316349 -0.13126868807 0.0012133645013
H 1.57987278312 0.445344558732 -0.766631905859
H 1.5818135125 0.455150160158 0.759419171902

-1.47310078978 0.0803823169263 -0.000737702308408 -0.834
-0.550834870487 0.0423753952678 0.000893453539404 0.417
-1.7769337891 -0.792983742062 0.00384361792166 0.417

-1.47310078978 0.0803823169263 -0.000737702308408 
-0.550834870487 0.0423753952678 0.000893453539404
-1.7769337891 -0.792983742062 0.00384361792166

