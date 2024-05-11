libname schizo "C:\Users\lucp10894\Documents\Hasselt_2020\Thesis\Github";

/*************************************************/
/* Code for replication of analysis in Chapter 3 */
/*************************************************/

/* Case study analysis */
proc import datafile='C:\Users\lucp10894\Documents\Hasselt_2020\Thesis\Github\schizo panss bprs ori.xlsx'
out = schizo_ori dbms=xlsx replace;
run;

/* Transpose the data to long format */
data schizo_long;
set schizo_ori;
endpoint = 1;
outcome = panss;
output;
endpoint = -1;
outcome = bprs;
output;
drop panss bprs;
run;

/* PROC MIXED (UN) */
proc mixed data = schizo_long covtest scoring = 5;
class endpoint id investid;
model outcome = endpoint endpoint*treat / solution noint corrb;
random endpoint endpoint*treat / subject = investid type = un g gcorr;
repeated endpoint / subject = id(investid) type = un r rcorr;
ods output ConvergenceStatus = ConvStatus_un IterHistory = IterHist_un solutionF = solutionF_un CovParms = CovParms_un gcorr = gcorr_un r = r_un rcorr = rcorr_un;
run;

/* PROC MIXED (FA0) */
proc mixed data = schizo_long covtest scoring = 5;
class endpoint id investid;
model outcome = endpoint endpoint*treat / solution noint corrb;
random endpoint endpoint*treat / subject = investid type = fa0(4) g gcorr;
repeated endpoint / subject = id(investid) type = fa0(2) r rcorr;
ods output ConvergenceStatus = ConvStatus_fa IterHistory = IterHist_fa solutionF = solutionF_fa CovParms = CovParms_fa g = g_fa gcorr = gcorr_fa r = r_fa rcorr = rcorr_fa;
run;

/* NOTE: change the name of the data where necessary ("_un" to "_fa") */
/* fixed effects */
data int_S;
set solutionF_un (where = (endpoint = -1 and effect = "endpoint"));
run;
data int_S;
set int_S (keep = estimate);
rename estimate = int_S;
run;

data int_T;
set solutionF_un (where = (endpoint = 1 and effect = "endpoint"));
run;
data int_T;
set int_T (keep = estimate);
rename estimate = int_T;
run;

data treat_S;
set solutionF_un (where = (endpoint = -1 and effect = "Treat*endpoint"));
run;
data treat_S;
set treat_S (keep = estimate);
rename estimate = treat_S;
run;

data treat_T;
set solutionF_un (where = (endpoint = 1 and effect = "Treat*endpoint"));
run;
data treat_T;
set treat_T (keep = estimate);
rename estimate = treat_T;
run;

/* D matrix components */
data D_11;
set CovParms_un (where = (CovParm = "UN(1,1)" and subject = "InvestId"));
run;
data D_11;
set D_11 (keep = estimate);
rename estimate = D_11;
run;

data D_21;
set CovParms_un (where = (CovParm = "UN(2,1)" and subject = "InvestId"));
run;
data D_21;
set D_21 (keep = estimate);
rename estimate = D_21;
run;

data D_22;
set CovParms_un (where = (CovParm = "UN(2,2)" and subject = "InvestId"));
run;
data D_22;
set D_22 (keep = estimate);
rename estimate = D_22;
run;

data D_31;
set CovParms_un (where = (CovParm = "UN(3,1)" and subject = "InvestId"));
run;
data D_31;
set D_31 (keep = estimate);
rename estimate = D_31;
run;

data D_32;
set CovParms_un (where = (CovParm = "UN(3,2)" and subject = "InvestId"));
run;
data D_32;
set D_32 (keep = estimate);
rename estimate = D_32;
run;

data D_33;
set CovParms_un (where = (CovParm = "UN(3,3)" and subject = "InvestId"));
run;
data D_33;
set D_33 (keep = estimate);
rename estimate = D_33;
run;

data D_41;
set CovParms_un (where = (CovParm = "UN(4,1)" and subject = "InvestId"));
run;
data D_41;
set D_41 (keep = estimate);
rename estimate = D_41;
run;

data D_42;
set CovParms_un (where = (CovParm = "UN(4,2)" and subject = "InvestId"));
run;
data D_42;
set D_42 (keep = estimate);
rename estimate = D_42;
run;

data D_43;
set CovParms_un (where = (CovParm = "UN(4,3)" and subject = "InvestId"));
run;
data D_43;
set D_43 (keep = estimate);
rename estimate = D_43;
run;

data D_44;
set CovParms_un (where = (CovParm = "UN(4,4)" and subject = "InvestId"));
run;
data D_44;
set D_44 (keep = estimate);
rename estimate = D_44;
run;

/* sigma matrix components */
data sigma_SS;
set r_un (where = (index = 1 and row = 2));
run;
data sigma_SS;
set sigma_SS (keep = Col2);
rename Col2 = sigma_SS;
run;

data sigma_ST;
set r_un (where = (index = 1 and row = 1));
run;
data sigma_ST;
set sigma_ST (keep = Col2);
rename Col2 = sigma_ST;
run;

data sigma_TT;
set r_un (where = (index = 1 and row = 1));
run;
data sigma_TT;
set sigma_TT (keep = Col1);
rename Col1 = sigma_TT;
run;

/* correlation for residual */
data corr_resid;
set rcorr_un (where = (index = 1 and row = 1));
run;
data corr_resid;
set corr_resid (keep = Col2);
rename Col2 = corr_resid;
run;

/* iteration */
data iter;
set IterHist_un end = last_iter;
if last_iter then output;
keep Iteration;
run;

/* table for all output */
data output;
option validvarname = any;
set ConvStatus_un (keep = Status pdG pdH);
rename Status = SAS_converged pdG = SAS_posdef_G pdH = SAS_posdef_H;
set iter;
set int_S;
set int_T;
set treat_S;
set treat_T;
set D_11;
set D_21;
set D_22;
set D_31;
set D_32;
set D_33;
set D_41;
set D_42;
set D_43;
set D_44;
set sigma_SS;
set sigma_ST;
set sigma_TT;
set corr_resid;
run;

/* export data */
proc export 
data = work.output
outfile = "C:\Users\lucp10894\Documents\Hasselt_2020\Thesis\Github\Schizo_ori_pb UN.txt"
DBMS = TAB replace;
putnames = yes;
run;

/* %contcontfull: two-stage analysis */

proc import datafile='C:\Users\lucp10894\Documents\Hasselt_2020\Thesis\Github\schizo panss bprs ori.xlsx'
out = schizo.schizo_pb dbms=xlsx replace;
run;

%inc "C:\Users\lucp10894\Documents\Hasselt_2020\Thesis\Github\CONTCONTFULL.sas";
%contcontfull(data=schizo.schizo_pb, true=panss, surrog=bprs, trt=treat, trial=investid, patid=id, weighted=1);


