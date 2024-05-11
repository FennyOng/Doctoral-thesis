*************************************************************************************************;
* Program          : CONTCONTFULL.SAS                                                     	 	                                        
                                              													
* Description      : Macro perfomrs surrogate validation analysis for two normal endpoints		 
*					 using Full fixed effects model as described in the book
					" Applied Surrogate Endpoint Evaluation Methods with SAS and R (Chapter 4 & 12)".              		
* SAS version      : Built and tested on SAS 9.4                                          		
* Programmer       : Theophile Bigirumurame                                  					

* Arguments:
* Data: a dataset containing one record per patient with measurements for both the true 		
  and the surrogate endpoints.
* true: a measurement of the true endpoint.  
* surrog: a measurement of the surrogate endpoint. 
* trt: treatment indicator variable (1= active, -1=control ). 
* trial: trial (or center) in which the patient was treated. This is the unit of the 
	study for which trial level surrogacy will be estimated.
* patid: patient's identification number.
* weighted: an option which allows to use weighted regression (weighted=1) in the 
	computation of the trial level surrogacy. The number of patients in the trial is used 
	as the weight.
* looa: an option to perform a leave one trial out analysis (looa=1). Both
	surrogacy measures are computed by leaving each one of the trial out in order to 
	check the influence of the left out trial on the overall surrogacy measures.

*******************************************************************************************;

%macro CONTCONTFULL(data=,
					true=,
					surrog=,
					trt=,
					trial=,
					patid=,
					weighted=,
					looa=,
					imagefmt=		/* Specify the format of plots(image) outputs . Possibilities are: PNG (default), BMP, DIB, EMF, EPSI, GIF, JFIF, JPEG, 
		                                     JPG, PBM, PDF, PS, SASEMF, STATIC, TIFF, WMF, XBM, XPM, PSL, SVG. For the description of these output formats, go to
		                                     http://support.sas.com/documentation/cdl/en/odsug/65308/HTML/default/viewer.htm#p0kroq43yu0lspn16hk1u4c65lti.htm        */ 
);


					%if %sysfunc(exist(&data))=1 %then %do;
		%put SPECIFIED DATASET IS :&data ;
%end;

%else %do;
		%put ERROR: WRONG OR NO DATASET IS SPECIFIED ;
		%put ERROR: PLEASE SPECIFY THE ARGUMENT [DATA] CORRECTLY;
		%goto endmac;
%end;


%if %sysfunc(varnum(%sysfunc(open(&data)),&true))=0 %then %do;
		%put ERROR: WRONG OR NO TRUE ENDPOINT IS SPECIFIED ;
		%put ERROR: PLEASE SPECIFY THE ARGUMENT [TRUE] CORRECTLY;
		%goto endmac;
%end;

%else %do;
		%put SPECIFIED TRUE ENDPOINT IS :&true ;
		
%end;



%if %sysfunc(varnum(%sysfunc(open(&data)),&surrog))=0 %then %do;
		%put ERROR: WRONG OR NO SURROGATE ENDPOINT IS SPECIFIED ;
		%put ERROR: PLEASE SPECIFY THE ARGUMENT [SURROG] CORRECTLY;
		%goto endmac;
%end;

%else %do;
		%put SPECIFIED TRUE ENDPOINT IS :&surrog ;
		
%end;


%if %sysfunc(varnum(%sysfunc(open(&data)),&trt))=0 %then %do;
		%put ERROR: WRONG OR NO TREATMENT IS SPECIFIED ;
		%put ERROR: PLEASE SPECIFY THE ARGUMENT [TRT] CORRECTLY;
		%goto endmac;
%end;

%else %do;
		%put SPECIFIED TREATMENT IS :&trt ;
		
%end;

%if %sysfunc(varnum(%sysfunc(open(&data)),&trial))=0 %then %do;
		%put ERROR: WRONG OR NO TRIAL IS SPECIFIED ;
		%put ERROR: PLEASE SPECIFY THE ARGUMENT [TRIAL] CORRECTLY;
		%goto endmac;
%end;

%else %do;
		%put SPECIFIED TRIAL IS :&trial ;
		
%end;


%if %sysfunc(varnum(%sysfunc(open(&data)),&patid))=0 %then %do;
		%put ERROR: WRONG OR NO PATIENT ID IS SPECIFIED ;
		%put ERROR: PLEASE SPECIFY THE ARGUMENT [PATIENT ID] CORRECTLY;
		%goto endmac;
%end;

%else %do;
		%put SPECIFIED PATIENT ID IS :&patid ;
		
%end;

data dataxx;
set &data(rename=(&true=truexx &surrog=surrogxx &trt=trtxx &trial=trialxx &patid=patidxx));
run;

data norm;
set dataxx;
response = truexx; 
endp   = 1;
output;
response = surrogxx;
endp    = -1;
output;
keep response endp patidxx trtxx  trialxx truexx surrogxx ;
run;
quit;

proc sql ;/*the weight data set containig the number of &subjects per&trial*/
create table weight2 as
select trialxx as trial,trtxx as treatment, count(patidxx)as n
from dataxx
group by trialxx,trtxx;
quit;

ods listing ;
ods graphics on/ reset=index  border=off imagename="expfullfix1" outputfmt=&imagefmt;
proc sgplot data=weight2 pad=(bottom=5%);                                                                                                       
vbar trial / response=n group=treatment  nostatlabel
groupdisplay=cluster;
xaxis label="Trials";
yaxis  label="Frequency";
title 'Number of subject per trial and per treatment arm ' ;
run; 
title;
ods graphics off;
ods listing close; 
quit;


proc sql ;
create table weight as
select trialxx , count(patidxx) as n, monotonic() as num
from dataxx
group by trialxx;
quit;



data norm;
merge norm weight;
by trialxx;
run;
quit;


ods listing ;
ods graphics on/ reset=index  border=off imagename="expfullfix2" outputfmt=&imagefmt;
proc sgplot data=dataxx;
scatter  x=surrogxx y=truexx;
title 'True endpoint vs. surrogate endpoint scatter plot';
xaxis label="Surrogate endpoint ";
yaxis label= "True endpoint ";
run;
title;
ods graphics off;
ods listing close;
quit;


data norm;
set norm;
if trtxx=1 then trtxx1=1; else trtxx1=-1;
outcome=response;
run;


ods select none;
PROC MIXED DATA=norm COVTEST;
CLASS endp patidxx trialxx;
MODEL outcome = endp*trialxx endp*trtxx1*trialxx / S NOINT cl outp=predic;
REPEATED endp / TYPE=un SUBJECT=patidxx(trialxx) ;
ods output solutionF=eb CovParms=covar ConvergenceStatus=convergence;
RUN;

*** INDIVIDUAL LEVEL SURROGACY  ***;
proc iml;
reset log;
use covar;
read all var{ESTIMATE} where(COVPARM="UN(2,1)") into d12;
read all var{ESTIMATE} where(COVPARM="UN(2,2)") into d22;
read all var{ESTIMATE} where(COVPARM="UN(1,1)") into d11;
close covar;
indiv=d12##2/(d11#d22);
create indiv var {indiv};
append;
quit;

data ebnu;
set eb;
if missing(trialxx) then delete;
run; 

 proc sql;
 create table true1 as
 select trialxx , Lower as Lower_t, Estimate as True, Upper as Upper_t 
 from ebnu
 where Effect="trtxx1*endp*trialxx" & endp=1  ;
 quit;

proc sql;
create table surr1 as
select trialxx, Lower as Lower_s, Estimate as Surr, Upper as Upper_s
from ebnu
where Effect="trtxx1*endp*trialxx" & endp=-1   ;
quit; 

proc sql;
create table surr12 as
select trialxx ,Lower as Lower_int, Estimate as Surrinterc, Upper as Upper_int
from ebnu
where Effect="endp*trialxx" & endp=-1  ;
quit;

data in;
merge true1 surr1 surr12;
by trialxx;
run;


***REMOVED****;
ods select all;
proc report data = in  nowd headline headskip split='$' 
style(report)=[cellspacing=1 borderwidth=1  bordercolor=black 
font_face="arial,helvetica"] style(header)={background=lightskyblue foreground=black};
column  trial ("True " Lower_t True Upper_t) ("Surrogate" Lower_s Surr Upper_s) ("Intercept" Lower_int Surrinterc Upper_int);
define trial / display "Trial" Style(column)=[ just=left  font_size=3] order=formatted width=96;
define Lower_t / display "Lower" Style(column)=[ just=left  font_size=3] order=formatted width=96;
define True / display "Estimate" Style(column)=[ just=left  font_size=3] order=formatted width=96;
define Upper_t / display "Upper" Style(column)=[ just=left  font_size=3] order=formatted width=96;
define Lower_s / display "Lower" Style(column)=[ just=left  font_size=3] order=formatted width=96;
define Surr / display "Estimate" Style(column)=[ just=left  font_size=3] order=formatted width=96;
define Upper_s / display "Upper" Style(column)=[ just=left  font_size=3] order=formatted width=96; 
define Lower_int / display "Lower" Style(column)=[ just=left  font_size=3] order=formatted width=96;
define Surrinterc / display "Estimate" Style(column)=[ just=left  font_size=3] order=formatted width=96;
define Upper_int / display "Upper" Style(column)=[ just=left  font_size=3] order=formatted width=96;
title 'Parameter estimates and confidence intervals for both endpoints';
format Lower_t True Upper_t Lower_s Surr Upper_s Lower_int Surrinterc Upper_int: comma9.4;
run;
title;


****TRIAL LEVEL SURROGACY  *****;
proc sql noprint;/*NUMBER of trial*/
select  count(trialxx) into :ntrial
from (select distinct trialxx from dataxx)
quit;

proc sql noprint ;/*NUMBER of subject*/
select  count(patidxx) into :ntotal
from  dataxx
quit;

data in2;
merge true1 surr1 surr12 weight;
by trialxx;
run; 
/* Calculating unadjusted R2.*/


%if &weighted=1 %then %do;
    ods select none;
	proc reg data=in2;
	model true=surrinterc surr;
	weight n;
	ods output FitStatistics=rsvv;
	run;
	quit;

    ods select all;
	proc iml;
	reset log;
	use rsvv;
	read all var{nValue2} where(Label2="R-Square") into R2;
	close rsvv;
	create trialv var {R2};
	append;
	quit;

    data mes;
    merge trialv indiv;
    run;

	data _null_;
	set mes;
	call symput('indiv',put(indiv,6.4));
	call symput('triall',put(R2,6.4));
	run;

    data surr_measure;
	set mes;
	indiv_sd=sqrt((4*indiv*(1-indiv)**2)/(&ntotal-3));
	LO_indiv=max(0, indiv + probit(0.025) *(indiv_sd));
	UP_indiv=min(1, indiv + probit(0.975)*(indiv_sd));
	R2_sd=sqrt((4*R2*(1-R2)**2)/(&ntrial-3));
	LOO_R2=max(0, R2 + probit(0.025) *(R2_sd));
	UP_R2=min(1, R2 + probit(0.975)*(R2_sd));
	format indiv R2 LOO_R2 UP_R2 LO_indiv UP_indiv: comma9.4;
	run;


	ods select all;
	options nodate pageno=1 linesize=64 pagesize=60 fmtsearch=(proclib);
	proc report data=surr_measure nowd headskip headline split='*'
	style(header)={background=lightskyblue foreground=black};
   	column  ("INDIVIDUAL" LO_indiv indiv UP_indiv)("TRIAL" LOO_R2 R2 UP_R2) ;
	define LO_indiv / display   'LOWER';
	define indiv / display 'Individual' order;
   	define UP_indiv / display 'UPPER '  ;
   	define LOO_R2 / display   'LOWER';
	define R2/ display 'R square' order;
   	define UP_R2 / display 'UPPER '  ;
   	format LO_indiv indiv UP_indiv LOO_R2 R2 UP_R2: comma9.4;
   	title "Surrogacy measures";
   	footnote 'Weighted regression was used ';
	run;
	title;

%end;

%if &weighted=0 %then %do;
	ods select none;
	proc reg data=in2;
	model true=surrinterc surr;
	ods output FitStatistics=rsvv;
	run;
	quit;


	ods select all;
	proc iml;
	reset log;
	use rsvv;
	read all var{nValue2} where(Label2="R-Square") into R2;
	close rsvv;
	create trialv var {R2};
	append;
	quit;

	data mes;
	merge trialv indiv;
	run;

	data _null_;
	set mes;
	call symput('indiv',put(indiv,6.4));
	call symput('triall',put(R2,6.4));
	run;

	data surr_measure;
	set mes;
	indiv_sd=sqrt((4*indiv*(1-indiv)**2)/(&ntotal-3));
	LO_indiv=max(0, indiv + probit(0.025) *(indiv_sd));
	UP_indiv=min(1, indiv + probit(0.975)*(indiv_sd));
	R2_sd=sqrt((4*R2*(1-R2)**2)/(&ntrial-3));
	LOO_R2=max(0, R2 + probit(0.025) *(R2_sd));
	UP_R2=min(1, R2 + probit(0.975)*(R2_sd));
	format indiv R2 LOO_R2 UP_R2 LO_indiv UP_indiv: comma9.4;
	run;


	ods select all;
	options nodate pageno=1 linesize=64 pagesize=60 fmtsearch=(proclib);
	proc report data=surr_measure nowd headskip headline split='*'
	style(header)={background=lightskyblue foreground=black};
    column  ("INDIVIDUAL" LO_indiv indiv UP_indiv)("TRIAL" LOO_R2 R2 UP_R2) ;
	define LO_indiv / display   'LOWER';
	define indiv / display 'Individual' order;
   	define UP_indiv / display 'UPPER '  ;
   	define LOO_R2 / display   'LOWER';
	define R2/ display 'R square' order;
	define UP_R2 / display 'UPPER '  ;
	format LO_indiv indiv UP_indiv LOO_R2 R2 UP_R2: comma9.4;
	title "Surrogacy measures";
	footnote "Unweighted regression was used";
	run;
	title;
%end;


%if &looa=1 %then %do;
%macro leaveone;
	 proc sql noprint;
	 select  count(trialxx) into :Nt
	 from (select distinct trialxx from norm)
	 quit;

 %do i = 1 %to &Nt;
	  	Data analysis holdout; 
	  	set norm;
	  	if num = &i then output holdout;
	  	else output analysis;
	  	run;

	/* REDUCED FIXED EFFECT*/;
		ods select none;
		PROC MIXED DATA=analysis COVTEST;
		CLASS endp patidxx trialxx;
		MODEL outcome = endp*trialxx endp*trtxx1*trialxx / S NOINT;
		REPEATED endp / TYPE=un SUBJECT=patidxx(trialxx) ;
		ods output solutionF=eb1 CovParms=covar1;
		RUN;

		*** INDIVIDUAL LEVEL SURROGACY  ***;
		proc iml;
		reset log;
		use covar1;
		read all var{ESTIMATE} where(COVPARM="UN(2,1)") into d12;
		read all var{ESTIMATE} where(COVPARM="UN(2,2)") into d22;
		read all var{ESTIMATE} where(COVPARM="UN(1,1)") into d11;
		close covar1;
		indiv=d12##2/(d11#d22);
		create indiv1 var {indiv};
		append;
		quit;

		data indiv1;
		set indiv1;
		num=&i;
		run;

		proc append base=rsquare12 data=indiv1;
		run;

		/* Trial level */;

		data surr3(rename=( estimate=surrog stderr=surrogstdr)) true3(rename=( estimate=true stderr=truestdr))
		surro3(rename=( estimate=surroginterc stderr=surrogintercstdr));
		set eb1;
		if missing(trialxx) then delete;
		keep trialxx endp  estimate stderr;
		if Effect="trtxx1*endp*trialxx" & endp=-1 then output surr3 ;
		if Effect="trtxx1*endp*trialxx" & endp=1 then output true3 ;
		if Effect="endp*trialxx" & endp=-1 then output surro3 ;
		run;

		proc sql ;/*the weight data set containig the number of &subjects per&trial*/;
		create table weight2 as
		select trialxx as trial, count(patidxx)/2 as n
		from analysis
		group by trialxx;
		quit;


		data est2;
		merge surr3 true3 surro3 weight2;
		keep trialxx surrog true surrogstdr truestdr surroginterc surrogintercstdr id n;
		id=_n_;
		run;

	
	%if &weighted=1 %then %do;
		ods select none;
		proc reg data=est2;
		model true=surroginterc surrog;
		weight n;
		ods output FitStatistics=rsv;
		run;
		quit;
	%end;

	%if &weighted=0 %then %do;
		ods select none;
		proc reg data=est2;
		model true=surroginterc surrog;
		ods output FitStatistics=rsv;
		run;
		quit;
	%end;

		ods select all;
		proc iml;
		reset log;
		use rsv;
		read all var{nValue2} where(Label2="Adj R-Sq") into R2;
		close rsv;
		create trial2 var {R2};
		append;
		quit;

			
		data trial2;
		set trial2;
		num=&i;
		run;

		proc append base=rsqu data=trial2;
		run;
%end;

		data norm1;
		set norm;
		keep trialxx num;
		attrib _all_ label='';
		run;

		proc sort data=norm1 nodupkey;
		by num;
		run;

		data rsquare21;
		merge rsquare12 rsqu norm1;
		by num;
		drop num;
		run;

		ods select all;
		proc report data = rsquare21  nowd headline headskip split='$' 
		style(report)=[cellspacing=1 borderwidth=1  bordercolor=black 
		font_face="arial,helvetica"] style(header)={background=lightskyblue foreground=black};
	    column   trialxx  indiv R2 ;  
	    define trialxx / display "Removed trial" Style(column)=[ just=left  font_size=3] order=formatted width=96;
	    define indiv / display "Indiv. level" Style(column)=[ just=left  font_size=3] order=formatted width=50;
	    define R2 / display "Trial level" Style(column)=[ just=left  font_size=3] order=formatted width=96;
	    title 'Leave one out analysis ';
		format INDIV R2: comma9.4;
		run;
		title;


		ods select all;
		ods listing  ;
		ods graphics on/ reset=index border=off imagename="expfullfix3" outputfmt=&imagefmt;
		proc sgplot data=rsquare21 noautolegend;
		scatter x=trialxx y=R2 /DATALABELATTRS=(Color=Green Family=Arial Size=10 Style=Italic Weight=Bold) ;
		refline &triall/axis=Y ;
		yaxis label="R(*ESC*){unicode '00b2'x}" max=1 min=0;
		xaxis TYPE  = DISCRETE label="Removed trial ";
		title 'Leave one out plot (trial level)';
		run;
		title;
		ods graphics off;
		ods listing close;
		quit;

		ods select all;
		ods listing  ;
		ods graphics on/ reset=index  border=off imagename="expfullfix4" outputfmt=&imagefmt;
		proc sgplot data=rsquare21 noautolegend;
		scatter x=trialxx y=indiv /DATALABELATTRS=(Color=Green Family=Arial Size=10 Style=Italic Weight=Bold) ;
		refline &indiv/axis=Y ;
		yaxis label="R(*ESC*){unicode '00b2'x}" max=1 min=0;
		xaxis type  = DISCRETE  label="Removed trial"; 
		title 'Leave one out plot (individual level)';
		run;
		title;
		ods graphics off;
		ods listing close;
		quit;
		
%mend;
%leaveone;
%end;
%endmac:
%mend;
