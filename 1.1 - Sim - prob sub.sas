proc datasets library=work kill nolist;
quit;
proc sql noprint;
	select name into: macro_vars separated by " " from dictionary.macros 
	where scope="GLOBAL" and not substr(name,1,3)="SYS" and not name="GRAPHTERM";
quit;
%symdel macro_vars &macro_vars;
title; footnote;
options nofmterr;

data one0;
	call streaminit(12345);
	do i=1 to 10000;
		x1=rand("Normal", 0, 1);
		x2=rand("Normal", 0, 1);
		x3=rand("Normal", 0, 1);
		x4=rand("Normal", 0, 1);
		x5=rand("Normal", 0, 1);
		x6=rand("Bernoulli", 0.3);
		x7=rand("Bernoulli", 0.4);
		x8=rand("Bernoulli", 0.5);
		x9=rand("Bernoulli", 0.6);
		x10=rand("Bernoulli", 0.7);
		x11=rand("Poisson", 1);
		x12=rand("Poisson", 2);
		x13=rand("Poisson", 3);
		x14=rand("Poisson", 4);
		x15=rand("Poisson", 5);
		x16_raw=rand("Uniform", 0, 1);
		x17_raw=rand("Uniform", 0, 1);
		x18_raw=rand("Uniform", 0, 1);
		x19_raw=rand("Uniform", 0, 1);
		x20_raw=rand("Uniform", 0, 1);

		U=rand("Uniform", 0, 1);
		T_censor=rand("WEIBULL", 2, 2000);
		group=rand("Bernoulli", 0.5);
		split=rand("Uniform", 0, 1);
		output;
	end;
run;

data one1;
	format subject T1_event T0_event T_censor T1 T0 E1 E0 better_effect BEST12.;
	set one0;
	if x16_raw<0.2 then x16=1;
	if x16_raw>=0.2 and x16_raw<0.4 then x16=2;
	if x16_raw>=0.4 and x16_raw<0.85 then x16=3;
	if x16_raw>=0.85 then x16=4;

	if x17_raw<0.4 then x17=1;
	if x17_raw>=0.4 and x17_raw<0.9 then x17=2;
	if x17_raw>=0.9 then x17=3;

	if x18_raw<0.5 then x18=1;
	if x18_raw>=0.5 and x18_raw<0.8 then x18=2;
	if x18_raw>=0.8 then x18=3;

	if x19_raw<0.15 then x19=1;
	if x19_raw>=0.15 and x19_raw<0.4 then x19=2;
	if x19_raw>=0.4 and x19_raw<0.7 then x19=3;
	if x19_raw>=0.7 and x19_raw<0.95 then x19=4;
	if x19_raw>=0.95 then x19=5;

	x20=x20_raw + 0.2*x1 - 0.1*x2;

	trt_effect=0.5;
	better_effect=0.2/(1+exp(-x1)) -0.2*x3 -0.2*x5*x5 -0.2*x7 -0.2*x9 +0.3*(x19=2) +0.5*(x19=3) +1.2*(x19=4 or x17=3);
	prob_better=exp(better_effect)/(1+exp(better_effect));
	call streaminit(12345);
	if prob_better>=0.7 then better_index=rand("Bernoulli", prob_better);
	if prob_better<0.7 then better_index=rand("Bernoulli", prob_better*0.3);
	
	m1=exp(-(-trt_effect) -(0.1/(1+exp(-x2)) +0.2*sin(x4) -0.1*x6*x6 -0.1*x8 -0.2*x10 +0.1*sin(3.14*x11*x12) +0.1*x13));
/*	m0=exp(0 -(-0.7*prob_better) -(0.1/(1+exp(-x2)) +0.2*sin(x4) -0.1*x6*x6 -0.1*x8 -0.2*x10 +0.1*sin(3.14*x11*x12) +0.1*x13));*/
	m0=exp(0 -(-1*better_index) -(0.1/(1+exp(-x2)) +0.2*sin(x4) -0.1*x6*x6 -0.1*x8 -0.2*x10 +0.1*sin(3.14*x11*x12) +0.1*x13));
	yita=2;
	d1=800;
	d0=800;
	T1_event = rand("WEIBULL", yita, d1 * m1);
	T0_event = rand("WEIBULL", yita, d0 * m0);
	T1 = min(T1_event, T_censor);    
	T0 = min(T0_event, T_censor);  
	E1 = (T1_event le T_censor);
	E0 = (T0_event le T_censor);
run;

proc freq data=one1;
	table better_index;
run;

proc format;
	value better_index
		0 = "Not Better"
		1 = "Better"
	;
run;

title;
footnote;
ods _all_ close;
ods listing gpath="C:\Z\1 - Stats\3 - Analyses\Survival Analysis\4 - Simulation\Subgroup\Output\1 - Describe\2 - Better Prob" image_dpi=300 style=mystyle;
ods graphics on / reset=all noscale border=off antialias=on antialiasmax=9800 width=5in height=5in imagename="0-1 - Explore Better Prob";
ods escapechar="^";
proc sgplot data=one1 noautolegend;
	histogram prob_better / name="b" legendlabel="Better Prob" transparency=0.5 binwidth=0.01 scale=count;
	keylegend "a" "b" / location=outside position=topright noborder across=1;
	xaxis display=(nolabel);
run;
ods html;

title;
footnote;
ods _all_ close;
ods listing gpath="C:\Z\1 - Stats\3 - Analyses\Survival Analysis\4 - Simulation\Subgroup\Output\1 - Describe\2 - Better Prob" image_dpi=300 style=mystyle;
ods graphics on / reset=all noscale border=off antialias=on antialiasmax=9800 width=5in height=5in imagename="0-1 - Explore Better Index";
ods escapechar="^";
proc sgplot data=one1 noautolegend;
	histogram prob_better / group=better_index name="b" transparency=0.5 binwidth=0.01 scale=count;
	keylegend "a" "b" / location=outside position=topright noborder across=1;
	format better_index better_index.;
	xaxis display=(nolabel);
run;
ods html;



title;
footnote;
ods _all_ close;
ods listing gpath="C:\Z\1 - Stats\3 - Analyses\Survival Analysis\4 - Simulation\Subgroup\Output\1 - Describe\2 - Better Prob" image_dpi=300 style=mystyle;
ods graphics on / reset=all noscale border=off antialias=on antialiasmax=9800 width=5in height=5in imagename="0-2 - Risk - prob";
ods escapechar="^";
proc sgplot data=one1 noautolegend;
	histogram m1 / name="a" legendlabel="M1" transparency=0.5 binwidth=0.1 scale=count fillattrs=(color=CX3A497E);
	histogram m0 / name="b" legendlabel="M0" transparency=0.5 binwidth=0.1 scale=count;
	keylegend "a" "b" / location=outside position=topright noborder across=1;
	xaxis display=(nolabel) values=(0 to 10);
run;
ods html;

data one2;
	set one1;
	if T1_event>3650 then do;
		E1=0; T1=3650;
	end;
	if T0_event>3650 then do;
		E0=0; T0=3650;
	end;
run;

proc sort data=one2;
	by split;
run;

data one3;
	set one2;
	if _N_<=7000 then train=1;
	else train=0;
run;

proc freq data=one3;
	table E1 E0;
run;

title;
footnote;
ods _all_ close;
ods listing gpath="C:\Z\1 - Stats\3 - Analyses\Survival Analysis\4 - Simulation\Subgroup\Output\1 - Describe\2 - Better Prob" image_dpi=300 style=mystyle;
ods graphics on / reset=all noscale border=off antialias=on antialiasmax=9800 width=5in height=5in imagename="0-3 - Time - prob";
ods escapechar="^";
proc sgplot data=one2 noautolegend;
	histogram T1 / name="a" legendlabel="Time 1" transparency=0.5 scale=count fillattrs=(color=CX3A497E);
	histogram T0 / name="b" legendlabel="Time 0" transparency=0.5 scale=count;
	keylegend "a" "b" / location=outside position=topright noborder across=1;
	where T1<=3650 and T0<3650;
	xaxis label="Time";
run;
ods html;

data sim1(keep=i X1 - X20 better_index prob_better group train Ptime Pstatus T1 T0 E1 E0 rename=(i=subject));
	format i Ptime Pstatus train T1 T0 E1 E0 prob_better better_index x1-x5 x20 BEST12.;
	set one3;
	if group=1 then do; 
		Ptime=T1;
		Pstatus=E1;
	end;
	if group=0 then do; 
		Ptime=T0;
		Pstatus=E0;
	end;
run;

title;
footnote;
ods _all_ close;
ods listing gpath="C:\Z\1 - Stats\3 - Analyses\Survival Analysis\4 - Simulation\Subgroup\Output\1 - Describe\2 - Better Prob" image_dpi=300 style=mystyle;
ods graphics on / reset=all noscale border=off antialias=on antialiasmax=9800 width=5in height=5in imagename="0-4 - KM";
ods escapechar="^";
ods select SurvivalPlot HomTests CensoredSummary;
proc lifetest data=sim1 plots=survival(atrisk=0 to 3650 by 365 test) maxtime=3650;
	time Ptime * Pstatus(0);
	strata group;
run;
ods html;



libname sim "C:\Z\1 - Stats\3 - Analyses\Survival Analysis\4 - Simulation\Subgroup\Data - Analysis";
data sim.probsub1;
	set sim1;
run;




