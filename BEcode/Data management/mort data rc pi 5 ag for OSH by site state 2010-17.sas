libname mydir "\\cdc\project\NCIRD_ID_EPI_Branch\ARM\Disease Burden\NCHS Mortality 2010-2017" ;
/*CREATE DATE (OF DEATH) VARIABLE*/
data a1; 
set mydir.mort2010 
mydir.mort2011
mydir.mort2012
mydir.mort2013
mydir.mort2014
mydir.mort2015
mydir.mort2016
mydir.mort2017
;
run;
***************************************************************************;
***************************************************************************;
data a1;
set a1;

%let monthendcode=10,11,12;
%let monthbegincode=1,2,3,4;

if (monthdth in (&monthendcode) AND year=2010) OR (monthdth in (&monthbegincode) AND year=2011) then season=1;
if (monthdth in (&monthendcode) AND year=2011) OR (monthdth in (&monthbegincode) AND year=2012) then season=2;
if (monthdth in (&monthendcode) AND year=2012) OR (monthdth in (&monthbegincode) AND year=2013) then season=3;
if (monthdth in (&monthendcode) AND year=2013) OR (monthdth in (&monthbegincode) AND year=2014) then season=4;
if (monthdth in (&monthendcode) AND year=2014) OR (monthdth in (&monthbegincode) AND year=2015) then season=5;
if (monthdth in (&monthendcode) AND year=2015) OR (monthdth in (&monthbegincode) AND year=2016) then season=6;
if (monthdth in (&monthendcode) AND year=2016) OR (monthdth in (&monthbegincode) AND year=2017) then season=7;

if season ne .;

/*DEAL WITH AGE VARIABLE*/
if age <2000 then age = age -1000;
if 2000<= age <7000 then age = 0;
if age =9999 then age = .;
* Define causes of death:

* ICD10 for influenza codes, all;
%let flutotcode='J09','J10','J11';

* ICD10 for respiratory system codes;
%let rcode='J00','J01','J02','J33','J04','J05','J06','J07','J08','J09','J10',
			'J11','J12','J13','J14','J15','J16','J17','J18','J19','J20',
			'J21','J22','J23','J24','J25','J26','J27','J28','J29','J30',
			'J31','J32','J33','J34','J35','J36','J37','J38','J39','J40',
			'J41','J42','J43','J44','J45','J46','J47','J48','J49','J50',
			'J51','J52','J53','J54','J55','J56','J57','J58','J59','J60',
			'J61','J62','J63','J64','J65','J66','J67','J68','J69','J70',
			'J71','J72','J73','J74','J75','J76','J77','J78','J79','J80',
			'J81','J82','J83','J84','J85','J86','J87','J88','J89','J90',
			'J91','J92','J93','J94','J95','J96','J97','J98','J99';

* ICD10 for circulatory system codes;
%let ccode='I00','I01','I02','I33','I04','I05','I06','I07','I08','I09','I10',
			'I11','I12','I13','I14','I15','I16','I17','I18','I19','I20',
			'I21','I22','I23','I24','I25','I26','I27','I28','I29','I30',
			'I31','I32','I33','I34','I35','I36','I37','I38','I39','I40',
			'I41','I42','I43','I44','I45','I46','I47','I48','I49','I50',
			'I51','I52','I53','I54','I55','I56','I57','I58','I59','I60',
			'I61','I62','I63','I64','I65','I66','I67','I68','I69','I70',
			'I71','I72','I73','I74','I75','I76','I77','I78','I79','I80',
			'I81','I82','I83','I84','I85','I86','I87','I88','I89','I90',
			'I91','I92','I93','I94','I95','I96','I97','I98','I99';

* ICD10 for pneumonia;
%let pcode='J12','J13','J14','J15','J16','J17','J18';

* ICD10 for Flu codes;
%let flutotcode='J10','J11';


* PULL DATA FROM US MORTALITY FOR AGED 15-49 (AGER27 9-15)   ;
* ALL "O" CODE, RESPIRATORY AND P&I CODES (J00-J99) (J09-J18);
if  substr(ucod,1,3) in (&rcode) or
	substr(enicon_1,1,3) in (&rcode) or substr(enicon_2,1,3) in (&rcode) or substr(enicon_3,1,3) in (&rcode) or
	substr(enicon_4,1,3) in (&rcode) or substr(enicon_5,1,3) in (&rcode) or substr(enicon_6,1,3) in (&rcode) or
	substr(enicon_7,1,3) in (&rcode) or substr(enicon_8,1,3) in (&rcode) or substr(enicon_9,1,3) in (&rcode) or
	substr(enicon_10,1,3) in (&rcode) or substr(enicon_11,1,3) in (&rcode) or substr(enicon_12,1,3) in (&rcode) or
	substr(enicon_13,1,3) in (&rcode) or substr(enicon_14,1,3) in (&rcode) or substr(enicon_15,1,3) in (&rcode) or
	substr(enicon_16,1,3) in (&rcode) or substr(enicon_17,1,3) in (&rcode) or substr(enicon_18,1,3) in (&rcode) or
	substr(enicon_19,1,3) in (&rcode) or substr(enicon_20,1,3) in (&rcode) or 
	substr(record_1,1,3) in (&rcode) or substr(record_2,1,3) in (&rcode) or substr(record_3,1,3) in (&rcode) or
	substr(record_4,1,3) in (&rcode) or substr(record_5,1,3) in (&rcode) or substr(record_6,1,3) in (&rcode) or
	substr(record_7,1,3) in (&rcode) or substr(record_8,1,3) in (&rcode) or substr(record_9,1,3) in (&rcode) or
	substr(record_10,1,3) in (&rcode) or substr(record_11,1,3) in (&rcode) or substr(record_12,1,3) in (&rcode) or
	substr(record_13,1,3) in (&rcode) or substr(record_14,1,3) in (&rcode) or substr(record_15,1,3) in (&rcode) or
	substr(record_16,1,3) in (&rcode) or substr(record_17,1,3) in (&rcode) or substr(record_18,1,3) in (&rcode) or
	substr(record_19,1,3) in (&rcode) or substr(record_20,1,3) in (&rcode) or

	substr(ucod,1,3) in (&ccode) or
	substr(enicon_1,1,3) in (&ccode) or substr(enicon_2,1,3) in (&ccode) or substr(enicon_3,1,3) in (&ccode) or
	substr(enicon_4,1,3) in (&ccode) or substr(enicon_5,1,3) in (&ccode) or substr(enicon_6,1,3) in (&ccode) or
	substr(enicon_7,1,3) in (&ccode) or substr(enicon_8,1,3) in (&ccode) or substr(enicon_9,1,3) in (&ccode) or
	substr(enicon_10,1,3) in (&ccode) or substr(enicon_11,1,3) in (&ccode) or substr(enicon_12,1,3) in (&ccode) or
	substr(enicon_13,1,3) in (&ccode) or substr(enicon_14,1,3) in (&ccode) or substr(enicon_15,1,3) in (&ccode) or
	substr(enicon_16,1,3) in (&ccode) or substr(enicon_17,1,3) in (&ccode) or substr(enicon_18,1,3) in (&ccode) or
	substr(enicon_19,1,3) in (&ccode) or substr(enicon_20,1,3) in (&ccode) or 
	substr(record_1,1,3) in (&ccode) or substr(record_2,1,3) in (&ccode) or substr(record_3,1,3) in (&ccode) or
	substr(record_4,1,3) in (&ccode) or substr(record_5,1,3) in (&ccode) or substr(record_6,1,3) in (&ccode) or
	substr(record_7,1,3) in (&ccode) or substr(record_8,1,3) in (&ccode) or substr(record_9,1,3) in (&ccode) or
	substr(record_10,1,3) in (&ccode) or substr(record_11,1,3) in (&ccode) or substr(record_12,1,3) in (&ccode) or
	substr(record_13,1,3) in (&ccode) or substr(record_14,1,3) in (&ccode) or substr(record_15,1,3) in (&ccode) or
	substr(record_16,1,3) in (&ccode) or substr(record_17,1,3) in (&ccode) or substr(record_18,1,3) in (&ccode) or
	substr(record_19,1,3) in (&ccode) or substr(record_20,1,3) in (&ccode) then rcu=1;
	else rcu=0;

if  substr(ucod,1,3) in (&pcode) or
	substr(enicon_1,1,3) in (&pcode) or substr(enicon_2,1,3) in (&pcode) or substr(enicon_3,1,3) in (&pcode) or
	substr(enicon_4,1,3) in (&pcode) or substr(enicon_5,1,3) in (&pcode) or substr(enicon_6,1,3) in (&pcode) or
	substr(enicon_7,1,3) in (&pcode) or substr(enicon_8,1,3) in (&pcode) or substr(enicon_9,1,3) in (&pcode) or
	substr(enicon_10,1,3) in (&pcode) or substr(enicon_11,1,3) in (&pcode) or substr(enicon_12,1,3) in (&pcode) or
	substr(enicon_13,1,3) in (&pcode) or substr(enicon_14,1,3) in (&pcode) or substr(enicon_15,1,3) in (&pcode) or
	substr(enicon_16,1,3) in (&pcode) or substr(enicon_17,1,3) in (&pcode) or substr(enicon_18,1,3) in (&pcode) or
	substr(enicon_19,1,3) in (&pcode) or substr(enicon_20,1,3) in (&pcode) or 
	substr(record_1,1,3) in (&pcode) or substr(record_2,1,3) in (&pcode) or substr(record_3,1,3) in (&pcode) or
	substr(record_4,1,3) in (&pcode) or substr(record_5,1,3) in (&pcode) or substr(record_6,1,3) in (&pcode) or
	substr(record_7,1,3) in (&pcode) or substr(record_8,1,3) in (&pcode) or substr(record_9,1,3) in (&pcode) or
	substr(record_10,1,3) in (&pcode) or substr(record_11,1,3) in (&pcode) or substr(record_12,1,3) in (&pcode) or
	substr(record_13,1,3) in (&pcode) or substr(record_14,1,3) in (&pcode) or substr(record_15,1,3) in (&pcode) or
	substr(record_16,1,3) in (&pcode) or substr(record_17,1,3) in (&pcode) or substr(record_18,1,3) in (&pcode) or
	substr(record_19,1,3) in (&pcode) or substr(record_20,1,3) in (&pcode) or

	substr(ucod,1,3) in (&flutotcode) or
	substr(enicon_1,1,3) in (&flutotcode) or substr(enicon_2,1,3) in (&flutotcode) or substr(enicon_3,1,3) in (&flutotcode) or
	substr(enicon_4,1,3) in (&flutotcode) or substr(enicon_5,1,3) in (&flutotcode) or substr(enicon_6,1,3) in (&flutotcode) or
	substr(enicon_7,1,3) in (&flutotcode) or substr(enicon_8,1,3) in (&flutotcode) or substr(enicon_9,1,3) in (&flutotcode) or
	substr(enicon_10,1,3) in (&flutotcode) or substr(enicon_11,1,3) in (&flutotcode) or substr(enicon_12,1,3) in (&flutotcode) or
	substr(enicon_13,1,3) in (&flutotcode) or substr(enicon_14,1,3) in (&flutotcode) or substr(enicon_15,1,3) in (&flutotcode) or
	substr(enicon_16,1,3) in (&flutotcode) or substr(enicon_17,1,3) in (&flutotcode) or substr(enicon_18,1,3) in (&flutotcode) or
	substr(enicon_19,1,3) in (&flutotcode) or substr(enicon_20,1,3) in (&flutotcode) or 
	substr(record_1,1,3) in (&flutotcode) or substr(record_2,1,3) in (&flutotcode) or substr(record_3,1,3) in (&flutotcode) or
	substr(record_4,1,3) in (&flutotcode) or substr(record_5,1,3) in (&flutotcode) or substr(record_6,1,3) in (&flutotcode) or
	substr(record_7,1,3) in (&flutotcode) or substr(record_8,1,3) in (&flutotcode) or substr(record_9,1,3) in (&flutotcode) or
	substr(record_10,1,3) in (&flutotcode) or substr(record_11,1,3) in (&flutotcode) or substr(record_12,1,3) in (&flutotcode) or
	substr(record_13,1,3) in (&flutotcode) or substr(record_14,1,3) in (&flutotcode) or substr(record_15,1,3) in (&flutotcode) or
	substr(record_16,1,3) in (&flutotcode) or substr(record_17,1,3) in (&flutotcode) or substr(record_18,1,3) in (&flutotcode) or
	substr(record_19,1,3) in (&flutotcode) or substr(record_20,1,3) in (&flutotcode) then pi=1;
	else pi=0;

** Restrict to counties in FSN;

if age <=4 then agecat=1;
if 5<= age <=17 then agecat=2;
if 18<= age <=49 then agecat=3;
if 50<= age <=64 then agecat=4;
if 65<= age then agecat=5;
if (placdth GE 2) AND (placdth LE 7) then osh=1;
else if placdth=1 then osh=0;
if osh ne .;
if (monthdth in (&monthendcode)) OR (monthdth in (&monthbegincode));
keep season agecat rcu pi stateoc countyoc osh ;

run;
*************************************************************************************************;
*** Restricting to right county-season combinations *********************************************;
data a2;
set a1;

if (season = 1 and stateoc = "CA") or
 (season = 1 and stateoc = "CO") or
 (season = 1 and stateoc = "CT") or
 (season = 1 and stateoc = "GA") or
 (season = 1 and stateoc = "ID") or
 (season = 1 and stateoc = "MD") or
 (season = 1 and stateoc = "MI") or
 (season = 1 and stateoc = "MN") or
 (season = 1 and stateoc = "NM") or
 (season = 1 and stateoc = "NY") or
 (season = 1 and stateoc = "OH") or
 (season = 1 and stateoc = "OK") or
 (season = 1 and stateoc = "OR") or
 (season = 1 and stateoc = "RI") or
 (season = 1 and stateoc = "TN") or
 (season = 1 and stateoc = "UT") or
 (season = 2 and stateoc = "CA") or
 (season = 2 and stateoc = "CO") or
 (season = 2 and stateoc = "CT") or
 (season = 2 and stateoc = "GA") or
 (season = 2 and stateoc = "MD") or
 (season = 2 and stateoc = "MI") or
 (season = 2 and stateoc = "MN") or
 (season = 2 and stateoc = "NM") or
 (season = 2 and stateoc = "NY") or
 (season = 2 and stateoc = "OH") or
 (season = 2 and stateoc = "OR") or
 (season = 2 and stateoc = "RI") or
 (season = 2 and stateoc = "TN") or
 (season = 2 and stateoc = "UT") or
 (season = 3 and stateoc = "CA") or
 (season = 3 and stateoc = "CO") or
 (season = 3 and stateoc = "CT") or
 (season = 3 and stateoc = "GA") or
 (season = 3 and stateoc = "IA") or
 (season = 3 and stateoc = "MD") or
 (season = 3 and stateoc = "MI") or
 (season = 3 and stateoc = "MN") or
 (season = 3 and stateoc = "NM") or
 (season = 3 and stateoc = "NY") or
 (season = 3 and stateoc = "OH") or
 (season = 3 and stateoc = "OR") or
 (season = 3 and stateoc = "RI") or
 (season = 3 and stateoc = "TN") or
 (season = 3 and stateoc = "UT") or
 (season = 4 and stateoc = "CA") or
 (season = 4 and stateoc = "CO") or
 (season = 4 and stateoc = "CT") or
 (season = 4 and stateoc = "GA") or
 (season = 4 and stateoc = "MD") or
 (season = 4 and stateoc = "MI") or
 (season = 4 and stateoc = "MN") or
 (season = 4 and stateoc = "NM") or
 (season = 4 and stateoc = "NY") or
 (season = 4 and stateoc = "OH") or
 (season = 4 and stateoc = "OR") or
 (season = 4 and stateoc = "TN") or
 (season = 4 and stateoc = "UT") or
 (season = 5 and stateoc = "CA") or
 (season = 5 and stateoc = "CO") or
 (season = 5 and stateoc = "CT") or
 (season = 5 and stateoc = "GA") or
 (season = 5 and stateoc = "MD") or
 (season = 5 and stateoc = "MI") or
 (season = 5 and stateoc = "MN") or
 (season = 5 and stateoc = "NM") or
 (season = 5 and stateoc = "NY") or
 (season = 5 and stateoc = "OH") or
 (season = 5 and stateoc = "OR") or
 (season = 5 and stateoc = "TN") or
 (season = 5 and stateoc = "UT") or
 (season = 6 and stateoc = "CA") or
 (season = 6 and stateoc = "CO") or
 (season = 6 and stateoc = "CT") or
 (season = 6 and stateoc = "GA") or
 (season = 6 and stateoc = "MD") or
 (season = 6 and stateoc = "MI") or
 (season = 6 and stateoc = "MN") or
 (season = 6 and stateoc = "NM") or
 (season = 6 and stateoc = "NY") or
 (season = 6 and stateoc = "OH") or
 (season = 6 and stateoc = "OR") or
 (season = 6 and stateoc = "TN") or
 (season = 6 and stateoc = "UT") or
 (season = 7 and stateoc = "CA") or
 (season = 7 and stateoc = "CO") or
 (season = 7 and stateoc = "CT") or
 (season = 7 and stateoc = "GA") or
 (season = 7 and stateoc = "MD") or
 (season = 7 and stateoc = "MI") or
 (season = 7 and stateoc = "MN") or
 (season = 7 and stateoc = "NM") or
 (season = 7 and stateoc = "NY") or
 (season = 7 and stateoc = "OH") or
 (season = 7 and stateoc = "OR") or
 (season = 7 and stateoc = "TN") or
 (season = 7 and stateoc = "UT") or
 (season = 8 and stateoc = "CA") or
 (season = 8 and stateoc = "CO") or
 (season = 8 and stateoc = "CT") or
 (season = 8 and stateoc = "GA") or
 (season = 8 and stateoc = "MD") or
 (season = 8 and stateoc = "MI") or
 (season = 8 and stateoc = "MN") or
 (season = 8 and stateoc = "NM") or
 (season = 8 and stateoc = "NY") or
 (season = 8 and stateoc = "OH") or
 (season = 8 and stateoc = "OR") or
 (season = 8 and stateoc = "TN") or
 (season = 8 and stateoc = "UT") ;
 run;
 *************************************************************************************************;
*************************************************************************************************;
proc summary data= a2;
class season agecat osh;
var rcu;
output out=out1 sum=;
run;

data mort2010_17_rcu_osh;
set out1; 
if _Type_=7;
keep season agecat osh rcu;
run;
*************************************************************************************************;
*************************************************************************************************;
proc summary data= a2;
class season agecat osh;
var pi;
output out=out2 sum=;
run;

data mort2010_17_pi_osh;
set out2; 
if _Type_=7;
keep season agecat osh pi;
run;
*************************************************************************************************;
*************************************************************************************************;
*************************************************************************************************;
*************************************************************************************************;
proc sort data=mort2010_17_rcu_osh out=rcu2;
by season agecat osh rcu;
run;
*************************************************************************************************;
*************************************************************************************************;
proc sort data=mort2010_17_pi_osh out=pi2;
by season agecat osh pi;
run;
*************************************************************************************************;
*************************************************************************************************;
data mort2010_17_osh;
merge rcu2 pi2;
by season agecat osh;
*************************************************************************************************;
run;
/***********************************************************************************
*//***********************************************************************************
*//***********************************************************************************
*//***********************************************************************************
*/
PROC EXPORT DATA= WORK.mort2010_17_osh 
            OUTFILE= "\\cdc.gov\private\L322\VOR1\Documents\mort2010_17_season_osh_site_state.csv" 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;
