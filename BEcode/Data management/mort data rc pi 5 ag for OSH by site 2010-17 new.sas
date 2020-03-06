libname mydir "\\cdc\project\NCIRD_ID_EPI_Branch\ARM\Various data\NCHS Mortality 2010-2017";
/*CREATE DATE (OF DEATH) VARIABLE*/
data a1; 
set mydir.mort2010 
;
run;
***************************************************************************;
***************************************************************************;
data a1;
set a1;

%let monthendcode=10,11,12;
%let monthbegincode=1,2,3,4;

if (monthdth in (&monthendcode)) OR (monthdth in (&monthbegincode));
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


if age <=4 then agecat=1;
if 5<= age <=17 then agecat=2;
if 18<= age <=49 then agecat=3;
if 50<= age <=64 then agecat=4;
if 65<= age then agecat=5;

if (placdth GE 2) AND (placdth LE 7) then osh=1;
else if placdth=1 then osh=0;
if osh ne .;
if (monthdth in (&monthendcode)) OR (monthdth in (&monthbegincode));
keep agecat rcu pi stateoc countyoc osh ;

run;
*************************************************************************************************;
*** Restricting to right county-season combinations *********************************************;
data a2;
set a1;

if (stateoc = "CA" and countyoc = 1) or 
(stateoc = "CA" and countyoc = 13) or 
(stateoc = "CA" and countyoc = 75) or 
(stateoc = "CO" and countyoc = 1) or 
(stateoc = "CO" and countyoc = 5) or 
(stateoc = "CO" and countyoc = 31) or 
(stateoc = "CO" and countyoc = 35) or 
(stateoc = "CO" and countyoc = 59) or 
(stateoc = "CT" and countyoc = 3) or 
(stateoc = "CT" and countyoc = 7) or 
(stateoc = "CT" and countyoc = 9) or 
(stateoc = "GA" and countyoc = 63) or 
(stateoc = "GA" and countyoc = 67) or 
(stateoc = "GA" and countyoc = 89) or 
(stateoc = "GA" and countyoc = 97) or 
(stateoc = "GA" and countyoc = 121) or 
(stateoc = "GA" and countyoc = 135) or 
(stateoc = "GA" and countyoc = 217) or 
(stateoc = "GA" and countyoc = 247) or 
(stateoc = "ID" and countyoc = 1) or 
(stateoc = "ID" and countyoc = 5) or 
(stateoc = "ID" and countyoc = 55) or 
(stateoc = "MD" and countyoc = 3) or 
(stateoc = "MD" and countyoc = 5) or 
(stateoc = "MD" and countyoc = 510) or 
(stateoc = "MD" and countyoc = 13) or 
(stateoc = "MD" and countyoc = 25) or 
(stateoc = "MD" and countyoc = 27) or 
(stateoc = "MI" and countyoc = 37) or 
(stateoc = "MI" and countyoc = 45) or 
(stateoc = "MI" and countyoc = 65) or 
(stateoc = "MN" and countyoc = 3) or 
(stateoc = "MN" and countyoc = 19) or 
(stateoc = "MN" and countyoc = 37) or 
(stateoc = "MN" and countyoc = 53) or 
(stateoc = "MN" and countyoc = 123) or 
(stateoc = "MN" and countyoc = 139) or 
(stateoc = "MN" and countyoc = 163) or 
(stateoc = "NM" and countyoc = 1) or 
(stateoc = "NM" and countyoc = 5) or 
(stateoc = "NM" and countyoc = 13) or 
(stateoc = "NM" and countyoc = 17) or 
(stateoc = "NM" and countyoc = 29) or 
(stateoc = "NM" and countyoc = 45) or 
(stateoc = "NM" and countyoc = 49) or 
(stateoc = "NY" and countyoc = 1) or 
(stateoc = "NY" and countyoc = 21) or 
(stateoc = "NY" and countyoc = 39) or 
(stateoc = "NY" and countyoc = 57) or 
(stateoc = "NY" and countyoc = 83) or 
(stateoc = "NY" and countyoc = 91) or 
(stateoc = "NY" and countyoc = 93) or 
(stateoc = "NY" and countyoc = 95) or 
(stateoc = "NY" and countyoc = 37) or 
(stateoc = "NY" and countyoc = 51) or 
(stateoc = "NY" and countyoc = 55) or 
(stateoc = "NY" and countyoc = 69) or 
(stateoc = "NY" and countyoc = 73) or 
(stateoc = "NY" and countyoc = 117) or 
(stateoc = "NY" and countyoc = 123) or 
(stateoc = "OH" and countyoc = 41) or 
(stateoc = "OH" and countyoc = 45) or 
(stateoc = "OH" and countyoc = 49) or 
(stateoc = "OH" and countyoc = 89) or 
(stateoc = "OH" and countyoc = 97) or 
(stateoc = "OH" and countyoc = 117) or 
(stateoc = "OH" and countyoc = 129) or 
(stateoc = "OH" and countyoc = 159) or 
(stateoc = "OK" and countyoc = 21) or 
(stateoc = "OK" and countyoc = 31) or 
(stateoc = "OK" and countyoc = 47) or 
(stateoc = "OK" and countyoc = 109) or 
(stateoc = "OK" and countyoc = 123) or 
(stateoc = "OR" and countyoc = 5) or 
(stateoc = "OR" and countyoc = 51) or 
(stateoc = "OR" and countyoc = 67) or 
(stateoc = "RI" and countyoc = 7) or 
(stateoc = "TN" and countyoc = 21) or 
(stateoc = "TN" and countyoc = 37) or 
(stateoc = "TN" and countyoc = 43) or 
(stateoc = "TN" and countyoc = 147) or 
(stateoc = "TN" and countyoc = 149) or 
(stateoc = "TN" and countyoc = 165) or 
(stateoc = "TN" and countyoc = 187) or 
(stateoc = "TN" and countyoc = 189) or 
(stateoc = "UT" and countyoc = 35) then incl1=1;
else incl1=0;
 
if (stateoc = "CA" and countyoc = 1) or 
(stateoc = "CA" and countyoc = 13) or 
(stateoc = "CA" and countyoc = 75) or 
(stateoc = "CO" and countyoc = 1) or 
(stateoc = "CO" and countyoc = 5) or 
(stateoc = "CO" and countyoc = 31) or 
(stateoc = "CO" and countyoc = 35) or 
(stateoc = "CO" and countyoc = 59) or 
(stateoc = "CT" and countyoc = 3) or 
(stateoc = "CT" and countyoc = 7) or 
(stateoc = "CT" and countyoc = 9) or 
(stateoc = "GA" and countyoc = 63) or 
(stateoc = "GA" and countyoc = 67) or 
(stateoc = "GA" and countyoc = 89) or 
(stateoc = "GA" and countyoc = 97) or 
(stateoc = "GA" and countyoc = 121) or 
(stateoc = "GA" and countyoc = 135) or 
(stateoc = "GA" and countyoc = 217) or 
(stateoc = "GA" and countyoc = 247) or 
(stateoc = "MD" and countyoc = 3) or 
(stateoc = "MD" and countyoc = 5) or 
(stateoc = "MD" and countyoc = 510) or 
(stateoc = "MD" and countyoc = 13) or 
(stateoc = "MD" and countyoc = 25) or 
(stateoc = "MD" and countyoc = 27) or 
(stateoc = "MI" and countyoc = 37) or 
(stateoc = "MI" and countyoc = 45) or 
(stateoc = "MI" and countyoc = 65) or 
(stateoc = "MN" and countyoc = 3) or 
(stateoc = "MN" and countyoc = 19) or 
(stateoc = "MN" and countyoc = 37) or 
(stateoc = "MN" and countyoc = 53) or 
(stateoc = "MN" and countyoc = 123) or 
(stateoc = "MN" and countyoc = 139) or 
(stateoc = "MN" and countyoc = 163) or 
(stateoc = "NM" and countyoc = 1) or 
(stateoc = "NM" and countyoc = 5) or 
(stateoc = "NM" and countyoc = 13) or 
(stateoc = "NM" and countyoc = 17) or 
(stateoc = "NM" and countyoc = 29) or 
(stateoc = "NM" and countyoc = 45) or 
(stateoc = "NM" and countyoc = 49) or 
(stateoc = "NY" and countyoc = 1) or 
(stateoc = "NY" and countyoc = 21) or 
(stateoc = "NY" and countyoc = 39) or 
(stateoc = "NY" and countyoc = 57) or 
(stateoc = "NY" and countyoc = 83) or 
(stateoc = "NY" and countyoc = 91) or 
(stateoc = "NY" and countyoc = 93) or 
(stateoc = "NY" and countyoc = 95) or 
(stateoc = "NY" and countyoc = 37) or 
(stateoc = "NY" and countyoc = 51) or 
(stateoc = "NY" and countyoc = 55) or 
(stateoc = "NY" and countyoc = 69) or 
(stateoc = "NY" and countyoc = 73) or 
(stateoc = "NY" and countyoc = 117) or 
(stateoc = "NY" and countyoc = 123) or 
(stateoc = "OH" and countyoc = 41) or 
(stateoc = "OH" and countyoc = 45) or 
(stateoc = "OH" and countyoc = 49) or 
(stateoc = "OH" and countyoc = 89) or 
(stateoc = "OH" and countyoc = 97) or 
(stateoc = "OH" and countyoc = 117) or 
(stateoc = "OH" and countyoc = 129) or 
(stateoc = "OH" and countyoc = 159) or 
(stateoc = "OR" and countyoc = 5) or 
(stateoc = "OR" and countyoc = 51) or 
(stateoc = "OR" and countyoc = 67) or 
(stateoc = "RI" and countyoc = 7) or 
(stateoc = "TN" and countyoc = 21) or 
(stateoc = "TN" and countyoc = 37) or 
(stateoc = "TN" and countyoc = 43) or 
(stateoc = "TN" and countyoc = 147) or 
(stateoc = "TN" and countyoc = 149) or 
(stateoc = "TN" and countyoc = 165) or 
(stateoc = "TN" and countyoc = 187) or 
(stateoc = "TN" and countyoc = 189) or 
(stateoc = "UT" and countyoc = 35)  then incl2=1;
else incl2=0;

if (stateoc = "CA" and countyoc = 1) or 
(stateoc = "CA" and countyoc = 13) or 
(stateoc = "CA" and countyoc = 75) or 
(stateoc = "CO" and countyoc = 1) or 
(stateoc = "CO" and countyoc = 5) or 
(stateoc = "CO" and countyoc = 31) or 
(stateoc = "CO" and countyoc = 35) or 
(stateoc = "CO" and countyoc = 59) or 
(stateoc = "CT" and countyoc = 7) or 
(stateoc = "CT" and countyoc = 9) or 
(stateoc = "GA" and countyoc = 63) or 
(stateoc = "GA" and countyoc = 67) or 
(stateoc = "GA" and countyoc = 89) or 
(stateoc = "GA" and countyoc = 97) or 
(stateoc = "GA" and countyoc = 121) or 
(stateoc = "GA" and countyoc = 135) or 
(stateoc = "GA" and countyoc = 217) or 
(stateoc = "GA" and countyoc = 247) or 
(stateoc = "IA" and countyoc = 49) or 
(stateoc = "IA" and countyoc = 77) or 
(stateoc = "IA" and countyoc = 121) or 
(stateoc = "IA" and countyoc = 153) or 
(stateoc = "IA" and countyoc = 181) or 
(stateoc = "MD" and countyoc = 3) or 
(stateoc = "MD" and countyoc = 5) or 
(stateoc = "MD" and countyoc = 510) or 
(stateoc = "MD" and countyoc = 13) or 
(stateoc = "MD" and countyoc = 25) or 
(stateoc = "MD" and countyoc = 27) or 
(stateoc = "MI" and countyoc = 37) or 
(stateoc = "MI" and countyoc = 45) or 
(stateoc = "MI" and countyoc = 49) or 
(stateoc = "MI" and countyoc = 65) or 
(stateoc = "MN" and countyoc = 3) or 
(stateoc = "MN" and countyoc = 19) or 
(stateoc = "MN" and countyoc = 37) or 
(stateoc = "MN" and countyoc = 53) or 
(stateoc = "MN" and countyoc = 123) or 
(stateoc = "MN" and countyoc = 139) or 
(stateoc = "MN" and countyoc = 163) or 
(stateoc = "NM" and countyoc = 1) or 
(stateoc = "NM" and countyoc = 5) or 
(stateoc = "NM" and countyoc = 13) or 
(stateoc = "NM" and countyoc = 17) or 
(stateoc = "NM" and countyoc = 29) or 
(stateoc = "NM" and countyoc = 45) or 
(stateoc = "NM" and countyoc = 49) or 
(stateoc = "NY" and countyoc = 1) or 
(stateoc = "NY" and countyoc = 21) or 
(stateoc = "NY" and countyoc = 39) or 
(stateoc = "NY" and countyoc = 57) or 
(stateoc = "NY" and countyoc = 83) or 
(stateoc = "NY" and countyoc = 91) or 
(stateoc = "NY" and countyoc = 93) or 
(stateoc = "NY" and countyoc = 95) or 
(stateoc = "NY" and countyoc = 37) or 
(stateoc = "NY" and countyoc = 51) or 
(stateoc = "NY" and countyoc = 55) or 
(stateoc = "NY" and countyoc = 69) or 
(stateoc = "NY" and countyoc = 73) or 
(stateoc = "NY" and countyoc = 117) or 
(stateoc = "NY" and countyoc = 123) or 
(stateoc = "OH" and countyoc = 41) or 
(stateoc = "OH" and countyoc = 45) or 
(stateoc = "OH" and countyoc = 49) or 
(stateoc = "OH" and countyoc = 89) or 
(stateoc = "OH" and countyoc = 97) or 
(stateoc = "OH" and countyoc = 117) or 
(stateoc = "OH" and countyoc = 129) or 
(stateoc = "OH" and countyoc = 159) or 
(stateoc = "OR" and countyoc = 5) or 
(stateoc = "OR" and countyoc = 51) or 
(stateoc = "OR" and countyoc = 67) or 
(stateoc = "RI" and countyoc = 7) or 
(stateoc = "TN" and countyoc = 21) or 
(stateoc = "TN" and countyoc = 37) or 
(stateoc = "TN" and countyoc = 43) or 
(stateoc = "TN" and countyoc = 147) or 
(stateoc = "TN" and countyoc = 149) or 
(stateoc = "TN" and countyoc = 165) or 
(stateoc = "TN" and countyoc = 187) or 
(stateoc = "TN" and countyoc = 189) or 
(stateoc = "UT" and countyoc = 35)  then incl3=1;
else incl3=0;

if (stateoc = "CA" and countyoc = 1) or 
(stateoc = "CA" and countyoc = 13) or 
(stateoc = "CA" and countyoc = 75) or 
(stateoc = "CO" and countyoc = 1) or 
(stateoc = "CO" and countyoc = 5) or 
(stateoc = "CO" and countyoc = 31) or 
(stateoc = "CO" and countyoc = 35) or 
(stateoc = "CO" and countyoc = 59) or 
(stateoc = "CT" and countyoc = 7) or 
(stateoc = "CT" and countyoc = 9) or 
(stateoc = "GA" and countyoc = 63) or 
(stateoc = "GA" and countyoc = 67) or 
(stateoc = "GA" and countyoc = 89) or 
(stateoc = "GA" and countyoc = 97) or 
(stateoc = "GA" and countyoc = 121) or 
(stateoc = "GA" and countyoc = 135) or 
(stateoc = "GA" and countyoc = 217) or 
(stateoc = "GA" and countyoc = 247) or 
(stateoc = "MD" and countyoc = 3) or 
(stateoc = "MD" and countyoc = 5) or 
(stateoc = "MD" and countyoc = 510) or 
(stateoc = "MD" and countyoc = 13) or 
(stateoc = "MD" and countyoc = 25) or 
(stateoc = "MD" and countyoc = 27) or 
(stateoc = "MI" and countyoc = 37) or 
(stateoc = "MI" and countyoc = 45) or 
(stateoc = "MI" and countyoc = 49) or 
(stateoc = "MI" and countyoc = 65) or 
(stateoc = "MN" and countyoc = 3) or 
(stateoc = "MN" and countyoc = 19) or 
(stateoc = "MN" and countyoc = 37) or 
(stateoc = "MN" and countyoc = 53) or 
(stateoc = "MN" and countyoc = 123) or 
(stateoc = "MN" and countyoc = 139) or 
(stateoc = "MN" and countyoc = 163) or 
(stateoc = "NM" and countyoc = 1) or 
(stateoc = "NM" and countyoc = 5) or 
(stateoc = "NM" and countyoc = 13) or 
(stateoc = "NM" and countyoc = 17) or 
(stateoc = "NM" and countyoc = 29) or 
(stateoc = "NM" and countyoc = 45) or 
(stateoc = "NM" and countyoc = 49) or 
(stateoc = "NY" and countyoc = 1) or 
(stateoc = "NY" and countyoc = 21) or 
(stateoc = "NY" and countyoc = 39) or 
(stateoc = "NY" and countyoc = 57) or 
(stateoc = "NY" and countyoc = 83) or 
(stateoc = "NY" and countyoc = 91) or 
(stateoc = "NY" and countyoc = 93) or 
(stateoc = "NY" and countyoc = 95) or 
(stateoc = "NY" and countyoc = 37) or 
(stateoc = "NY" and countyoc = 51) or 
(stateoc = "NY" and countyoc = 55) or 
(stateoc = "NY" and countyoc = 69) or 
(stateoc = "NY" and countyoc = 73) or 
(stateoc = "NY" and countyoc = 117) or 
(stateoc = "NY" and countyoc = 123) or 
(stateoc = "OH" and countyoc = 41) or 
(stateoc = "OH" and countyoc = 45) or 
(stateoc = "OH" and countyoc = 49) or 
(stateoc = "OH" and countyoc = 73) or 
(stateoc = "OH" and countyoc = 89) or 
(stateoc = "OH" and countyoc = 97) or 
(stateoc = "OH" and countyoc = 117) or 
(stateoc = "OH" and countyoc = 127) or 
(stateoc = "OH" and countyoc = 129) or 
(stateoc = "OH" and countyoc = 159) or 
(stateoc = "OR" and countyoc = 5) or 
(stateoc = "OR" and countyoc = 51) or 
(stateoc = "OR" and countyoc = 67) or 
(stateoc = "TN" and countyoc = 21) or 
(stateoc = "TN" and countyoc = 37) or 
(stateoc = "TN" and countyoc = 43) or 
(stateoc = "TN" and countyoc = 147) or 
(stateoc = "TN" and countyoc = 149) or 
(stateoc = "TN" and countyoc = 165) or 
(stateoc = "TN" and countyoc = 187) or 
(stateoc = "TN" and countyoc = 189) or 
(stateoc = "UT" and countyoc = 35)  then incl4=1;
else incl4=0;

if (stateoc = "CA" and countyoc = 1) or 
(stateoc = "CA" and countyoc = 13) or 
(stateoc = "CA" and countyoc = 75) or 
(stateoc = "CO" and countyoc = 1) or 
(stateoc = "CO" and countyoc = 5) or 
(stateoc = "CO" and countyoc = 31) or 
(stateoc = "CO" and countyoc = 35) or 
(stateoc = "CO" and countyoc = 59) or 
(stateoc = "CT" and countyoc = 7) or 
(stateoc = "CT" and countyoc = 9) or 
(stateoc = "GA" and countyoc = 63) or 
(stateoc = "GA" and countyoc = 67) or 
(stateoc = "GA" and countyoc = 89) or 
(stateoc = "GA" and countyoc = 97) or 
(stateoc = "GA" and countyoc = 121) or 
(stateoc = "GA" and countyoc = 135) or 
(stateoc = "GA" and countyoc = 217) or 
(stateoc = "GA" and countyoc = 247) or 
(stateoc = "MD" and countyoc = 3) or 
(stateoc = "MD" and countyoc = 5) or 
(stateoc = "MD" and countyoc = 510) or 
(stateoc = "MD" and countyoc = 13) or 
(stateoc = "MD" and countyoc = 25) or 
(stateoc = "MD" and countyoc = 27) or 
(stateoc = "MI" and countyoc = 37) or 
(stateoc = "MI" and countyoc = 45) or 
(stateoc = "MI" and countyoc = 49) or 
(stateoc = "MI" and countyoc = 65) or 
(stateoc = "MN" and countyoc = 3) or 
(stateoc = "MN" and countyoc = 19) or 
(stateoc = "MN" and countyoc = 37) or 
(stateoc = "MN" and countyoc = 53) or 
(stateoc = "MN" and countyoc = 123) or 
(stateoc = "MN" and countyoc = 139) or 
(stateoc = "MN" and countyoc = 163) or 
(stateoc = "NM" and countyoc = 1) or 
(stateoc = "NM" and countyoc = 5) or 
(stateoc = "NM" and countyoc = 13) or 
(stateoc = "NM" and countyoc = 17) or 
(stateoc = "NM" and countyoc = 29) or 
(stateoc = "NM" and countyoc = 45) or 
(stateoc = "NM" and countyoc = 49) or 
(stateoc = "NY" and countyoc = 1) or 
(stateoc = "NY" and countyoc = 21) or 
(stateoc = "NY" and countyoc = 39) or 
(stateoc = "NY" and countyoc = 57) or 
(stateoc = "NY" and countyoc = 83) or 
(stateoc = "NY" and countyoc = 91) or 
(stateoc = "NY" and countyoc = 93) or 
(stateoc = "NY" and countyoc = 95) or 
(stateoc = "NY" and countyoc = 37) or 
(stateoc = "NY" and countyoc = 51) or 
(stateoc = "NY" and countyoc = 55) or 
(stateoc = "NY" and countyoc = 69) or 
(stateoc = "NY" and countyoc = 73) or 
(stateoc = "NY" and countyoc = 117) or 
(stateoc = "NY" and countyoc = 123) or 
(stateoc = "OH" and countyoc = 41) or 
(stateoc = "OH" and countyoc = 45) or 
(stateoc = "OH" and countyoc = 49) or 
(stateoc = "OH" and countyoc = 73) or 
(stateoc = "OH" and countyoc = 89) or 
(stateoc = "OH" and countyoc = 97) or 
(stateoc = "OH" and countyoc = 117) or 
(stateoc = "OH" and countyoc = 127) or 
(stateoc = "OH" and countyoc = 129) or 
(stateoc = "OH" and countyoc = 159) or 
(stateoc = "OR" and countyoc = 5) or 
(stateoc = "OR" and countyoc = 51) or 
(stateoc = "OR" and countyoc = 67) or 
(stateoc = "TN" and countyoc = 21) or 
(stateoc = "TN" and countyoc = 37) or 
(stateoc = "TN" and countyoc = 43) or 
(stateoc = "TN" and countyoc = 147) or 
(stateoc = "TN" and countyoc = 149) or 
(stateoc = "TN" and countyoc = 165) or 
(stateoc = "TN" and countyoc = 187) or 
(stateoc = "TN" and countyoc = 189) or 
(stateoc = "UT" and countyoc = 35)  then incl5=1;
else incl5=0;

if (stateoc = "CA" and countyoc = 1) or 
(stateoc = "CA" and countyoc = 13) or 
(stateoc = "CA" and countyoc = 75) or 
(stateoc = "CO" and countyoc = 1) or 
(stateoc = "CO" and countyoc = 5) or 
(stateoc = "CO" and countyoc = 31) or 
(stateoc = "CO" and countyoc = 35) or 
(stateoc = "CO" and countyoc = 59) or 
(stateoc = "CT" and countyoc = 7) or 
(stateoc = "CT" and countyoc = 9) or 
(stateoc = "GA" and countyoc = 63) or 
(stateoc = "GA" and countyoc = 67) or 
(stateoc = "GA" and countyoc = 89) or 
(stateoc = "GA" and countyoc = 97) or 
(stateoc = "GA" and countyoc = 121) or 
(stateoc = "GA" and countyoc = 135) or 
(stateoc = "GA" and countyoc = 217) or 
(stateoc = "GA" and countyoc = 247) or 
(stateoc = "MD" and countyoc = 3) or 
(stateoc = "MD" and countyoc = 5) or 
(stateoc = "MD" and countyoc = 510) or 
(stateoc = "MD" and countyoc = 13) or 
(stateoc = "MD" and countyoc = 25) or 
(stateoc = "MD" and countyoc = 27) or 
(stateoc = "MI" and countyoc = 37) or 
(stateoc = "MI" and countyoc = 45) or 
(stateoc = "MI" and countyoc = 49) or 
(stateoc = "MI" and countyoc = 65) or 
(stateoc = "MN" and countyoc = 3) or 
(stateoc = "MN" and countyoc = 19) or 
(stateoc = "MN" and countyoc = 37) or 
(stateoc = "MN" and countyoc = 53) or 
(stateoc = "MN" and countyoc = 123) or 
(stateoc = "MN" and countyoc = 139) or 
(stateoc = "MN" and countyoc = 163) or 
(stateoc = "NM" and countyoc = 1) or 
(stateoc = "NM" and countyoc = 5) or 
(stateoc = "NM" and countyoc = 13) or 
(stateoc = "NM" and countyoc = 17) or 
(stateoc = "NM" and countyoc = 29) or 
(stateoc = "NM" and countyoc = 45) or 
(stateoc = "NM" and countyoc = 49) or 
(stateoc = "NY" and countyoc = 1) or 
(stateoc = "NY" and countyoc = 21) or 
(stateoc = "NY" and countyoc = 39) or 
(stateoc = "NY" and countyoc = 57) or 
(stateoc = "NY" and countyoc = 83) or 
(stateoc = "NY" and countyoc = 91) or 
(stateoc = "NY" and countyoc = 93) or 
(stateoc = "NY" and countyoc = 95) or 
(stateoc = "NY" and countyoc = 37) or 
(stateoc = "NY" and countyoc = 51) or 
(stateoc = "NY" and countyoc = 55) or 
(stateoc = "NY" and countyoc = 69) or 
(stateoc = "NY" and countyoc = 73) or 
(stateoc = "NY" and countyoc = 117) or 
(stateoc = "NY" and countyoc = 123) or 
(stateoc = "OH" and countyoc = 41) or 
(stateoc = "OH" and countyoc = 45) or 
(stateoc = "OH" and countyoc = 49) or 
(stateoc = "OH" and countyoc = 73) or 
(stateoc = "OH" and countyoc = 89) or 
(stateoc = "OH" and countyoc = 97) or 
(stateoc = "OH" and countyoc = 117) or 
(stateoc = "OH" and countyoc = 127) or 
(stateoc = "OH" and countyoc = 129) or 
(stateoc = "OH" and countyoc = 159) or 
(stateoc = "OR" and countyoc = 5) or 
(stateoc = "OR" and countyoc = 51) or 
(stateoc = "OR" and countyoc = 67) or 
(stateoc = "TN" and countyoc = 21) or 
(stateoc = "TN" and countyoc = 37) or 
(stateoc = "TN" and countyoc = 43) or 
(stateoc = "TN" and countyoc = 147) or 
(stateoc = "TN" and countyoc = 149) or 
(stateoc = "TN" and countyoc = 165) or 
(stateoc = "TN" and countyoc = 187) or 
(stateoc = "TN" and countyoc = 189) or 
(stateoc = "UT" and countyoc = 35)  then incl6=1;
else incl6=0;

if (stateoc = "CA" and countyoc = 1) or 
(stateoc = "CA" and countyoc = 13) or 
(stateoc = "CA" and countyoc = 75) or 
(stateoc = "CO" and countyoc = 1) or 
(stateoc = "CO" and countyoc = 5) or 
(stateoc = "CO" and countyoc = 31) or 
(stateoc = "CO" and countyoc = 35) or 
(stateoc = "CO" and countyoc = 59) or 
(stateoc = "CT" and countyoc = 7) or 
(stateoc = "CT" and countyoc = 9) or 
(stateoc = "GA" and countyoc = 63) or 
(stateoc = "GA" and countyoc = 67) or 
(stateoc = "GA" and countyoc = 89) or 
(stateoc = "GA" and countyoc = 97) or 
(stateoc = "GA" and countyoc = 121) or 
(stateoc = "GA" and countyoc = 135) or 
(stateoc = "GA" and countyoc = 217) or 
(stateoc = "GA" and countyoc = 247) or 
(stateoc = "MD" and countyoc = 3) or 
(stateoc = "MD" and countyoc = 5) or 
(stateoc = "MD" and countyoc = 510) or 
(stateoc = "MD" and countyoc = 13) or 
(stateoc = "MD" and countyoc = 25) or 
(stateoc = "MD" and countyoc = 27) or 
(stateoc = "MI" and countyoc = 37) or 
(stateoc = "MI" and countyoc = 45) or 
(stateoc = "MI" and countyoc = 49) or 
(stateoc = "MI" and countyoc = 65) or 
(stateoc = "MN" and countyoc = 3) or 
(stateoc = "MN" and countyoc = 19) or 
(stateoc = "MN" and countyoc = 37) or 
(stateoc = "MN" and countyoc = 53) or 
(stateoc = "MN" and countyoc = 123) or 
(stateoc = "MN" and countyoc = 139) or 
(stateoc = "MN" and countyoc = 163) or 
(stateoc = "NM" and countyoc = 1) or 
(stateoc = "NM" and countyoc = 5) or 
(stateoc = "NM" and countyoc = 13) or 
(stateoc = "NM" and countyoc = 17) or 
(stateoc = "NM" and countyoc = 29) or 
(stateoc = "NM" and countyoc = 45) or 
(stateoc = "NM" and countyoc = 49) or 
(stateoc = "NY" and countyoc = 1) or 
(stateoc = "NY" and countyoc = 21) or 
(stateoc = "NY" and countyoc = 39) or 
(stateoc = "NY" and countyoc = 57) or 
(stateoc = "NY" and countyoc = 83) or 
(stateoc = "NY" and countyoc = 91) or 
(stateoc = "NY" and countyoc = 93) or 
(stateoc = "NY" and countyoc = 95) or 
(stateoc = "NY" and countyoc = 37) or 
(stateoc = "NY" and countyoc = 51) or 
(stateoc = "NY" and countyoc = 55) or 
(stateoc = "NY" and countyoc = 69) or 
(stateoc = "NY" and countyoc = 73) or 
(stateoc = "NY" and countyoc = 117) or 
(stateoc = "NY" and countyoc = 123) or 
(stateoc = "OH" and countyoc = 41) or 
(stateoc = "OH" and countyoc = 45) or 
(stateoc = "OH" and countyoc = 49) or 
(stateoc = "OH" and countyoc = 73) or 
(stateoc = "OH" and countyoc = 89) or 
(stateoc = "OH" and countyoc = 97) or 
(stateoc = "OH" and countyoc = 117) or 
(stateoc = "OH" and countyoc = 127) or 
(stateoc = "OH" and countyoc = 129) or 
(stateoc = "OH" and countyoc = 159) or 
(stateoc = "OR" and countyoc = 5) or 
(stateoc = "OR" and countyoc = 51) or 
(stateoc = "OR" and countyoc = 67) or 
(stateoc = "TN" and countyoc = 21) or 
(stateoc = "TN" and countyoc = 37) or 
(stateoc = "TN" and countyoc = 43) or 
(stateoc = "TN" and countyoc = 147) or 
(stateoc = "TN" and countyoc = 149) or 
(stateoc = "TN" and countyoc = 165) or 
(stateoc = "TN" and countyoc = 187) or 
(stateoc = "TN" and countyoc = 189) or 
(stateoc = "UT" and countyoc = 35)  then incl7=1;
else incl7=0;

if (stateoc = "CA" and countyoc = 1) or 
(stateoc = "CA" and countyoc = 13) or 
(stateoc = "CA" and countyoc = 75) or 
(stateoc = "CO" and countyoc = 1) or 
(stateoc = "CO" and countyoc = 5) or 
(stateoc = "CO" and countyoc = 31) or 
(stateoc = "CO" and countyoc = 35) or 
(stateoc = "CO" and countyoc = 59) or 
(stateoc = "CT" and countyoc = 7) or 
(stateoc = "CT" and countyoc = 9) or 
(stateoc = "GA" and countyoc = 63) or 
(stateoc = "GA" and countyoc = 67) or 
(stateoc = "GA" and countyoc = 89) or 
(stateoc = "GA" and countyoc = 97) or 
(stateoc = "GA" and countyoc = 121) or 
(stateoc = "GA" and countyoc = 135) or 
(stateoc = "GA" and countyoc = 217) or 
(stateoc = "GA" and countyoc = 247) or 
(stateoc = "MD" and countyoc = 3) or 
(stateoc = "MD" and countyoc = 5) or 
(stateoc = "MD" and countyoc = 510) or 
(stateoc = "MD" and countyoc = 13) or 
(stateoc = "MD" and countyoc = 25) or 
(stateoc = "MD" and countyoc = 27) or 
(stateoc = "MI" and countyoc = 37) or 
(stateoc = "MI" and countyoc = 45) or 
(stateoc = "MI" and countyoc = 49) or 
(stateoc = "MI" and countyoc = 65) or 
(stateoc = "MI" and countyoc = 161) or 
(stateoc = "MN" and countyoc = 3) or 
(stateoc = "MN" and countyoc = 19) or 
(stateoc = "MN" and countyoc = 37) or 
(stateoc = "MN" and countyoc = 53) or 
(stateoc = "MN" and countyoc = 123) or 
(stateoc = "MN" and countyoc = 139) or 
(stateoc = "MN" and countyoc = 163) or 
(stateoc = "NM" and countyoc = 1) or 
(stateoc = "NM" and countyoc = 5) or 
(stateoc = "NM" and countyoc = 13) or 
(stateoc = "NM" and countyoc = 17) or 
(stateoc = "NM" and countyoc = 29) or 
(stateoc = "NM" and countyoc = 45) or 
(stateoc = "NM" and countyoc = 49) or 
(stateoc = "NY" and countyoc = 1) or 
(stateoc = "NY" and countyoc = 21) or 
(stateoc = "NY" and countyoc = 39) or 
(stateoc = "NY" and countyoc = 57) or 
(stateoc = "NY" and countyoc = 83) or 
(stateoc = "NY" and countyoc = 91) or 
(stateoc = "NY" and countyoc = 93) or 
(stateoc = "NY" and countyoc = 95) or 
(stateoc = "NY" and countyoc = 37) or 
(stateoc = "NY" and countyoc = 51) or 
(stateoc = "NY" and countyoc = 55) or 
(stateoc = "NY" and countyoc = 69) or 
(stateoc = "NY" and countyoc = 73) or 
(stateoc = "NY" and countyoc = 117) or 
(stateoc = "NY" and countyoc = 123) or 
(stateoc = "OH" and countyoc = 41) or 
(stateoc = "OH" and countyoc = 45) or 
(stateoc = "OH" and countyoc = 49) or 
(stateoc = "OH" and countyoc = 73) or 
(stateoc = "OH" and countyoc = 89) or 
(stateoc = "OH" and countyoc = 97) or 
(stateoc = "OH" and countyoc = 117) or 
(stateoc = "OH" and countyoc = 127) or 
(stateoc = "OH" and countyoc = 129) or 
(stateoc = "OH" and countyoc = 159) or 
(stateoc = "OR" and countyoc = 5) or 
(stateoc = "OR" and countyoc = 51) or 
(stateoc = "OR" and countyoc = 67) or 
(stateoc = "TN" and countyoc = 21) or 
(stateoc = "TN" and countyoc = 37) or 
(stateoc = "TN" and countyoc = 43) or 
(stateoc = "TN" and countyoc = 147) or 
(stateoc = "TN" and countyoc = 149) or 
(stateoc = "TN" and countyoc = 165) or 
(stateoc = "TN" and countyoc = 187) or 
(stateoc = "TN" and countyoc = 189) or 
(stateoc = "UT" and countyoc = 35 )  then incl8=1;
else incl8=0;

if incl1=1 or incl2=1 or incl3=1 or incl4=1 or incl5=1 or incl6=1 or incl7=1 or incl8=1;

keep agecat rcu pi osh incl1 incl2 incl3 incl4 incl5 incl6 incl7 incl8; 
run ;
*************************************************************************************************;
*************************************************************************************************;
data a3; 
set a2;
if incl1=1;
season=1;
run;
proc summary data= a3;
class season agecat osh;
var rcu;
output out=out1a sum=;
run;

data mort2010_11_rcu_osh;
set out1a; 
if _Type_=7;
keep season agecat osh rcu;
run;
*************************************************************************************************;
proc summary data= a3;
class season agecat osh;
var pi;
output out=out1b sum=;
run;

data mort2010_11_pi_osh;
set out1b; 
if _Type_=7;
keep season agecat osh pi;
run;
*************************************************************************************************;
*************************************************************************************************;
data a3; 
set a2;
if incl2=1;
season=2;
run;
proc summary data= a3;
class season agecat osh;
var rcu;
output out=out2a sum=;
run;

data mort2011_12_rcu_osh;
set out2a; 
if _Type_=7;
keep season agecat osh rcu;
run;
*************************************************************************************************;
proc summary data= a3;
class season agecat osh;
var pi;
output out=out2b sum=;
run;

data mort2011_12_pi_osh;
set out2b; 
if _Type_=7;
keep season agecat osh pi;
run;
*************************************************************************************************;
*************************************************************************************************;
data a3; 
set a2;
if incl3=1;
season=3;
run;
proc summary data= a3;
class season agecat osh;
var rcu;
output out=out3a sum=;
run;

data mort2012_13_rcu_osh;
set out3a; 
if _Type_=7;
keep season agecat osh rcu;
run;
*************************************************************************************************;
proc summary data= a3;
class season agecat osh;
var pi;
output out=out3b sum=;
run;

data mort2012_13_pi_osh;
set out3b; 
if _Type_=7;
keep season agecat osh pi;
run;
*************************************************************************************************;
*************************************************************************************************;
data a3; 
set a2;
if incl4=1;
season=4;
run;
proc summary data= a3;
class season agecat osh;
var rcu;
output out=out4a sum=;
run;

data mort2013_14_rcu_osh;
set out4a; 
if _Type_=7;
keep season agecat osh rcu;
run;
*************************************************************************************************;
proc summary data= a3;
class season agecat osh;
var pi;
output out=out4b sum=;
run;

data mort2013_14_pi_osh;
set out4b; 
if _Type_=7;
keep season agecat osh pi;
run;
*************************************************************************************************;
*************************************************************************************************;
data a3; 
set a2;
if incl5=1;
season=5;
run;
proc summary data= a3;
class season agecat osh;
var rcu;
output out=out5a sum=;
run;

data mort2014_15_rcu_osh;
set out5a; 
if _Type_=7;
keep season agecat osh rcu;
run;
*************************************************************************************************;
proc summary data= a3;
class season agecat osh;
var pi;
output out=out5b sum=;
run;

data mort2014_15_pi_osh;
set out5b; 
if _Type_=7;
keep season agecat osh pi;
run;
*************************************************************************************************;
*************************************************************************************************;
data a3; 
set a2;
if incl6=1;
season=6;
run;
proc summary data= a3;
class season agecat osh;
var rcu;
output out=out6a sum=;
run;

data mort2015_16_rcu_osh;
set out6a; 
if _Type_=7;
keep season agecat osh rcu;
run;
*************************************************************************************************;
proc summary data= a3;
class season agecat osh;
var pi;
output out=out6b sum=;
run;

data mort2015_16_pi_osh;
set out6b; 
if _Type_=7;
keep season agecat osh pi;
run;
*************************************************************************************************;
*************************************************************************************************;
data a3; 
set a2;
if incl7=1;
season=7;
run;
proc summary data= a3;
class season agecat osh;
var rcu;
output out=out7a sum=;
run;

data mort2016_17_rcu_osh;
set out7a; 
if _Type_=7;
keep season agecat osh rcu;
run;
*************************************************************************************************;
proc summary data= a3;
class season agecat osh;
var pi;
output out=out7b sum=;
run;

data mort2016_17_pi_osh;
set out7b; 
if _Type_=7;
keep season agecat osh pi;
run;
*************************************************************************************************;
*************************************************************************************************;
*************************************************************************************************;
proc sort data=mort2010_11_rcu_osh out=rcu1;
by season agecat osh;
run;
proc sort data=mort2010_11_pi_osh out=pi1;
by season agecat osh;
run;
*************************************************************************************************;
*************************************************************************************************;
proc sort data=mort2011_12_rcu_osh out=rcu2;
by season agecat osh;
run;
proc sort data=mort2011_12_pi_osh out=pi2;
by season agecat osh;
run;
*************************************************************************************************;
*************************************************************************************************;
proc sort data=mort2012_13_rcu_osh out=rcu3;
by season agecat osh;
run;
proc sort data=mort2012_13_pi_osh out=pi3;
by season agecat osh;
run;
*************************************************************************************************;
*************************************************************************************************;
proc sort data=mort2013_14_rcu_osh out=rcu4;
by season agecat osh;
run;
proc sort data=mort2013_14_pi_osh out=pi4;
by season agecat osh;
run;
*************************************************************************************************;
*************************************************************************************************;
proc sort data=mort2014_15_rcu_osh out=rcu5;
by season agecat osh;
run;
proc sort data=mort2014_15_pi_osh out=pi5;
by season agecat osh;
run;
*************************************************************************************************;
*************************************************************************************************;
proc sort data=mort2015_16_rcu_osh out=rcu6;
by season agecat osh;
run;
proc sort data=mort2015_16_pi_osh out=pi6;
by season agecat osh;
run;
*************************************************************************************************;
*************************************************************************************************;
proc sort data=mort2016_17_rcu_osh out=rcu7;
by season agecat osh;
run;
proc sort data=mort2016_17_pi_osh out=pi7;
by season agecat osh;
run;
*************************************************************************************************;
*************************************************************************************************;
data mort2010_17_osh;
merge rcu1 pi1 rcu2 pi2 rcu3 pi3 rcu4 pi4 rcu5 pi5 rcu6 pi6 rcu7 pi7;
by season agecat osh ;
run;
/***********************************************************************************
*//***********************************************************************************
*//***********************************************************************************
*//***********************************************************************************
*/
PROC EXPORT DATA= WORK.mort2010_17_osh 
            OUTFILE= "\\cdc.gov\private\L322\VOR1\Documents\Influenza Research\Projects\Excess mortality\EMAnalysis\EMdata\mort2010_17_season_osh.csv" 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;
