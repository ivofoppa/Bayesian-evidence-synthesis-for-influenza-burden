#### Processing data by site
### Ivo M. Foppa, 12/2018
### The directory containing all datasets has to be defined here:
bfolder <- "//cdc.gov/project/NCIRD_ID_EPI_Branch/ARM/Disease Burden/_Detection/Data/" ## project folder
setwd(paste0(bfolder,''))
#########################################################################################
#########################################################################################
#########################################################################################
library(readxl)
#########################################################################################
#########################################################################################
#########################################################################################
PROC IMPORT OUT= WORK.ALL0 	
read_excel('FULL_EIPBURDEN_DATE_v2.xls')
DBMS=XLS REPLACE/
SHEET="PartA"/ 
GETNAMES=YES/
RUN/

# ### 2011-12 season (from CA1, CA2, NM, NY, OR) ### #/ 
PROC IMPORT OUT= WORK.ALL1 	
read_excel('FULL_EIPBURDEN_201112.xls')
DBMS=XLS REPLACE/
SHEET="PartA"/ 
GETNAMES=YES/
RUN/

# ############################################################################################# #/

# ############################################################################################# #/
# ### California ### #/
# ### 2012-13 & 2013-14 season ### #/
PROC IMPORT OUT= WORK.CA2 	
read_excel('CA_InfluenzaBurden_012916.xlsx')
DBMS=XLSX REPLACE/
SHEET="Sheet1"/ 
GETNAMES=YES/
RUN/
# ### 2012-13 season ### #/
PROC IMPORT OUT= WORK.CAK1 	
read_excel('CA_Kaiser_Flu Burden_2012_13.xlsx')
DBMS=XLSX REPLACE/
SHEET="Sheet1"/ 
GETNAMES=YES/
RUN/
# ### 2013-14 season ### #/
PROC IMPORT OUT= WORK.CAK2 	
read_excel('CA_Kaiser_Flu Burden_2013_14.xlsx')
DBMS=XLSX REPLACE/
SHEET="Sheet1"/ 
GETNAMES=YES/
RUN/
# ### 2014-15 season ### #/
PROC IMPORT OUT= WORK.CAK3 	
read_excel('CA_Kaiser flu burden_2014_15.xlsx')
DBMS=XLSX REPLACE/
SHEET="Sheet1"/ 
GETNAMES=YES/
RUN/

# ### 2015-16 season ### #/
PROC IMPORT OUT= WORK.CAK4 	
read_excel('Kaiser_1516_burden_to CDC.xlsx')
DBMS=XLSX REPLACE/
SHEET="Sheet1"/ 
GETNAMES=YES/
RUN/

PROC IMPORT OUT= WORK.CA3 	
read_excel('CA_1516_BurdenDatabase_Final.xlsx')
DBMS=XLSX REPLACE/
SHEET="Patient data"/ 
GETNAMES=YES/
RUN/
# ############################################################################################# #/

# ############################################################################################# #/
# ### Colorado ### #/
# ### 2012-13 & 2013-14 season ### #/
PROC IMPORT OUT= WORK.CO2 	
read_excel('CO_BURDEN_12032014.xlsx')
DBMS=XLSX REPLACE/
SHEET="Sheet1"/ 
GETNAMES=YES/
RUN/
# ### 2014-15 season ### #/
PROC IMPORT OUT= WORK.CO3 	
read_excel('CO_BURDEN_201415.xlsx')
DBMS=XLSX REPLACE/
SHEET="Sheet1"/ 
GETNAMES=YES/
RUN/

# ### 2015-16 season ### #/
PROC IMPORT OUT= WORK.CO4 	
read_excel('CO_Burden_1516_10032017.xlsx')
DBMS=XLSX REPLACE/
SHEET="Patient data"/ 
GETNAMES=YES/
RUN/
# ############################################################################################# #/

# ############################################################################################# #/
# ### Connecticut ### #/
# ### 2012-13 season ### #/
PROC IMPORT OUT= WORK.CT2 	
read_excel('CT_BURDEN_04272015.xlsx')
DBMS=XLSX REPLACE/
SHEET="Sheet1"/ 
GETNAMES=YES/
RUN/
# ### 2013-14 season ### #/
PROC IMPORT OUT= WORK.CT3 	
read_excel('CT_EIPBURDEN_12072015.xlsx')
DBMS=XLSX REPLACE/
SHEET="Sheet1"/ 
GETNAMES=YES/
RUN/
# ### 2014-15 season ### #/
PROC IMPORT OUT= WORK.CT4 	
read_excel('CT_BURDEN_1415_11222016.xlsx')
DBMS=XLSX REPLACE/
SHEET="Sheet1"/ 
GETNAMES=YES/
RUN/

# ### 2015-16 season ### #/
PROC IMPORT OUT= WORK.CT5 	
read_excel('CT_BURDEN_1516_10312017.xlsx')
DBMS=XLSX REPLACE/
SHEET="Patient data"/ 
GETNAMES=YES/
RUN/
# ############################################################################################# #/

# ############################################################################################# #/
# ### Georgia ### #/
# ### 2010-11 & 2011-12 & 2012-13 seasons ### #/
PROC IMPORT OUT= WORK.GA 	
read_excel('GA_EIPBURDEN_04052016.xlsx')
DBMS=XLSX REPLACE/
SHEET="Sheet1"/ 
GETNAMES=YES/
RUN/

# ### 2013-14 & 2014-15 seasons ### #/
PROC IMPORT OUT= WORK.GA2 	
read_excel('GA Influenza Burden 2013-14_2014-15 CDC Submission 9-14-2017.xlsx')
DBMS=XLSX REPLACE/
SHEET="Sheet1"/ 
GETNAMES=YES/
RUN/

# ### 2015-16 seasons ### #/
PROC IMPORT OUT= WORK.GA3 	
read_excel('FluBurdenSample_2015_2016_CDC Submission.xls')
DBMS=XLS REPLACE/
SHEET="Final CDC Format"/ 
GETNAMES=YES/
RUN/
# ############################################################################################# #/


# ############################################################################################# #/
# ### Michigan ### #/
# ### 2015-16 season ### #/
PROC IMPORT OUT= WORK.MI 	
read_excel('MI_BURDEN_1516_06132018.xlsx')
DBMS=XLSX REPLACE/
SHEET="Patient data"/ 
GETNAMES=YES/
RUN/

# ############################################################################################# #/
# ### Minnesota ### #/
# ### 2012-13 & 2013-14 season ### #/
PROC IMPORT OUT= WORK.MN 	
read_excel('MN_BURDEN_03092016.xlsx')
DBMS=XLSX REPLACE/
SHEET="Sheet1"/ 
GETNAMES=YES/
RUN/

# ### 2014-15 season ### #/
PROC IMPORT OUT= WORK.MN2 	
read_excel('STATE_EIPBURDEN_DATE_v2014_MN_1415Data.xlsx')
DBMS=XLSX REPLACE/
SHEET="Sheet1"/ 
GETNAMES=YES/
RUN/

# ### 2015-16 season ### #/
PROC IMPORT OUT= WORK.MN3 	
read_excel('MN_BURDEN_1516_12192017.xlsx')
DBMS=XLSX REPLACE/
SHEET="Patient data"/ 
GETNAMES=YES/
RUN/
# ############################################################################################# #/

# ############################################################################################# #/
# ### New Mexico ### #/
# ### 2012-13 & 2013-14 season ### #/
PROC IMPORT OUT= WORK.NM2 	
read_excel('NM_EIPBURDEN_2015.05.05.xlsx')
DBMS=XLSX REPLACE/
SHEET="Sheet1"/ 
GETNAMES=YES/
RUN/
# ### 2014-15 season ### #/
PROC IMPORT OUT= WORK.NM3	
read_excel('2014-15 EIPBURDEN NM Final.xlsx')
DBMS=XLSX REPLACE/
SHEET="Sheet1"/ 
GETNAMES=YES/
RUN/
# ### 2015-16 season ### #/
PROC IMPORT OUT= WORK.NM4	
read_excel('NM_Burden_1516_12212017.xls')
DBMS=XLS REPLACE/
SHEET="Patient data"/ 
GETNAMES=YES/
RUN/
# ############################################################################################# #/

# ############################################################################################# #/
# ### New York ### #/
# ### 2012-13 & 2013-14 season ### #/
PROC IMPORT OUT= WORK.NY2 	
read_excel('STATE_EIPBURDEN_DATE_v2014_9.23.14.xlsx')
DBMS=XLSX REPLACE/
SHEET="Sheet1"/ 
GETNAMES=YES/
RUN/
# ### 2014-15 season ### #/
PROC IMPORT OUT= WORK.NY3 	
read_excel('STATE_EIPBURDEN_NYR_2014-2015.xlsx')
DBMS=XLSX REPLACE/
SHEET="Sheet1"/ 
GETNAMES=YES/
RUN/
# ### 2015-16 season ### #/
## NYR #/
  PROC IMPORT OUT= WORK.NY4 	
read_excel('NY_BURDEN_1516_03012018.xlsx')
DBMS=XLSX REPLACE/
SHEET="Patient data"/ 
GETNAMES=YES/
RUN/

## NYA #/
  PROC IMPORT OUT= WORK.NY5 	
read_excel('NY_BURDEN_1516_06072018.xlsx')
DBMS=XLSX REPLACE/
SHEET="Patient data"/ 
GETNAMES=YES/
RUN/
# ############################################################################################# #/

# ############################################################################################# #/
# ### Oregon ### #/
# ### 2012-13 season ### #/
PROC IMPORT OUT= WORK.OR2 	
read_excel('OR_Burden_12122014.xlsx')
DBMS=XLSX REPLACE/
SHEET="Sheet1"/ 
GETNAMES=YES/
RUN/
# ### 2013-14 season ### #/
PROC IMPORT OUT= WORK.OR3 	
read_excel('OR_EIPBURDEN_08242015_v2014.xlsx')
DBMS=XLSX REPLACE/
SHEET="Sheet1"/ 
GETNAMES=YES/
RUN/
# ### 2014-15 season ### #/
PROC IMPORT OUT= WORK.OR4 	
read_excel('STATE_EIPBURDEN_DATE_OR_v2014.xlsx')
DBMS=XLSX REPLACE/
SHEET="Sheet1"/ 
GETNAMES=YES/
RUN/

# ### 2015-16 season ### #/
PROC IMPORT OUT= WORK.OR5 	
read_excel('OR_BURDEN_1516_04102018.xlsx')
DBMS=XLSX REPLACE/
SHEET="Patient data"/ 
GETNAMES=YES/
RUN/
# ############################################################################################# #/

# ############################################################################################# #/
# ### Tennessee ### #/
# ### 2012-13 & 2013-14 & 2014-15 season ### #/
PROC IMPORT OUT= WORK.TN 	
read_excel('TN_EIPBURDEN_20170112.xlsx')
DBMS=XLSX REPLACE/
SHEET="Sheet1"/ 
GETNAMES=YES/
RUN/
# ### 2015-16 season ### #/
PROC IMPORT OUT= WORK.TN2 	
read_excel('TN_EIPBURDEN_20180710.xlsx')
DBMS=XLSX REPLACE/ 
GETNAMES=YES/
RUN/
# ############################################################################################# #/


# #### Data cleaning: in preparation for combining datasets #### #/
Data CA2b (DROP = __:)/
set CA2(rename=
          (TestType =__TestType
           TestResult  =__TestResult
           TestType2 = __TestType2
           TestResult2 =__TestResult2)
)/
TestType = input(__TestType,1.)/
TestResult = input(__TestResult,1.)/
TestType2 = input(__TestType2,1.)/
TestResult2 = input(__TestResult2,1.)/
run/

# ### in 2014/15, California-Kaiser did not include test type in their dataset. 
They only use PCR, so if the patient was tested it was using PCR (testtype = 1) ### #/
data cak3b/ set cak3/ 
if testedflu = 1 then testtype = 1/
run/

data MNb /
set MN/ #(rename=
            (DateHosp =__DateHosp
             TestType =__TestType
             TestResult  =__TestResult
             TestType2 = __TestType2
             TestResult2 =__TestResult2
             DateDeath =__DateDeath
             ICD9_1 =__ICD9_1)
)/
# DateHosp = input(__DateHosp,best.)/
# TestType = input(__TestType,1.)/
# TestResult = input(__TestResult,1.)/
# TestType2 = input(__TestType2,1.)/
# TestResult2 = input(__TestResult2,1.)/
# DateDeath = input(__DateDeath,best.)/
# ICD9_1 = input(__ICD9_1,$char6.)/
# format DateDeath date9./
run/


data GAb ## (DROP = __:)/
set GA/ #(rename=
            (DateHosp =__DateHosp
             TestType =__TestType
             TestResult =__TestResult
             ICD9_1 =__ICD9_1)
)/

#  DateHosp = input(__DateHosp,best.)/
#  TestType = input(__TestType,1.)/
#  TestResult = input(__TestResult,1.)/
#  ICD9_1 = input(__ICD9_1,$char6.)/
run/

data GA2b/
set GA2 (rename=
           (ICD9_1 = __ICD9_1
            ICU = __ICU
            Died = __Died)
)/
ICD9_1 = input(__ICD9_1, $char6.)/

state = "GA"/

if Flu_Test_Type_1 = "Rapid" 	then TestType = 2/
if Flu_Test_Type_2 = "Rapid" 	then TestType2 = 2/ 
if Flu_Test_Type_1 = "RT-PCR" 	then TestType = 1/
if Flu_Test_Type_2 = "RT-PCR" 	then TestType2 = 1/
if Flu_Test_Type_1 = "PCR" 		then TestType = 1/
if Flu_Test_Type_2 = "PCR" 		then TestType2 = 1/

if FluTest_1__Yes_No_ = "Yes" 	then testedflu = 1/
if FluTest_1__Yes_No_ = "No"  	then testedflu = 2/

if Flu_Test_Result_1 in ("Positive A", "Positive B", "Positive Unknown", "Positive A / H1N1")
then TestResult = 1/
if Flu_Test_Result_2 in ("Positive A", "Positive B", "Positive Unknown", "Positive A / H1N1")
then TestResult2 = 1/
if Flu_Test_Result_1 = "Negative" then TestResult = 2/
if Flu_Test_Result_2 = "Negative" then TestResult = 2/

if __ICU = "Yes" 				then ICU = 1/
if __ICU = "No" 				then ICU = 2/
if __ICU = "Unk"				then ICU = 9/

if __Died = "Yes" 				then Died = 1/
if __Died = "No" 				then Died = 2/
if __Died = "Unk"				then Died = 9/

age = Age___Years_/
datehosp = admit_date/
month=put(datehosp, monname3.)/
drop __ICD9_1 __ICU __Died Age___Years_/
run/

data NM3b (DROP = __:)/
set NM3 (rename=
           (TestType2 = __TestType2
            TestResult2 =__TestResult2)
)/

TestType2 = input(__TestType2,1.)/
TestResult2 = input(__TestResult2,1.)/
run/

data TNb (DROP = __:)/
set TN (rename=
          (PatientNo_ = __PatientNo_
           TestType2 = __TestType2
           TestResult2 =__TestResult2
           ICD9_1 =__ICD9_1)
)/ 
PatientNo_ = input(__PatientNo_, $char6.)/
TestType2 = input(__TestType2,1.)/
TestResult2 = input(__TestResult2,1.)/
ICD9_1 = input(__ICD9_1,$char6.)/
run/


## A OHalloran - 07/10/2018 - rename variables so they are consistent with previous seasons 
create age, state, and year variables
For now, only keeping essential variables because variables across sites have different
types which won't allow concatenation #/

data ca3b (keep=testedflu testtype testtype2 testresult testresult2 died hospid datehosp
dob year age state)/

set ca3 (rename= (FluTst=TestedFlu
FluTstTyp1=TestType
FluTstTyp2=TestType2
FluResult1=TestResult
FluResult2=TestResult2
Outcome=Died
hosp_tx=hospid
dob=_dob))/
dob=input(dob,mmddyy10.)/
datehosp=input(admdate,mmddyy10.)/
year=year(datehosp)/
age = floor((intck("month", dob, datehosp,'C')) / 12)/
state = "CA"/
format datehosp date9. dob date9./
run/

## in 2015-16, California-Kaiser did not include test type in their dataset. 
They only use PCR, so if the patient was tested it was using PCR (testtype = 1) #/
data cak4b (keep=testedflu testtype testtype2 testresult testresult2 died hospid datehosp
dob year age state)/
set cak4 (rename= (FluTst=TestedFlu
##FluTstTyp1=TestType
FluTstTyp2=TestType2#/
FluResult1=TestResult
##FluResult2=TestResult2#/
Outcome=Died
hosp_tx=hospid
admdate=datehosp
))/
TestResult2=input(FluResult2,8.)/
year=year(datehosp)/
age = floor((intck("month", dob, datehosp,'C')) / 12)/
state = "CA2"/
if testedflu = 1 then do/
if testresult ne . then testtype=1/
if testresult2 ne . then testtype2=1/
end/
format datehosp date9. dob date9./
run/

data co4b (keep=testedflu testtype testtype2 testresult testresult2 died hospid datehosp
dob year age state)/
set co4 (rename= (FluTst=TestedFlu
FluTstTyp1=TestType
##FluTstTyp2=TestType2#/
FluResult1=TestResult
##FluResult2=TestResult2#/
Outcome=Died
##hosp_tx=hospid#/
admdate=datehosp))/
hospid=put(hosp_tx,5.)/
TestResult2=input(FluResult2,8.)/
TestType2=input(FluTstTyp2,8.)/
year=year(datehosp)/
age = floor((intck("month", dob, datehosp,'C')) / 12)/
state = "CO"/
run/

## admdates and dob not read in for 4 records 
will manually enter these #/

## check with CT about case hospitalized before DOB #/
data ct5b (keep=testedflu testtype testtype2 testresult testresult2 died hospid datehosp
dob year age state)/
set ct5 (rename= (FluTst=TestedFlu
FluTstTyp1=TestType
FluTstTyp2=TestType2
FluResult1=TestResult
FluResult2=TestResult2
Outcome=Died
hosp_tx=hospid
admdate=datehosp))/
if caseid=1 then do/
datehosp='16dec2015'd/
dob='19dec2011'd/
end/

if caseid=2 then do/
datehosp='29dec2015'd/
dob='18nov2011'd/
end/

if caseid=3 then do/
datehosp='22dec2015'd/
dob='10oct2011'd/
end/

## A O'Halloran 07/10/2018 
Confirmed with Kim that the DOB should be jan 22 1989#/
  if caseid=105 then do/
datehosp='23dec2015'd/
dob='22jan1989'd/
end/

## A O'Halloran 07/10/2018 
Confirmed with Kim that the DOB should be sep 20 1944#/
if caseid=213 then do/
dob='20sep1944'd/
end/

year=year(datehosp)/
age = floor((intck("month", dob, datehosp,'C')) / 12)/
state = "CT"/
run/

data ga3b (keep=testedflu testtype testtype2 testresult testresult2 died hospid datehosp
dob year age state)/
set ga3(rename= (FluTst=TestedFlu
FluTstTyp1=TestType
##FluTstTyp2=TestType2#/
FluResult1=TestResult
##FluResult2=TestResult2#/
Outcome=Died
##hosp_tx=hospid#/
admdate=datehosp))/
hospid=put(hosp_tx,5.)/
TestResult2=input(FluResult2,8.)/
TestType2=input(FluTstType2,8.)/
year=year(datehosp)/
age = floor((intck("month", dob, datehosp,'C')) / 12)/
state = "GA"/
run/

data mib (keep=testedflu testtype testtype2 testresult testresult2 died hospid datehosp
dob year age state)/
set mi(rename= (FluTst=TestedFlu
FluTstTyp1=TestType
FluTstTyp2=TestType2
FluResult1=TestResult
FluResult2=TestResult2
Outcome=Died
hosp_tx=hospid
dob=_dob
##admdate=datehosp#/))/
dob=input(_dob,mmddyy10.)/
datehosp=input(admdate,mmddyy10.)/
year=year(datehosp)/
age = floor((intck("month", dob, datehosp,'C')) / 12)/
state = "MI"/
run/

## check with MN about case MN151600290 born in 2020 #/
data mn3b (keep=testedflu testtype testtype2 testresult testresult2 died hospid datehosp
dob year age state)/
set mn3(rename= (FluTst=TestedFlu
FluTstTyp1=TestType
##FluTstTyp2=TestType2#/
FluResult1=TestResult
##FluResult2=TestResult2#/
Outcome=Died
hosp_tx=hospid
admdate=datehosp))/

TestType2=input(FluTstTyp2,8.)/
TestResult2=input(FluResult2,8.)/

## A O'Halloran 07/10/2018 
Confirmed with Melissa McMahon that the DOB should be sep 20 1944#/
  if caseid="MN151600290" then do/
dob='05nov1920'd/
end/

year=year(datehosp)/
age = floor((intck("month", dob, datehosp,'C')) / 12)/
state = "MN"/

run/

## A OHalloran 07/10/2018
nm4 read in 4 empty records - deleting these #/
  
  data nm4b (keep=testedflu testtype testtype2 testresult testresult2 died hospid datehosp
             dob year age state)/
set nm4(rename= (FluTst=TestedFlu
                 FluTstTyp1=TestType
                 FluTstTyp2=TestType2
                 FluResult1=TestResult
                 FluResult2=TestResult2
                 ##Outcome=Died#/
                   hosp_tx=hospid
                 admdate=datehosp))/
if caseid ne " "/
year=year(datehosp)/
age = floor((intck("month", dob, datehosp,'C')) / 12)/
state = "NM"/
run/

#proc contents data=ny4/
#run/

data ny4b (keep=testedflu testtype testtype2 testresult testresult2 died hospid datehosp
           dob year age state)/
set ny4(rename= (##FluTst=TestedFlu
                 FluTstTyp1=TestType
                 FluTstTyp2=TestType2
                 FluResult1=TestResult
                 FluResult2=TestResult2
                 Outcome=Died
                 hosp_tx=hospid#/
                   admdate=datehosp))/
if caseid ne " "/

hospid=put(hosp_tx,5.)/

if FluTst="1=Yes" then testedflu=1/
else if flutst="2=No" then testedflu=2/

if FluResult1="1=Positive" then TestResult=1/
else if FluResult1="2=Negative" then TestResult=2/

if FluResult2="1=Positive" then TestResult2=1/
else if FluResult2="2=Negative" then TestResult2=2/

if FluTstTyp1 = "1=Rapid Antigen" 	then TestType = 1/
if FluTstTyp2 = "1=Rapid Antigen" 	then TestType2 = 1/ 
if FluTstTyp1 = "2=Molecular Assay" 	then TestType = 2/
if FluTstTyp2 = "2=Molecular Assay" 	then TestType2 = 2/
if FluTstTyp1 = "3=Viral Culture" 		then TestType = 3/
if FluTstTyp2 = "3=Viral Culture" 		then TestType2 = 3/

if outcome="1=Deceased" then died=1/
else if outcome="2=Alive" then died=2/

year=year(datehosp)/
age = floor((intck("month", dob, datehosp,'C')) / 12)/
state = "NY"/
run/

data ny5b (keep=testedflu testtype testtype2 testresult testresult2 died hospid datehosp
           dob year age state)/
set ny5(rename= (FluTst=TestedFlu
                 FluTstTyp1=TestType
                 ##FluTstTyp2=TestType2#/
                   FluResult1=TestResult
                 ##FluResult2=TestResult2#/
                   Outcome=Died
                 hosp_tx=hospid
                 admdate=datehosp))/


if caseid ne " "/

TestType2=input(FluTstTyp2,8.)/
TestResult2=input(FluResult2,8.)/
year=year(datehosp)/
age = floor((intck("month", dob, datehosp,'C')) / 12)/
state = "NY"/
run/

data or5b(keep=testedflu testtype testtype2 testresult testresult2 died hospid datehosp
          dob year age state)/
set or5(rename= (FluTst=TestedFlu
                 FluTstTyp1=TestType
                 ##FluTstTyp2=TestType2#/
                   FluResult1=TestResult
                 ##FluResult2=TestResult2#/
                   Outcome=Died
                 hosp_tx=hospid
                 admdate=datehosp))/
if caseid ne " "/
TestType2=input(FluTstTyp2,8.)/
TestResult2=input(FluResult2,8.)/
year=year(datehosp)/
age = floor((intck("month", dob, datehosp,'C')) / 12)/
state = "OR"/
run/

## A O'Halloran 07/11/2018 
TN submitted new data file with dob and age added
using this file instead#/

##proc contents data=tn2/
run/ #/

##
data tn2b/
set tn2(rename= (FluTst=TestedFlu
FluTstTyp1=TestType
FluTstTyp2=TestType2
FluResult1=TestResult
FluResult2=TestResult2
Outcome=Died
hosp_tx=hospid
))/
# Per Danielle Ndi on 7/12/2018/
if caseid = "TNFB15160350" then admdate="4/19/2016"/ 


datehosp=input(admdate,mmddyy10.)/

if caseid = "TNFB15160159" then datehosp=37049/
if caseid ne " "/
year=year(datehosp)/
age = floor((intck("month", dob, datehosp,'C')) / 12)/
state = "TN"/
run/
#/


# ### Combining datasets from 2010-11, 2011-12, 2012-13, 2013-14, and 2014-15 seasons ### #/
## A OHalloran - 07/12/2018 - Adding 1516 data #/
data all00 / 

set 	##ALL0 ALL1#/
CA2B CAK1 CAK2 CAK3b CA3B CAK4B
CO2 CO3 CO4B
CT2 CT3 CT4 CT5B
GAb GA2b GA3B
MIB
MNb MN2 MN3B
NM2 NM3b NM4B
NY2 NY3 NY4B NY5B
OR2 OR3 OR4 OR5B
TNb/
month1516=month(datehosp)/ ## A OHalloran - month for 1516 data #/
month=upcase(month)/
year=year(DateHosp)/
run/

