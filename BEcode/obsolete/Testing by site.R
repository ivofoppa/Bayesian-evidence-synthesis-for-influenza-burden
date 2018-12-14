bfolder <- "C:/Users/VOR1/Documents/GitHub/Bayesian-evidence-synthesis-for-influenza-burden/BEdata" ## project folder
setwd(paste0(bfolder,''))

load(file = "my_work_space.RData")

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
CAK3$testtype <- sapply(seq_along(CAK3[,1]), function(k) ifelse(CAK3$testedflu[k]==1,1,NA))

MN/ #(rename=
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
### GA2
GA2$State <- 'GA'
for (col in 1:length(GA2[1,])){
  GA2[,col] <- as.vector(GA2[,col])
}

for (k in 1:length(GA2$HospID)){
  if (GA2$Flu.Test.Type.1[k] =='Rapid') {GA2$TestType[k] <- 2
  } else if (GA2$Flu.Test.Type.1[k] =='RT-PCR') {GA2$TestType[k] <- 1
  } else if (GA2$Flu.Test.Type.1[k] =='PCR') {GA2$TestType[k] <- 1
  } else {GA2$TestType <- 9}

  if (is.na(FluTest.1..Yes.No.[k])){FluTest.1..Yes.No.[k] <- 9
  } else if (GA2$Flu.Test.Type.2[k] =='Rapid') {GA2$TestType2[k] <- 2
  } else if (GA2$Flu.Test.Type.2[k] =='RT-PCR') {GA2$TestType2[k] <- 1
  } else if (GA2$Flu.Test.Type.2[k] =='PCR') {GA2$TestType2[k] <- 1
  } else {GA2$TestType <- 9}


  if (is.na(FluTest.1..Yes.No.[k])){FluTest.1..Yes.No.[k] <- 9
  } else if (GA2$FluTest.1..Yes.No.[k] =='Yes') {GA2$testedflu[k] <- 1
 } else if (GA2$FluTest.1..Yes.No.[k] =='No') {GA2$testedflu[k] <- 2
 } else if (GA2$Flu.Test.Result.1[k] %in% c('Positive A', 'Positive B', 'Positive Unknown', 'Positive A H1N1') ) {GA2$TestResult[k] <- 1
 } else if (GA2$Flu.Test.Result.2[k] %in% c('Positive A', 'Positive B', 'Positive Unknown', 'Positive A H1N1')) {GA2$TestResult2[k] <- 1
 } else if (GA2$Flu.Test.Result.1[k] =='Negative') {GA2$TestResult[k] <- 2
 } else if (GA2$Flu.Test.Result.2[k] =='Negative') {GA2$TestResult[k] <- 2}

  if (is.na(GA2$ICU[k])){GA2$ICU[k] <- 9
  } else if (GA2$ICU[k] =='Yes'){GA2$ICU[k] <- 1
  } else if (GA2$ICU[k] =='No'){GA2$ICU[k] <- 2
  } else {GA2$ICU[k] <- 9}
}
if (__Died ==Yes)				 Died = 1/
if (__Died ==No)				 Died = 2/
if (__Died ==Unk"				 Died = 9/

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
if testedflu = 1  do/
if testresult ne .  testtype=1/
if testresult2 ne .  testtype2=1/
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
if caseid=1  do/
datehosp='16dec2015'd/
dob='19dec2011'd/
end/

if caseid=2  do/
datehosp='29dec2015'd/
dob='18nov2011'd/
end/

if caseid=3  do/
datehosp='22dec2015'd/
dob='10oct2011'd/
end/

## A O'Halloran 07/10/2018
Confirmed with Kim that the DOB should be jan 22 1989#/
  if caseid=105  do/
datehosp='23dec2015'd/
dob='22jan1989'd/
end/

## A O'Halloran 07/10/2018
Confirmed with Kim that the DOB should be sep 20 1944#/
if caseid=213  do/
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
  if caseid="MN151600290"  do/
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

if FluTst="1=Yes"  testedflu=1/
else if flutst="2=No"  testedflu=2/

if FluResult1="1=Positive"  TestResult=1/
else if FluResult1="2=Negative"  TestResult=2/

if FluResult2="1=Positive"  TestResult2=1/
else if FluResult2="2=Negative"  TestResult2=2/

if FluTstTyp1 = "1=Rapid Antigen" 	 TestType = 1/
if FluTstTyp2 = "1=Rapid Antigen" 	 TestType2 = 1/
if FluTstTyp1 = "2=Molecular Assay" 	 TestType = 2/
if FluTstTyp2 = "2=Molecular Assay" 	 TestType2 = 2/
if FluTstTyp1 = "3=Viral Culture" 		 TestType = 3/
if FluTstTyp2 = "3=Viral Culture" 		 TestType2 = 3/

if outcome="1=Deceased"  died=1/
else if outcome="2=Alive"  died=2/

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
if caseid = "TNFB15160350"  admdate="4/19/2016"/


datehosp=input(admdate,mmddyy10.)/

if caseid = "TNFB15160159"  datehosp=37049/
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
