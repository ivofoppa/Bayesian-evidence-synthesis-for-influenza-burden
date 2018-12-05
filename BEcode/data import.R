ALL0 <- data.frame(read_excel('FULL_EIPBURDEN_DATE_v2.xls', sheet='PartA'))
ALL1 <- data.frame(read_excel('FULL_EIPBURDEN_201112.xls', sheet='PartA'))

CA2 <- data.frame(read.xlsx('CA_InfluenzaBurden_012916.xlsx', sheetIndex=1))
CAK1 <- data.frame(read.xlsx('CA_Kaiser_Flu Burden_2012_13.xlsx', sheetIndex=1))
CAK2 <- data.frame(read.xlsx('CA_Kaiser_Flu Burden_2013_14.xlsx', sheetIndex=1))
CAK3 <- data.frame(read.xlsx('CA_Kaiser flu burden_2014_15.xlsx', sheetIndex = 1))
CAK4 <- data.frame(read.xlsx('Kaiser_1516_burden_to CDC.xlsx', sheetIndex=1))
CA3 <- data.frame(read.xlsx('CA_1516_BurdenDatabase_Final.xlsx', sheetIndex=1))

#### Colorado ### 
#### 2012-13 & 2013-14 season ####
CO2 <- data.frame(read.xlsx('CO_BURDEN_12032014.xlsx', sheetIndex=3))
#### 2014-15 season ####
CO3 <- data.frame(read.xlsx('CO_BURDEN_201415.xlsx', sheetIndex=3))
#### 2015-16 season ####
CO4 <- data.frame(read.xlsx('CO_Burden_1516_10032017.xlsx', sheetIndex=3))
#### Connecticut ####
#### 2012-13 season ####
CT2 <- data.frame(read.xlsx('CT_BURDEN_04272015.xlsx', sheetIndex=3))
#### 2013-14 season ####
CT3 <- data.frame(read.xlsx('CT_EIPBURDEN_12072015.xlsx', sheetIndex=3))
## CT
#### 2014-15 season ####
CT4 <- data.frame(read.xlsx('CT_BURDEN_1415_11222016.xlsx', sheetIndex=2))
#### 2015-16 season ####
CT5 <- data.frame(read.xlsx('CT_BURDEN_1516_10312017.xlsx', sheetIndex=3))
#### Georgia ####
#### 2010-11 & 2011-12 & 2012-13 seasons ####
GA <- data.frame(read.xlsx('GA_EIPBURDEN_04052016.xlsx', sheetIndex=4))
#### 2013-14 & 2014-15 seasons ####
GA2 <- data.frame(read.xlsx('GA Influenza Burden 2013-14_2014-15 CDC Submission 9-14-2017.xlsx', sheetIndex=1))
#### 2015-16 seasons ####
GA3 <- data.frame(read_excel('FluBurdenSample_2015_2016_CDC Submission.xls', sheet=' Final CDC Format'))

#### Michigan ####
#### 2015-16 season ####
MI <- data.frame(read.xlsx('MI_BURDEN_1516_06132018.xlsx', sheetIndex=3))

#### Minnesota ####
#### 2012-13 & 2013-14 season ####
MN <- data.frame(read.xlsx('MN_BURDEN_03092016.xlsx', sheetIndex=3))
#### 2014-15 season ####
MN2 <- data.frame(read.xlsx('STATE_EIPBURDEN_DATE_v2014_MN_1415Data.xlsx', sheetIndex=3))
#### 2015-16 season ####
MN3 <- data.frame(read.xlsx('MN_BURDEN_1516_12192017.xlsx', sheetIndex=3))

#### New Mexico ####
#### 2012-13 & 2013-14 season ####
NM2 <- data.frame(read.xlsx('NM_EIPBURDEN_2015.05.05.xlsx', sheetIndex=1))
#### 2014-15 season ####
NM3 <- data.frame(read.xlsx('2014-15 EIPBURDEN NM Final.xlsx', sheetIndex=1))
#### 2015-16 season ####
NM4 <- data.frame(read_excel('NM_Burden_1516_12212017.xls', sheet=' Patient data'))

#### New York ####
#### 2012-13 & 2013-14 season ####
NY2 <- data.frame(read.xlsx('STATE_EIPBURDEN_DATE_v2014_9.23.14.xlsx', sheetIndex=4))
#### 2014-15 season ####
NY3 <- data.frame(read.xlsx('STATE_EIPBURDEN_NYR_2014-2015.xlsx', sheetIndex=3))
#### 2015-16 season ####
## NYR #
NY4 <- data.frame(read.xlsx('NY_BURDEN_1516_03012018.xlsx', sheetIndex=1))
## NYA #
NY5 <- data.frame(read.xlsx('NY_BURDEN_1516_06072018.xlsx', sheetIndex=1))

#### Oregon ####
#### 2012-13 season ####
OR2 <- data.frame(read.xlsx('OR_Burden_12122014.xlsx', sheetIndex=1))
#### 2013-14 season ####
OR3 <- data.frame(read.xlsx('OR_EIPBURDEN_08242015_v2014.xlsx', sheetIndex=3))
#### 2014-15 season ####
OR4 <- data.frame(read.xlsx('STATE_EIPBURDEN_DATE_OR_v2014.xlsx', sheetIndex=3))
#### 2015-16 season ####
OR5 <- data.frame(read.xlsx('OR_BURDEN_1516_04102018.xlsx', sheetIndex=3))

#### Tennessee ####
#### 2012-13 & 2013-14 & 2014-15 season ####
TN <- data.frame(read.xlsx('TN_EIPBURDEN_20170112.xlsx', sheetIndex=3))
#### 2015-16 season ####
TN2 <- data.frame(read.xlsx('TN_EIPBURDEN_20180710.xlsx', sheetIndex=1))
