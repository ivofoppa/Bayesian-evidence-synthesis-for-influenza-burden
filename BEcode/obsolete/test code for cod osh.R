ls <- unique(dataset[which(dataset$Died==1 & (is.na(dataset$ICD9_1) | dataset$ICD9_1!='') & dataset$TestResult==1),])
### only flu pos
ls <- unique(dataset[which(dataset$Died==1 & (is.na(dataset$ICD9_1) | dataset$ICD9_1!='') & dataset$TestResult==1),])
length(ls[,1])

ls2 <- ls$ICD9_1
ls3 <- sapply(as.vector(ls2),function(x) ifelse(is.character(x) & !is.na(x),substr(x,1,3),'0'),USE.NAMES = F)

ls4 <- sort(unique(ls3))

### ICD-9 codes:
# "038" Septicemia

# "162" "198" exclude: cancer

# "410" ischemic heart disease; "415" diseases of pulm. circulation; "424" other forms of heart disease; 
# "431" exclude? intracerebral hemorrhage

# "466" acute bronchitis/bronchiolitis 

# "480" viral pneumonia;  "481" pneumococcal pneumonia; "482" other bacterial pneumonia; 
# "483" pneumonia due to other specified organisms; "484" Pneumonia in infectious diseases classified elsewhere;
# "485" bronchopenumonia, "486" unspecified; pneumonia, unspecified; 

# "487" influenza; "488" avian influenza 
# "510" empyema; "511" pleurisy; "518" other diseases of the lung; 

# exclude "820" fracture of lower limb; "965" poisoning
