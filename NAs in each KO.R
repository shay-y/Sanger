## read 0\1 data:
zoo <- read.csv(file = "Oct29_2014_Xray_B6N_MGPSelect_Oor1_cleaned.csv") # 8176,79 

## 46 variables of interest:
XrayVariableList <- c("Number.Of.Thoracic.Vertebrae","Number.Of.Lumbar.Vertebrae","Number.Of.Pelvic.Vertebrae","Number.Of.Caudal.Vertebrae", "Transitional.Vertebrae" , "Shape.Of.Vertebrae","Fusion.Of.Vertebrae", "Processes.On.Vertebrae","Maxilla","Zygomatic.Bone","Number.Of.Cervical.Vertebrae","Skull.Shape","Number.Of.Ribs.Right","Number.Of.Ribs.Left","Shape.Of.Ribcage","Shape.Of.Ribs","Rib.Fusions","Clavicle","Scapula"  ,"Humerus","Radius","Ulna","Pelvis","Femur","Tibia","Fibula","Joints","Shape.Of.Spine","Teeth","Mandible","Number.Of.Digits","Digit.Integrity","Syndactylism","Polysyndactylism","Brachydactylism","Kyphosis","Lordosis","Scoliosis","Spinous.Processes","Transverse.Processes","Fusion.Processes","Caudal.Processes","Cervical.Processes","Lumbar.Processes","Sacral.Processes","Thoracic.Processes") 

## remove 3 variables with no values:
XrayVariableList_c <- XrayVariableList[!(XrayVariableList %in% c("Spinous.Processes","Transverse.Processes","Processes.On.Vertebrae"))]

## create data with variables of interest only:
zoo1 <- zoo[c("Colony.Prefix","Genotype2","Genotype",XrayVariableList_c)]

sum(is.na(zoo1[which(apply(zoo1,1,function(x) sum(is.na(x)))>14),]))
sum(is.na(zoo1))
## summarize number of NA occurences in each KO group and each outcome:
n_na <- aggregate(. ~ Genotype2 + Colony.Prefix,data = zoo1,subset = Genotype2 !="WT",
                  function(x) sum(is.na(x)),na.action=NULL)[-(1:3)]

## NAs in each KO:
rowSums(n_na)

## NAs in each outcome:
colSums(n_na)

## number of NA occurences in WT group and each outcome:
n_na <- aggregate(. ~ 1,data = zoo1,subset = Genotype2 =="WT",
                  function(x) sum(is.na(x)),na.action=NULL)[-(1:3)]
n_na
