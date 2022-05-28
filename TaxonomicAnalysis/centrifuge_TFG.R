# TFG centrifuge analysis

setwd("/media/sequentia/visitors/visitor8/TFG")

library(phyloseq)
library(vegan)

# read files


SRR15853012 <- read.delim("read-based/centrifuge_report/SRR15853012-report.txt", 
                          stringsAsFactors = F, sep="\t")
SRR15852941 <- read.delim("read-based/centrifuge_report/SRR15852941-report.txt", 
                          stringsAsFactors = F, sep="\t")
SRR15852948 <- read.delim("read-based/centrifuge_report/SRR15852948-report.txt", 
                         stringsAsFactors = F, sep="\t")
SRR15853009 <- read.delim("read-based/centrifuge_report/SRR15853009-report.txt", 
                          stringsAsFactors = F, sep="\t")
SRR15853002 <- read.delim("read-based/centrifuge_report/SRR15853002-report.txt", 
                          stringsAsFactors = F, sep="\t")
SRR15852986 <- read.delim("read-based/centrifuge_report/SRR15852986-report.txt", 
                          stringsAsFactors = F, sep="\t")
SRR15852987 <- read.delim("read-based/centrifuge_report/SRR15852987-report.txt", 
                          stringsAsFactors = F, sep="\t")
SRR15852998 <- read.delim("read-based/centrifuge_report/SRR15852998-report.txt", 
                         stringsAsFactors = F, sep="\t")
SRR15853011 <- read.delim("read-based/centrifuge_report/SRR15853011-report.txt", 
                                   stringsAsFactors = F, sep="\t")
SRR15853000 <- read.delim("read-based/centrifuge_report/SRR15853000-report.txt", 
                          stringsAsFactors = F, sep="\t")
SRR15853017 <- read.delim("read-based/centrifuge_report/SRR15853017-report.txt", 
                          stringsAsFactors = F, sep="\t")
SRR15853005 <- read.delim("read-based/centrifuge_report/SRR15853005-report.txt", 
                          stringsAsFactors = F, sep="\t")
SRR15852988 <- read.delim("read-based/centrifuge_report/SRR15852988-report.txt", 
                          stringsAsFactors = F, sep="\t")
SRR15852994 <- read.delim("read-based/centrifuge_report/SRR15852994-report.txt", 
                          stringsAsFactors = F, sep="\t")
SRR15852995 <- read.delim("read-based/centrifuge_report/SRR15852995-report.txt", 
                          stringsAsFactors = F, sep="\t") 
SRR15852997 <- read.delim("read-based/centrifuge_report/SRR15852997-report.txt", 
                          stringsAsFactors = F, sep="\t")  
SRR15852989 <- read.delim("read-based/centrifuge_report/SRR15852989-report.txt", 
                          stringsAsFactors = F, sep="\t") 
SRR15852990 <- read.delim("read-based/centrifuge_report/SRR15852990-report.txt", 
                          stringsAsFactors = F, sep="\t")
SRR15853001 <- read.delim("read-based/centrifuge_report/SRR15853001-report.txt", 
                          stringsAsFactors = F, sep="\t") 
SRR15852992 <- read.delim("read-based/centrifuge_report/SRR15852992-report.txt", 
                          stringsAsFactors = F, sep="\t")
SRR15852951 <- read.delim("read-based/centrifuge_report/SRR15852951-report.txt", 
                          stringsAsFactors = F, sep="\t")
SRR15852922 <- read.delim("read-based/centrifuge_report/SRR15852922-report.txt", 
                          stringsAsFactors = F, sep="\t")
SRR15852955 <- read.delim("read-based/centrifuge_report/SRR15852955-report.txt", 
                          stringsAsFactors = F, sep="\t") 
SRR15852956 <- read.delim("read-based/centrifuge_report/SRR15852956-report.txt", 
                          stringsAsFactors = F, sep="\t") 
SRR15852962 <- read.delim("read-based/centrifuge_report/SRR15852962-report.txt", 
                          stringsAsFactors = F, sep="\t")
SRR15852963 <- read.delim("read-based/centrifuge_report/SRR15852963-report.txt", 
                          stringsAsFactors = F, sep="\t") 
SRR15852982 <- read.delim("read-based/centrifuge_report/SRR15852982-report.txt", 
                          stringsAsFactors = F, sep="\t") 
SRR15853026 <- read.delim("read-based/centrifuge_report/SRR15853026-report.txt", 
                          stringsAsFactors = F, sep="\t") 
SRR15853027 <- read.delim("read-based/centrifuge_report/SRR15853027-report.txt", 
                          stringsAsFactors = F, sep="\t")
SRR15852932 <- read.delim("read-based/centrifuge_report/SRR15852932-report.txt", 
                          stringsAsFactors = F, sep="\t")
SRR15852950 <- read.delim("read-based/centrifuge_report/SRR15852950-report.txt", 
                          stringsAsFactors = F, sep="\t") 
SRR15852952 <- read.delim("read-based/centrifuge_report/SRR15852952-report.txt", 
                          stringsAsFactors = F, sep="\t")
SRR15852954 <- read.delim("read-based/centrifuge_report/SRR15852954-report.txt", 
                                 stringsAsFactors = F, sep="\t")
SRR15852957 <- read.delim("read-based/centrifuge_report/SRR15852957-report.txt", 
                          stringsAsFactors = F, sep="\t")
SRR15852960 <- read.delim("read-based/centrifuge_report/SRR15852960-report.txt", 
                          stringsAsFactors = F, sep="\t")
SRR15852993 <-read.delim("read-based/centrifuge_report/SRR15852993-report.txt", 
                         stringsAsFactors = F, sep="\t")
SRR15852934 <- read.delim("read-based/centrifuge_report/SRR15852934-report.txt", 
                          stringsAsFactors = F, sep="\t")
SRR15852968 <- read.delim("read-based/centrifuge_report/SRR15852968-report.txt", 
                          stringsAsFactors = F, sep="\t")
SRR15852969 <- read.delim("read-based/centrifuge_report/SRR15852969-report.txt", 
                          stringsAsFactors = F, sep="\t") 
SRR15853028 <- read.delim("read-based/centrifuge_report/SRR15853028-report.txt", 
                          stringsAsFactors = F, sep="\t")

# selection of columns V1, V3 and V5


SRR15853012 <- SRR15853012[, c("name", "taxRank","numReads")]  
SRR15852941 <- SRR15852941[, c("name", "taxRank","numReads")]
SRR15852948 <- SRR15852948[, c("name", "taxRank","numReads")]
SRR15853009 <- SRR15853009[, c("name", "taxRank","numReads")]
SRR15853002 <- SRR15853002[, c("name", "taxRank","numReads")]
SRR15852986 <- SRR15852986[, c("name", "taxRank","numReads")]
SRR15852987 <- SRR15852987[, c("name", "taxRank","numReads")]
SRR15852998 <- SRR15852998[, c("name", "taxRank","numReads")]
SRR15853011 <- SRR15853011[, c("name", "taxRank","numReads")]
SRR15853000 <- SRR15853000[, c("name", "taxRank","numReads")]
SRR15853017 <- SRR15853017[, c("name", "taxRank","numReads")]
SRR15853005 <- SRR15853005[, c("name", "taxRank","numReads")]
SRR15852988 <- SRR15852988[, c("name", "taxRank","numReads")]
SRR15852994 <- SRR15852994[, c("name", "taxRank","numReads")]
SRR15852995 <- SRR15852995[, c("name", "taxRank","numReads")]  
SRR15852997 <- SRR15852997[, c("name", "taxRank","numReads")]
SRR15852989 <- SRR15852989[, c("name", "taxRank","numReads")] 
SRR15852990 <- SRR15852990[, c("name", "taxRank","numReads")]
SRR15853001 <- SRR15853001[, c("name", "taxRank","numReads")] 
SRR15852992 <- SRR15852992[, c("name", "taxRank","numReads")]
SRR15852951 <- SRR15852951[, c("name", "taxRank","numReads")] 
SRR15852922 <- SRR15852922[, c("name", "taxRank","numReads")] 
SRR15852955 <- SRR15852955[, c("name", "taxRank","numReads")] 
SRR15852956 <- SRR15852956[, c("name", "taxRank","numReads")] 
SRR15852962 <- SRR15852962[, c("name", "taxRank","numReads")] 
SRR15852963 <- SRR15852963[, c("name", "taxRank","numReads")] 
SRR15852982 <- SRR15852982[, c("name", "taxRank","numReads")] 
SRR15853026 <- SRR15853026[, c("name", "taxRank","numReads")] 
SRR15853027 <- SRR15853027[, c("name", "taxRank","numReads")] 
SRR15852932 <- SRR15852932[, c("name", "taxRank","numReads")] 
SRR15852950 <- SRR15852950[, c("name", "taxRank","numReads")] 
SRR15852952 <- SRR15852952[, c("name", "taxRank","numReads")] 
SRR15852954 <- SRR15852954[, c("name", "taxRank","numReads")] 
SRR15852957 <- SRR15852957[, c("name", "taxRank","numReads")] 
SRR15852960 <- SRR15852960[, c("name", "taxRank","numReads")] 
SRR15852993 <- SRR15852993[, c("name", "taxRank","numReads")]  
SRR15852934 <- SRR15852934[, c("name", "taxRank","numReads")] 
SRR15852968 <- SRR15852968[, c("name", "taxRank","numReads")] 
SRR15852969 <- SRR15852969[, c("name", "taxRank","numReads")]      
SRR15853028 <- SRR15853028[, c("name", "taxRank","numReads")] 

# want to work with species level

SRR15853012 <- SRR15853012[grep("species", SRR15853012$taxRank),]
SRR15852941 <- SRR15852941[grep("species", SRR15852941$taxRank),]
SRR15852948 <- SRR15852948[grep("species", SRR15852948$taxRank),]
SRR15853009 <- SRR15853009[grep("species", SRR15853009$taxRank),]
SRR15853002 <- SRR15853002[grep("species", SRR15853002$taxRank),]
SRR15852986 <- SRR15852986[grep("species", SRR15852986$taxRank),]
SRR15852987 <- SRR15852987[grep("species", SRR15852987$taxRank),]
SRR15852998 <- SRR15852998[grep("species", SRR15852998$taxRank),]
SRR15853011 <- SRR15853011[grep("species", SRR15853011$taxRank),]
SRR15853000 <- SRR15853000[grep("species", SRR15853000$taxRank),]
SRR15853017 <- SRR15853017[grep("species", SRR15853017$taxRank),]
SRR15853005 <- SRR15853005[grep("species", SRR15853005$taxRank),]
SRR15852988 <- SRR15852988[grep("species", SRR15852988$taxRank),]
SRR15852994 <- SRR15852994[grep("species", SRR15852994$taxRank),]
SRR15852995 <- SRR15852995[grep("species", SRR15852995$taxRank),]
SRR15852997 <- SRR15852997[grep("species", SRR15852997$taxRank),]  
SRR15852989 <- SRR15852989[grep("species", SRR15852989$taxRank),] 
SRR15852990 <- SRR15852990[grep("species", SRR15852990$taxRank),]
SRR15853001 <- SRR15853001[grep("species", SRR15853001$taxRank),] 
SRR15852992 <- SRR15852992[grep("species", SRR15852992$taxRank),]
SRR15852951 <- SRR15852951[grep("species", SRR15852951$taxRank),] 
SRR15852922 <- SRR15852922[grep("species", SRR15852922$taxRank),] 
SRR15852955 <- SRR15852955[grep("species", SRR15852955$taxRank),] 
SRR15852956 <- SRR15852956[grep("species", SRR15852956$taxRank),] 
SRR15852962 <- SRR15852962[grep("species", SRR15852962$taxRank),] 
SRR15852963 <- SRR15852963[grep("species", SRR15852963$taxRank),] 
SRR15852982 <- SRR15852982[grep("species", SRR15852982$taxRank),] 
SRR15853026 <- SRR15853026[grep("species", SRR15853026$taxRank),] 
SRR15853027 <- SRR15853027[grep("species", SRR15853027$taxRank),] 
SRR15852932 <- SRR15852932[grep("species", SRR15852932$taxRank),] 
SRR15852950 <- SRR15852950[grep("species", SRR15852950$taxRank),] 
SRR15852952 <- SRR15852952[grep("species", SRR15852952$taxRank),] 
SRR15852954 <- SRR15852954[grep("species", SRR15852954$taxRank),] 
SRR15852957 <- SRR15852957[grep("species", SRR15852957$taxRank),] 
SRR15852960 <- SRR15852960[grep("species", SRR15852960$taxRank),] 
SRR15852993 <- SRR15852993[grep("species", SRR15852993$taxRank),] 
SRR15852934 <- SRR15852934[grep("species", SRR15852934$taxRank),] 
SRR15852968 <- SRR15852968[grep("species", SRR15852968$taxRank),] 
SRR15852969 <- SRR15852969[grep("species", SRR15852969$taxRank),]      
SRR15853028 <- SRR15853028[grep("species", SRR15853028$taxRank),] 

# Drop columns whose name starts with "tax
SRR15853012 <- SRR15853012[,!grepl("^tax",names(SRR15853012))]
SRR15852941 <- SRR15852941[,!grepl("^tax",names(SRR15852941))]
SRR15852948 <- SRR15852948[,!grepl("^tax",names(SRR15852948))]
SRR15853009 <- SRR15853009[,!grepl("^tax",names(SRR15853009))]
SRR15853002 <- SRR15853002[,!grepl("^tax",names(SRR15853002))]
SRR15852986 <- SRR15852986[,!grepl("^tax",names(SRR15852986))]
SRR15852987 <- SRR15852987[,!grepl("^tax",names(SRR15852987))]
SRR15852998 <- SRR15852998[,!grepl("^tax",names(SRR15852998))]
SRR15853011 <- SRR15853011[,!grepl("^tax",names(SRR15853011))]
SRR15853000 <- SRR15853000[,!grepl("^tax",names(SRR15853000))]
SRR15853017 <- SRR15853017[,!grepl("^tax",names(SRR15853017))]
SRR15853005 <- SRR15853005[,!grepl("^tax",names(SRR15853005))]
SRR15852988 <- SRR15852988[,!grepl("^tax",names(SRR15852988))]
SRR15852994 <- SRR15852994[,!grepl("^tax",names(SRR15852994))]
SRR15852995 <- SRR15852995[,!grepl("^tax",names(SRR15852995))]
SRR15852997 <- SRR15852997[,!grepl("^tax",names(SRR15852997))]  
SRR15852989 <- SRR15852989[,!grepl("^tax",names(SRR15852989))] 
SRR15852990 <- SRR15852990[,!grepl("^tax",names(SRR15852990))]
SRR15853001 <- SRR15853001[,!grepl("^tax",names(SRR15853001))] 
SRR15852992 <- SRR15852992[,!grepl("^tax",names(SRR15852992))]
SRR15852951 <- SRR15852951[,!grepl("^tax",names(SRR15852951))] 
SRR15852922 <- SRR15852922[,!grepl("^tax",names(SRR15852922))] 
SRR15852955 <- SRR15852955[,!grepl("^tax",names(SRR15852955))] 
SRR15852956 <- SRR15852956[,!grepl("^tax",names(SRR15852956))] 
SRR15852962 <- SRR15852962[,!grepl("^tax",names(SRR15852962))] 
SRR15852963 <- SRR15852963[,!grepl("^tax",names(SRR15852963))] 
SRR15852982 <- SRR15852982[,!grepl("^tax",names(SRR15852982))] 
SRR15853026 <- SRR15853026[,!grepl("^tax",names(SRR15853026))] 
SRR15853027 <- SRR15853027[,!grepl("^tax",names(SRR15853027))] 
SRR15852932 <- SRR15852932[,!grepl("^tax",names(SRR15852932))] 
SRR15852950 <- SRR15852950[,!grepl("^tax",names(SRR15852950))] 
SRR15852952 <- SRR15852952[,!grepl("^tax",names(SRR15852952))] 
SRR15852954 <- SRR15852954[,!grepl("^tax",names(SRR15852954))] 
SRR15852957 <- SRR15852957[,!grepl("^tax",names(SRR15852957))] 
SRR15852960 <- SRR15852960[,!grepl("^tax",names(SRR15852960))] 
SRR15852993 <- SRR15852993[,!grepl("^tax",names(SRR15852993))] 
SRR15852934 <- SRR15852934[,!grepl("^tax",names(SRR15852934))] 
SRR15852968 <- SRR15852968[,!grepl("^tax",names(SRR15852968))] 
SRR15852969 <- SRR15852969[,!grepl("^tax",names(SRR15852969))]      
SRR15853028 <- SRR15853028[,!grepl("^tax",names(SRR15853028))] 



# change the read column name by their accesion number

colnames(SRR15853012)[2] <- "SRR15853012"
colnames(SRR15852941)[2] <- "SRR15852941"
colnames(SRR15852948)[2] <- "SRR15852948"
colnames(SRR15853009)[2] <- "SRR15853009"
colnames(SRR15853002)[2] <- "SRR15853002"
colnames(SRR15852986)[2] <- "SRR15852986"
colnames(SRR15852987)[2] <- "SRR15852987"
colnames(SRR15852998)[2] <- "SRR15852998"
colnames(SRR15853011)[2] <- "SRR15853011"
colnames(SRR15853000)[2] <- "SRR15853000"
colnames(SRR15853017)[2] <- "SRR15853017"
colnames(SRR15853005)[2] <- "SRR15853005"
colnames(SRR15852988)[2] <- "SRR15852988"
colnames(SRR15852994)[2] <- "SRR15852994"
colnames(SRR15852995)[2] <- "SRR15852995"
colnames(SRR15852997)[2] <- "SRR15852997"
colnames(SRR15852989)[2] <- "SRR15852989"
colnames(SRR15852990)[2] <- "SRR15852990"
colnames(SRR15853001)[2] <- "SRR15853001"
colnames(SRR15852992)[2] <- "SRR15852992"
colnames(SRR15852951)[2] <- "SRR15852951"
colnames(SRR15852922)[2] <- "SRR15852922"
colnames(SRR15852955)[2] <- "SRR15852955"
colnames(SRR15852956)[2] <- "SRR15852956"
colnames(SRR15852962)[2] <- "SRR15852962"
colnames(SRR15852963)[2] <- "SRR15852963"
colnames(SRR15852982)[2] <- "SRR15852982"
colnames(SRR15853026)[2] <- "SRR15853026"
colnames(SRR15853027)[2] <- "SRR15853027"
colnames(SRR15852932)[2] <- "SRR15852932"
colnames(SRR15852950)[2] <- "SRR15852950"
colnames(SRR15852952)[2] <- "SRR15852952"
colnames(SRR15852954)[2] <- "SRR15852954"
colnames(SRR15852957)[2] <- "SRR15852957"
colnames(SRR15852960)[2] <- "SRR15852960"
colnames(SRR15852993)[2] <- "SRR15852993"
colnames(SRR15852934)[2] <- "SRR15852934"
colnames(SRR15852968)[2] <- "SRR15852968"
colnames(SRR15852969)[2] <- "SRR15852969"
colnames(SRR15853028)[2] <- "SRR15853028"

library(data.table) 

# merge the files
merged3 <- merge(as.data.table(SRR15853012), as.data.table(SRR15852941), by = "name", all = T)
merged3 <- merge(merged3, as.data.table(SRR15852948), by = "name", all = T)
merged3 <- merge(merged3, as.data.table(SRR15853009), by = "name", all = T)
merged3 <- merge(merged3, as.data.table(SRR15853002), by = "name", all = T)
merged3 <- merge(merged3, as.data.table(SRR15852986), by = "name", all = T)
merged3 <- merge(merged3, as.data.table(SRR15852987), by = "name", all = T)
merged3 <- merge(merged3, as.data.table(SRR15852998), by = "name", all = T)
merged3 <- merge(merged3, as.data.table(SRR15853011), by = "name", all = T)
merged3 <- merge(merged3,as.data.table(SRR15853000), by = "name", all = T)
merged3 <- merge(merged3,as.data.table(SRR15853017), by = "name", all = T)
merged3 <- merge(merged3,as.data.table(SRR15853005), by = "name", all = T)
merged3 <- merge(merged3,as.data.table(SRR15852988), by = "name", all = T)
merged3 <- merge(merged3,as.data.table(SRR15852994), by = "name", all = T)
merged3 <- merge(merged3,as.data.table(SRR15852995), by = "name", all = T)
merged3 <- merge(merged3,as.data.table(SRR15852997), by = "name", all = T)   
merged3 <- merge(merged3,as.data.table(SRR15852989), by = "name", all = T)
merged3 <- merge(merged3,as.data.table(SRR15852990), by = "name", all = T)
merged3 <- merge(merged3,as.data.table(SRR15853001), by = "name", all = T)  
merged3 <- merge(merged3,as.data.table(SRR15852992), by = "name", all = T)
merged3 <- merge(merged3,as.data.table(SRR15852951), by = "name", all = T)
merged3 <- merge(merged3,as.data.table(SRR15852922), by = "name", all = T)
merged3 <- merge(merged3,as.data.table(SRR15852955), by = "name", all = T)
merged3 <- merge(merged3,as.data.table(SRR15852956), by = "name", all = T)  
merged3 <- merge(merged3,as.data.table(SRR15852962), by = "name", all = T)  
merged3 <- merge(merged3,as.data.table(SRR15852963), by = "name", all = T)  
merged3 <- merge(merged3,as.data.table(SRR15852982), by = "name", all = T)  
merged3 <- merge(merged3,as.data.table(SRR15853026), by = "name", all = T)  
merged3 <- merge(merged3,as.data.table(SRR15853027), by = "name", all = T)  
merged3 <- merge(merged3,as.data.table(SRR15852932), by = "name", all = T)  
merged3 <- merge(merged3,as.data.table(SRR15852950), by = "name", all = T)  
merged3 <- merge(merged3,as.data.table(SRR15852952), by = "name", all = T)  
merged3 <- merge(merged3,as.data.table(SRR15852954), by = "name", all = T)  
merged3 <- merge(merged3,as.data.table(SRR15852957), by = "name", all = T)  
merged3 <- merge(merged3,as.data.table(SRR15852960), by = "name", all = T)  
merged3 <- merge(merged3,as.data.table(SRR15852993), by = "name", all = T)  
merged3 <- merge(merged3,as.data.table(SRR15852934), by = "name", all = T)  
merged3 <- merge(merged3,as.data.table(SRR15852968), by = "name", all = T)  
merged3 <- merge(merged3,as.data.table(SRR15852969), by = "name", all = T)       
merged3 <- merge(merged3,as.data.table(SRR15853028), by = "name", all = T)
                 

merged <- merged3

# replace NA by 0's
merged[is.na(merged)] <- 0

# convert from data table to data frame
library("data.table")
setDF(merged)

# remove Homo sapiens, as it is a host contamination 
merged <-merged[!grepl("Homo sapiens", merged$name),]

# species as row names
rownames(merged) <- merged$name

# remove first column, as it is a duplicate
merged <- merged[,2:41] # we have 40 samples


# select PD and CONTROL samples


general <- merged[, c("SRR15853012", "SRR15852941", "SRR15852948","SRR15853009", "SRR15853002",
                      "SRR15852986", "SRR15852987", "SRR15852998", "SRR15853011","SRR15853000",
                      "SRR15853017","SRR15853005", "SRR15852988", "SRR15852994","SRR15852995", 
                      "SRR15852997", "SRR15852989","SRR15852990","SRR15853001","SRR15852992",
                      "SRR15852951", "SRR15852922", "SRR15852955", "SRR15852956","SRR15852962",
                      "SRR15852963", "SRR15852982","SRR15853026","SRR15853027","SRR15852932",
                      "SRR15852950", "SRR15852952","SRR15852954","SRR15852957","SRR15852960",
                      "SRR15852993", "SRR15852934","SRR15852968","SRR15852969","SRR15853028")]
                      
PD <- merged[, c("SRR15853012", "SRR15852941", "SRR15852948","SRR15853009", "SRR15853002",
             "SRR15852986", "SRR15852987", "SRR15852998", "SRR15853011","SRR15853000",
             "SRR15853017","SRR15853005", "SRR15852988", "SRR15852994","SRR15852995", 
             "SRR15852997", "SRR15852989","SRR15852990","SRR15853001","SRR15852992")]


control <- merged[, c("SRR15852951", "SRR15852922", "SRR15852955", "SRR15852956","SRR15852962",
                      "SRR15852963", "SRR15852982","SRR15853026","SRR15853027","SRR15852932",
                      "SRR15852950", "SRR15852952","SRR15852954","SRR15852957","SRR15852960",
                      "SRR15852993", "SRR15852934","SRR15852968","SRR15852969","SRR15853028")]



# Filter species with only one read support (removing FP)
PD <- PD[rowSums(PD) > 1,]
control <- control[rowSums(control) >1,]
general <- general[rowSums(general) > 1,]



# compute abundances
species_percentage_PD <- t(t(PD)/colSums(PD))*100
species_percentage_control <- t(t(control)/colSums(control))*100
species_percentage_general <- t(t(general)/colSums(general))*100

## Taxonomic composition
# remove species with low abundance, below 1%
species_percentage_filter_PD <- species_percentage_PD[!apply(species_percentage_PD,1,function(x){sum(x<1)}) ==
                                                  20,] # 20 samples

species_percentage_filter_control <- species_percentage_control[!apply(species_percentage_control,1,function(x){sum(x<1)}) ==
                                                        20,] # 20 samples

species_percentage_filter_general <- species_percentage_general[!apply(species_percentage_general,1,function(x){sum(x<1)}) ==
                                                                  40,] # 40 samples
# Summarizing % species removed
species_percentage_filter_PD <- rbind(species_percentage_filter_PD, 100-colSums(species_percentage_filter_PD))
rownames(species_percentage_filter_PD)[nrow(species_percentage_filter_PD)] <- "Other"

species_percentage_filter_control <- rbind(species_percentage_filter_control, 100-colSums(species_percentage_filter_control))
rownames(species_percentage_filter_control)[nrow(species_percentage_filter_control)] <- "Other"

species_percentage_filter_general <- rbind(species_percentage_filter_general, 100-colSums(species_percentage_filter_general))
rownames(species_percentage_filter_general)[nrow(species_percentage_filter_general)] <- "Other"

# save merge data.frame
save(merged,file="merged.Rda")
load("data.Rda")

## Diversity analysis
# create phyloseq object
species_phyloseq <- otu_table(merged, taxa_are_rows = T)
rarecurve(t(species_phyloseq), step = 10000, label = T, ylab = "species")

# alpha diversity
plot_richness(species_phyloseq)

# get alpha diversity numbers
alpha <- estimate_richness(species_phyloseq)
alpha

# get Bray-Curtis beta-diversity numbers
beta <- distance(species_phyloseq, "bray")
beta
sort(beta)

# show order distance matrix
m <- as.matrix(beta)
m[upper.tri(m)] <- 0
library(reshape2)
mm <- subset(melt(m), value!=0)
mm[order(mm$value), 3:1]

# ordered samples -->

## Shannon diversity cmparison (PD and CONTROL)

library(ggpubr)
library(ggplot2)



#########
# Plots #
#########

# convert matrix into data.frame
data_PD <- data.frame(species_percentage_filter_PD)
data_control <- data.frame(species_percentage_filter_control)
data_general <- data.frame(species_percentage_filter_general)

# convert the first column into species name
setDT(data_PD, keep.rownames = TRUE)
colnames(data_PD)[1] <- "specie"

setDT(data_control, keep.rownames = TRUE)
colnames(data_control)[1] <- "specie"

setDT(data_general, keep.rownames = TRUE)
colnames(data_general)[1] <- "specie"
View(data_general)


## stacked barplot preparation
# reshape data into long format
data_long_PD <- melt(data_PD, id.vars = "specie", variable.name = "sample")
data_long_control <- melt(data_control, id.vars = "specie", variable.name = "sample")
data_long_general <- melt(data_general, id.vars = "specie", variable.name = "sample")
View(data_long_general)




## make the stacked barplot
library(forcats)
PD_plot <- ggplot(data_long_PD, aes(x = sample, y = value, fill = specie)) +
  geom_bar(stat = "identity", col = "black", position = "stack") +
  ylab("Abundance (%)") + coord_flip() + scale_fill_discrete(name = "Species")+theme(axis.title.y=element_blank())+ ggtitle("PD samples")

control_plot <- ggplot(data_long_control, aes(x = sample, y = value, fill = specie)) +
  geom_bar(stat = "identity", col = "black", position = "stack") +
  ylab("Abundance (%)") + coord_flip() + scale_fill_discrete(name = "Species")+theme(axis.title.y=element_blank())+ ggtitle("control samples")

general_plot <- ggplot(data_long_general, aes(x = sample, y = value, fill = specie)) +
  geom_bar(stat = "identity", col = "black", position = "stack") +
  ylab("Abundance (%)") + coord_flip() + scale_fill_discrete(name = "Species")+theme(axis.title.y=element_blank())

PD_plot
control_plot
general_plot

# make interactive barplot
library(plotly)
ggplotly(PD_plot)
ggplotly(control_plot)
  

# remove this if it gives problems
## add sample Status
library(dplyr)

# convert sample column factor into characters

dataFactor <- data_long_general
dataFactor <- data.frame(lapply(dataFactor, as.character), stringsAsFactors = FALSE)
# adding column status based on other column
dataFactor <- dataFactor %>%
  mutate(status = case_when(
    endsWith(sample, "3012") ~ "PD",
    endsWith(sample, "2941") ~ "PD",
    endsWith(sample, "2948") ~ "PD",
    endsWith(sample, "3009") ~ "PD",
    endsWith(sample, "3002") ~ "PD",
    endsWith(sample, "2986") ~ "PD",
    endsWith(sample, "2987") ~ "PD",
    endsWith(sample, "2998") ~ "PD",
    endsWith(sample, "3011") ~ "PD",
    endsWith(sample, "3000") ~ "PD",
    endsWith(sample, "3017") ~ "PD",
    endsWith(sample, "3005") ~ "PD",
    endsWith(sample, "2988") ~ "PD",
    endsWith(sample, "2994") ~ "PD",
    endsWith(sample, "2995") ~ "PD",
    endsWith(sample, "2997") ~ "PD",
    endsWith(sample, "2989") ~ "PD",
    endsWith(sample, "2990") ~ "PD",
    endsWith(sample, "3001") ~ "PD",
    endsWith(sample, "2992") ~ "PD",
    endsWith(sample, "2951") ~ "control",
    endsWith(sample, "2922") ~ "control",
    endsWith(sample, "2955") ~ "control",
    endsWith(sample, "2956") ~ "control",
    endsWith(sample, "2962") ~ "control",
    endsWith(sample, "2963") ~ "control",
    endsWith(sample, "2982") ~ "control",
    endsWith(sample, "3026") ~ "control",
    endsWith(sample, "3027") ~ "control",
    endsWith(sample, "2932") ~ "control",
    endsWith(sample, "2950") ~ "control",
    endsWith(sample, "2952") ~ "control",
    endsWith(sample, "2954") ~ "control",
    endsWith(sample, "2957") ~ "control",
    endsWith(sample, "2960") ~ "control",
    endsWith(sample, "2993") ~ "control",
    endsWith(sample, "2934") ~ "control",
    endsWith(sample, "2968") ~ "control",
    endsWith(sample, "2969") ~ "control",
    endsWith(sample, "3028") ~ "control"
  ))
View(dataFactor)  

# from character to factor

dataFinal1 <- dataFactor
dataFinal1$specie <- as.factor(dataFinal1$specie)
dataFinal1$sample <- as.factor(dataFinal1$sample)

# from character to numeric
dataFinal1$value <- as.numeric(dataFinal1$value)

# plot
library(ggplot2)
library(ggplotly)
level_order <- factor(dataFinal1$sample, level = c("SRR15853012", "SRR15852941", "SRR15852948","SRR15853009", "SRR15853002",
                                                   "SRR15852986", "SRR15852987", "SRR15852998", "SRR15853011","SRR15853000",
                                                   "SRR15853017","SRR15853005", "SRR15852988", "SRR15852994","SRR15852995", 
                                                   "SRR15852997", "SRR15852989","SRR15852990","SRR15853001","SRR15852992",
                                                   "SRR15852951", "SRR15852922", "SRR15852955", "SRR15852956","SRR15852962",
                                                   "SRR15852963", "SRR15852982","SRR15853026","SRR15853027","SRR15852932",
                                                   "SRR15852950", "SRR15852952","SRR15852954","SRR15852957","SRR15852960",
                                                   "SRR15852993", "SRR15852934","SRR15852968","SRR15852969","SRR15853028"))
status_plot <- ggplot(dataFinal1, aes(x = level_order, y = value, fill = specie)) +
  geom_bar(stat = "identity", col = "black", position = "stack") +
  ylab("Abundance (%)") + scale_fill_discrete(name = "Species")+
  facet_grid(.~status, scales="free", space = "free_x") + xlab("samples") + ylab("Abundance (%)")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  guides(color = guide_legend(override.aes = list(size = 5)))+
  theme(legend.title = element_text(size = 15), 
        legend.text = element_text(size = 11))
  


status_plot

ggplotly(status_plot)




## creation of the ML taxonomic matrix

############
# Sum the columns to verify that all numbers are 100
colSums(species_percentage_general[ , 1:40])

# transform matrix into a dataframe
taxo_data <- as.data.frame(species_percentage_general)

# convert to matrix
taxo_matrix <- as.matrix(taxo_data)

## metadata
# To aggregate the metadata and to be allow to use the machine learning algorithmns we need
# to have our samples as rows and the species name as columns
# so we have to transpose our matrix

# transpose matrix
taxo_matrix_t <- data.frame(t(taxo_matrix))



##

# metadata vector
condition <- c("PD","PD","PD","PD","PD","PD","PD","PD","PD","PD","PD","PD","PD","PD","PD",
               "PD","PD","PD","PD","PD","control","control","control","control","control",
               "control","control","control","control","control","control","control",
               "control","control","control","control","control","control","control","control")

# 1 is PD and 0 is control
condition2 <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
                0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)

# attach the metadata to the dataframe
taxonomic_matrix <- cbind(taxo_matrix_t,condition)
taxonomic_matrix2 <- cbind(taxo_matrix_t,condition2)


# Download data frame to CSV file
write.csv(taxonomic_matrix2, file = '/media/sequentia/visitors/visitor8/TFG/taxonomic_matrix2.csv')


