# make and evaluate protein abundance predictions based on models M1, M2, M3 and M4
######################################### the main configurable parameters in section below ###############################
# input and output file names (the input ones should be consistent with the available data files)
#infile_rna <- "nci60_rna.xlsx"               # rna data file
infile_rna <- "nci60_rna.tsv"                   # rna data file
#infile_proteome <- "nci60_prot.xlsx"          # proteome data file
infile_proteome <- "nci60_prot.tsv"             # proteome data file
is_name_col <- 1                                # set the value to 1 if the first column contains protein names, set to 0 therwise
# outfile_correlations <- "out_correlations.xlsx"   # output file with performance assessments
outfile_correlations_tsv <- "out_correlations.tsv"  # output file with performance assessments
# outfile_imputed <- "out_imputed.xlsx"         # output file with imputed data replacing 0-es and NA
outfile_imputed_tsv <- "out_imputed.tsv"        # output file with imputed data replacing 0-es and NA
# select model used for output file of imputed data 
# models: M1 - all rna, M2 all prot (execpt itself), M3 - own RNA and all prot (except itself), M4 - all RNA and all prot (except itself)
# valid values: 11,21,22,23,31,32,33,41,42,43 (the first digit = number of model, the second digit = number of iterations) 
imp_mod <- 42;
######################################### the main configurable parameters in section above ###############################
######################################### function definitions in section below #W#########################################
# construct T/F vector with T positions for rows were both p and r are non-zero
nonzero_rows <- function(p,r)
{
s1 <- apply(p,1, function(row) all (row !=0))
s2 <- apply(r,1, function(row) all (row !=0))
s <- s1 & s2
return(s)
}
#replace 0-os with averages per rows
row_awg <- function(df)
{
df1 <- df
rc <- nrow(df)
nz_count <- rowSums(df != 0)
for (i in 1:rc) if (nz_count[i] > 0) 
	{
	sum <- sum(df[i,])
	awg <- sum/nz_count[i]
	dc <- df[i,]
	dc[dc == 0] <- awg
	df1[i,] <- dc
	}
return(df1)
}
#replace 0-os with averages per columns
col_awg <- function(df)
{
df1 <- df
cc <- ncol(df)
nz_count <- colSums(df != 0)
for (i in 1:cc) if (nz_count[[i]] > 0) 
	{
	sum <- sum(df[,i])
	awg <- sum/nz_count[[i]]
	dc <- df[,i]
	dc[dc == 0] <- awg
	df1[,i] <- dc
	}
return(df1)
}
# construct M1 formula
m1_model <- function(i,cn,name_str)
{
f1 <- paste0(name_str,i:i,"] ~ ")
f2 <- paste0("r0[,",1:cn,"]")
f3 <- paste0(f2,collapse="+")
f <- paste0(paste(f1),paste(f3))
}
# construct M2 formula
m2_model <- function(i,cn,name_str)
{
f1 <- paste0(name_str,i:i,"] ~ ")
if (i == 1) f2 <- paste0("p_work[,",2:cn,"]")
if (i == cn) f2 <- paste0("p_work[,",1:(cn-1),"]")
if (i > 1 && i < cn) 
	{
	f2a <-paste0("p_work[,",1:(i-1),"]")
	f2b <-paste0("p_work[,",(i+1):cn,"]")
	f2 <- append(f2a,f2b)
	}
f3 <- paste0(f2,collapse="+")
f <- paste0(paste(f1),paste(f3))	
}
# construct M3 formula
m3_model <- function(i,cn,name_str)
{
f1 <- paste0(name_str,i:i,"] ~ ")
f2 <- paste0("p_work[,",1:cn,"]")
f3 <- paste0(f2,collapse="+")
f <- paste0(paste(f1),paste(f3))
}
# construct M4 formula
m4_model <- function(i,cn,name_str)
{
f1 <- paste0(name_str,i:i,"] ~ ")
if (i == 1) f2c <- paste0("p_work[,",2:cn,"]")
if (i == cn) f2c <- paste0("p_work[,",1:(cn-1),"]")
if (i > 1 && i < cn) 
	{
	f2a <-paste0("p_work[,",1:(i-1),"]")
	f2b <-paste0("p_work[,",(i+1):cn,"]")
	f2c <- append(f2a,f2b)
	}
f2d <- paste0("r0[,",1:cn,"]")
f2 <- append(f2c,f2d)
f3 <- paste0(f2,collapse="+")
f <- paste0(paste(f1),paste(f3))	
}
######################################### function definitions in section below ###########################################
######################################### executable code follows #########################################################
# read input files (NA are treated as 0-es for original data sets)
library(readxl)
# proteomics file
# p_in1 <- read_excel(infile_proteome)
p_in1 <- read.table(infile_proteome, sep = '\t', header = TRUE)
p_cc <- ncol(p_in1)
if (is_name_col == 1) {p_in2 <- p_in1[,2:p_cc]} else {p_in2 <- p_in1}
p_in <- suppressWarnings(as.data.frame(sapply(p_in2, as.numeric))) # convert to numeric just in case...
p_in[is.na(p_in)] <- 0 # replace all NA with 0
p0 <- log(p_in+1)
p_out <- p_in1 # the file for imputed data (all 0-es will be replaced with the predicted values)
# transcriptomics file
# r_in1 <- read_excel(infile_rna)
r_in1 <- read.table(infile_rna, sep = '\t', header = TRUE)
r_cc <- ncol(r_in1)
if (is_name_col == 1) {r_in2 <- r_in1[,2:r_cc]} else {r_in2 <- r_in1}
r_in <- suppressWarnings(as.data.frame(sapply(r_in2, as.numeric))) # convert to numeric just in case...
r_in[is.na(r_in)] <- 0 # replace all NA with 0
r0 <- log(r_in+1)
# data sets with row average substitions instead of proteome 0-es
p0r <- row_awg(p0)
# data sets with col average substitions instead of proteome 0-es
p0c <- col_awg(p0)
#
# create data frame for results - reuse column names of proteome file, number of rows should be incresed as needed 
out_cc <- ncol(p0)
out_rc <- 23
out <- data.frame(matrix(nrow = out_rc, ncol = out_cc))
colnames(out) = colnames(p0)
for (i in 1:out_rc) for (j in 1:out_cc) out[i,j] <- " "
# 
# create placeholder for 1st_2nd iteration predictions (not needed for m1 where more iterations does not make sense)
p1a_m2 <- p0 # all predicted values
p1_m2 <- p0 # only 0-es replaced by predictions
p2a_m2 <- p0 # all predicted values
p2_m2 <- p0 # only 0-es replaced by predictions
p1a_m3 <- p0 # all predicted values
p1_m3 <- p0 # only 0-es replaced by predictions
p2a_m3 <- p0 # all predicted values
p2_m3 <- p0 # only 0-es replaced by predictions
p1a_m4 <- p0 # all predicted values
p1_m4 <- p0 # only 0-es replaced by predictions
p2a_m4 <- p0 # all predicted values
p2_m4 <- p0 # only 0-es replaced by predictions
# build M3 linear models
d_numb <- 4 # number of digits to be printed
p0_cn <- ncol(p0)
p0_rn <- nrow(p0)
for (i in 1:p0_cn) {
# try to duplicate proteome matrix and replace own proteome column with rna column
p_work <- p0
p_work[,i] <- r0[,i]
# p - rna correlation + also correlations with averages instead of 0-es
c1 <- cor(p0[,i],r0[,i])
out[2,i] <- round(c1,digits=d_numb)
print(c1)
c1r <- cor(p0r[,i],r0[,i])
out[3,i] <- round(c1r,digits=d_numb)
print(c1r)
c1c <- cor(p0c[,i],r0[,i])
out[4,i] <- round(c1c,digits=d_numb)
print(c1c)
print(" ")
# p -  predicted p correlation (m1)
f <- m1_model(i,p0_cn,"p0[,")
lm <- lm(as.formula(f))
fc <- fitted(lm)
c2 <- cor(p0[,i],fc)
out[6,i] <- round(c2,digits=d_numb)
print(c2)
for (j in 1:p0_rn) if (p0[j,i] == 0) if (imp_mod == 11) p_out[j,i+is_name_col] <- fc[j]
# p -  predicted p correlation (m2)
f <- m2_model(i,p0_cn,"p0[,")
lm <- lm(as.formula(f))
fc <- fitted(lm)
c2 <- cor(p0[,i],fc)
out[7,i] <- round(c2,digits=d_numb)
out[13,i] <- round(c2,digits=d_numb)
print(c2)
# apply these predictions for 0 values in p1
for (j in 1:p0_rn) if (p1_m2[j,i] == 0) 
	{
	p1_m2[j,i] <- fc[j] 
	if (imp_mod == 21) p_out[j,i+is_name_col] <- fc[j]
	}
# p -  predicted p correlation (m3)
f <- m3_model(i,p0_cn,"p0[,")
lm <- lm(as.formula(f))
fc <- fitted(lm)
c2 <- cor(p0[,i],fc)
out[8,i] <- round(c2,digits=d_numb)
out[17,i] <- round(c2,digits=d_numb)
print(c2)
# apply these predictions for 0 values in p1
for (j in 1:p0_rn) if (p1_m3[j,i] == 0) 
	{
	p1_m3[j,i] <- fc[j] 
	if (imp_mod == 31) p_out[j,i+is_name_col] <- fc[j]
	}
# p -  predicted p correlation (m4)
f <- m4_model(i,p0_cn,"p0[,")
lm <- lm(as.formula(f))
fc <- fitted(lm)
c2 <- cor(p0[,i],fc)
out[9,i] <- round(c2,digits=d_numb)
out[21,i] <- round(c2,digits=d_numb)
print(c2)
# apply these predictions for 0 values in p1
fc <- fitted(lm)
for (j in 1:p0_rn)
	{
	if (p1_m4[j,i] == 0) p1_m4[j,i] <- fc[j]
	}
for (j in 1:p0_rn) if (p1_m4[j,i] == 0) 
	{
	p1_m4[j,i] <- fc[j] 
	if (imp_mod == 41) p_out[j,i+is_name_col] <- fc[j]
	}
# now repeat with row averages instead of 0-es
p_work <- p0r
p_work[,i] <- r0[,i]
f <- m3_model(i,p0_cn,"p0r[,")
# p -  predicted p correlation
lm <- lm(as.formula(f))
fc <- fitted(lm)
c3 <- cor(p0[,i],fc)
out[10,i] <- round(c3,digits=d_numb)
print(c3)
# now repeat with col averages instead of 0-es
p_work <- p0c
p_work[,i] <- r0[,i]
f <- m3_model(i,p0_cn,"p0c[,")
# p -  predicted p correlation
lm <- lm(as.formula(f))
fc <- fitted(lm)
c4 <- cor(p0[,i],fc)
out[11,i] <- round(c4,digits=d_numb)
print(c4)
# print column name
print(colnames(p0)[i])
print(" ")
}
#
#
# the 2nd iteration - apply M3 to p1
p1_cn <- ncol(p1_m3)
p1_rn <- nrow(p1_m3)
for (i in 1:p1_cn) {
# try to duplicate proteome matrix and replace own proteome column with rna column
# model m2
p_work <- p1_m2
p_work[,i] <- r0[,i]
f <-m2_model(i,p1_cn,"p1_m2[,")
# print(f)
lm <- lm(as.formula(f))
# p - rna correlation + also correlations with averages instead of 0-es
fc <- fitted(lm)
c1 <- cor(p0[,i],fc)
out[14,i] <- round(c1,digits=d_numb)
print(c1)
# apply these predictions for 0 values in p2
for (j in 1:p1_rn) if (p2_m2[j,i] == 0) 
	{
	p2_m2[j,i] <- fc[j] 
	if (imp_mod == 22) p_out[j,i+is_name_col] <- fc[j]
	}
# model m3
p_work <- p1_m3
p_work[,i] <- r0[,i]
f <-m3_model(i,p1_cn,"p1_m3[,")
# print(f)
lm <- lm(as.formula(f))
# p - rna correlation + also correlations with averages instead of 0-es
fc <- fitted(lm)
c1 <- cor(p0[,i],fc)
out[18,i] <- round(c1,digits=d_numb)
print(c1)
# apply these predictions for 0 values in p2
for (j in 1:p1_rn) if (p2_m3[j,i] == 0) 
	{
	p2_m3[j,i] <- fc[j] 
	if (imp_mod == 32) p_out[j,i+is_name_col] <- fc[j]
	}
# model m4
p_work <- p1_m4
p_work[,i] <- r0[,i]
f <-m4_model(i,p1_cn,"p1_m4[,")
# print(f)
lm <- lm(as.formula(f))
# p - rna correlation + also correlations with averages instead of 0-es
fc <- fitted(lm)
c1 <- cor(p0[,i],fc)
out[22,i] <- round(c1,digits=d_numb)
print(" ")
print(c1)
# apply these predictions for 0 values in p2
for (j in 1:p1_rn) if (p2_m4[j,i] == 0) 
	{
	p2_m4[j,i] <- fc[j] 
	if (imp_mod == 42) p_out[j,i+is_name_col] <- fc[j]
	}
}
#
#
# the 3rd iteration - apply M3 to p2
p2_cn <- ncol(p2_m3)
p2_rn <- nrow(p2_m3)
for (i in 1:p2_cn) {
# try to duplicate proteome matrix and replace own proteome column with rna column
# model m2
p_work <- p2_m2
p_work[,i] <- r0[,i]
f <-m2_model(i,p2_cn,"p2_m2[,")
# print(f)
lm <- lm(as.formula(f))
# p - rna correlation + also correlations with averages instead of 0-es
fc <- fitted(lm)
c2 <- cor(p0[,i],fc)
out[15,i] <- round(c2,digits=d_numb)
print(c2)
for (j in 1:p2_rn) if (p0[j,i] == 0) if (imp_mod == 23) p_out[j,i+is_name_col] <- fc[j]
# model m3
p_work <- p2_m3
p_work[,i] <- r0[,i]
f <-m3_model(i,p2_cn,"p2_m3[,")
# print(f)
lm <- lm(as.formula(f))
# p - rna correlation + also correlations with averages instead of 0-es
fc <- fitted(lm)
c2 <- cor(p0[,i],fc)
out[19,i] <- round(c2,digits=d_numb)
print(c2)
for (j in 1:p2_rn) if (p0[j,i] == 0) if (imp_mod == 33) p_out[j,i+is_name_col] <- fc[j]
# model m4
p_work <- p2_m4
p_work[,i] <- r0[,i]
f <-m4_model(i,p2_cn,"p2_m4[,")
# print(f)
lm <- lm(as.formula(f))
# p - rna correlation + also correlations with averages instead of 0-es
fc <- fitted(lm)
c2 <- cor(p0[,i],fc)
out[23,i] <- round(c2,digits=d_numb)
print(c2)
print(" ")
for (j in 1:p2_rn) if (p0[j,i] == 0) if (imp_mod == 43) p_out[j,i+is_name_col] <- fc[j]
}
# write output files
library(xlsx)
# write correlation file - add extra column with explanation what is shown within each row
out2 <- data.frame(matrix(nrow = out_rc, ncol = out_cc+1))
colnames(out2) <- append(" ",colnames(out))
out2[,2:(out_cc+1)] <- out[,1:out_cc]
for (i in 1:out_rc) out2[1,i] <- " "
out2[2,1] <- "prot-rna correlations, NA = 0" 
out2[3,1] <- "prot-rna correlations, NA = row average"
out2[4,1] <- "prot-rna correlations, NA = col average"
out2[5,1] <- " "
out2[6,1] <- "prot-M1 correlations"
out2[7,1] <- "prot-1xM2 correlations, NA = 0"
out2[8,1] <- "prot-1xM3 correlations, NA = 0" 
out2[9,1] <- "prot-1xM4 correlations, NA = 0" 
out2[10,1] <- "prot-1xM3 correlations, NA = row average"
out2[11,1] <- "prot-1xM3 correlations, NA = col average"
out2[12,1] <- " "
out2[13,1] <- "prot-1xM2 correlations, NA = 0" 
out2[14,1] <- "prot-2xM2 correlations, NA = 0" 
out2[15,1] <- "prot-3xM2 correlations, NA = 0" 
out2[16,1] <- " "
out2[17,1] <- "prot-1xM3 correlations, NA = 0" 
out2[18,1] <- "prot-2xM3 correlations, NA = 0" 
out2[19,1] <- "prot-3xM3 correlations, NA = 0" 
out2[20,1] <- " "
out2[21,1] <- "prot-1xM4 correlations, NA = 0" 
out2[22,1] <- "prot-2xM4 correlations, NA = 0" 
out2[23,1] <- "prot-3xM4 correlations, NA = 0" 
# for imputation output file replace negative values with 0-es 
c_out <- ncol(p_out)
for (i in (1+is_name_col):c_out) 
	{
	p_outc <- p_out[,i]
	p_outc[p_outc < 0] <- 0
	p_out[,i] <- p_outc
	}
# write files themselves
write.table(out2, outfile_correlations_tsv, sep ="\t", row.names=F)
write.table(p_out, outfile_imputed_tsv, sep ="\t", row.names=F)
# write.xlsx(out2,outfile_correlations)
# write.xlsx(p_out,outfile_imputed)
print("Done!")