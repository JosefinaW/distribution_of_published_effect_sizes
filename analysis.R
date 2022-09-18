#Necessary packages for the analysis
install.packages("plotly")
install.packages("gridExtra")
install.packages("ggExtra")
install.packages("RColorBrewer")
install.packages("pwr")
install.packages("PEIP")
install.packages("plyr")
install.packages("dplyr")
install.packages("bayestestR")
install.packages("hexbin")
install.packages("ggExtra")
install.packages("spearmanCI")
install.packages("matlab")

library(plotly)
library(ggExtra)
library(plyr)
library(dplyr)
library(pwr)
library(PEIP)
library(tidyverse)
library(plyr)
library(gridExtra)
library(ggExtra)
library(RColorBrewer)
library(hexbin)
library(bayestestR)
library(spearmanCI)
library(matlab)



##FUNCTIONS
#Function to find p-value for r and n to check against the reported p-value
CorToP <- function(r,n){ 
  T_value <- r*sqrt((n-2)/(1-r^2))
  P_value <- 2*pt(T_value, n-2, lower=FALSE)
  return(P_value)}

#the same function but one sided
OSCorToP <- function(r,n){ 
  T_value <- r*sqrt((n-2)/(1-r^2))
  P_value <- pt(T_value, n-2, lower=FALSE)
  return(P_value)}



##SAMPLE1 data
#data frame preparation
#the data are loaded from two csv files
data.p2 <- read_csv("Data260121p2.csv")
data.p1 <- read_csv("Data260121p1.csv")

#changing column names for easier manipulation
colnames(data.p2)[4] <- ("SignificanceLevel")
colnames(data.p2)[3] <- ("Direction")
colnames(data.p2)[13] <- ("InPaper")
colnames(data.p2)[14] <- ("InTable")
colnames(data.p2)[15] <- ("type_p")
colnames(data.p2)[17] <- ("CorrType")
colnames(excl_non_spec)[10] <-("tpf")


#changing corr and signif columns into numeric values
data.p2$Correlation <- as.numeric(as.character(data.p2$Correlation))
data.p2$SignificanceLevel <- as.numeric(as.character(data.p2$SignificanceLevel))
#applying the fuction to get exact p-values given certain sample size and r for two-sided and one-sided variant
data.p2$Exact_P <- mapply(CorToP, r=data.p2$Correlation, n=data.p2$N)
data.p2$OSExact_P <- mapply(OSCorToP, r=data.p2$Correlation, n=data.p2$N)

##CREATING MERGED DF, merged by StudyID
df <- merge(data.p1, data.p2, by = "StudyID")

#replacing NAs with 0 in column specifying "reported as nonsignificant"
df$Nonsignificant[is.na(df$Nonsignificant)] <- 0

#adding a column specifying subfield to the df
targetdev <- c("ChildDev","JCPP","DevPsy")
targetsoc <- c("JPSP", "SPPS","EJP")
df$subfield <- ifelse(df$Journal %in% targetdev, "dev",
                      ifelse(df$Journal %in% targetsoc, "soc", "somethingelse"))

#r values by field for only records where significance is specified
excl_non_spec <- df  %>% filter((type_p > 0) | (Nonsignificant > 0))


#creating df containing only values reported as nonsignificant
nonsig_only <- excl_non_spec %>% filter(Nonsignificant ==1)

#creating df containing only non-specified values
non_spec <- df %>%
  filter(Nonsignificant == 0) %>%
  filter(type_p ==0)

#creating df containing only reportedly significant values
signific_only <- excl_non_spec %>% filter (Nonsignificant==0)

#misreporting of significant values - this was not part of the manuscript but was used to check the reason for over the significance boundary values in Figure 3 A and B.
P_larger_Sig <- df %>%
  filter(type_p != 3) %>% filter(Nonsignificant != 1) %>%
  filter(SignificanceLevel < OSExact_P) #or filter(SignificanceLevel < Exact_P)

#how many corrs reported as nonsignificant with reported p-value are significant at 0.05 level when exact p calculated
nonsig_p <- excl_non_spec %>%
  filter(Nonsignificant==1) %>% filter(SignificanceLevel>Exact_P) %>% filter (SignificanceLevel>0.05&Exact_P <= 0.05)

#how many correlations reported as nonsignificant without p-value reported are significant at 0.05
nonsig_p2 <- excl_non_spec%>%
  filter(Nonsignificant==1)%>%
  filter(is.na(SignificanceLevel)) %>% filter(OSExact_P < 0.05)  

#values reported with p=
peq <- df%>%filter(type_p==1)

#Manuscript figures
#percentile and cumulative percentile plot for Sample 1 (= Figure 2 A )

dfS1 <- df%>%select("Correlation","N")%>%filter(N>0)
dfS1 <- mutate(dfS1, df = N-2) #creating degrees of freedom column from N 
colnames(dfS1)<-c("rabs","N","df") #change column names
chdf3 <- dfS1 %>% arrange(df) #creating dataframe arranged by df
#loop to calculate cumulative percentiles so that for each df, the percentiles were calculated for all r-values with that df and smaller df.
qTable3 <- data.frame()

for (x in unique(chdf3$df)) {
  chdfRow <- chdf3[chdf3$df<=x,c("df","rabs")]
  qRow <- t(data.frame(quantile(chdfRow[["rabs"]], na.rm = TRUE)))
  qRow["df"] <- max(chdfRow["df"])
  qTable3 <- rbind(qTable3,qRow)
}
colnames(qTable3) <- c("q0","q1","q2","q3","q4","df") #change column names for quantiles and df
#calculating median values for Figure 2 A
chdf3 <- chdf3 %>% mutate(id = row_number()) 
chdf3$category <- cut(chdf3$id, seq(0,12406,1000), labels=c(1:12))
chdf_median <- chdf3%>%
  group_by(category)%>% 
  summarise(median=median(rabs), dfmedian=median(df))
chdf_median <- chdf_median%>%select(dfmedian, median)

#plot for the cumulative percentiles (Figure 2 A) for Sample 1
perplot1 <- ggplot(qTable3, aes(x=q2, y=N, color="q2")) + geom_point() + scale_y_log10() + geom_point(data=qTable3, aes(x=q1,y=N, color="q1")) + geom_point(data=qTable3, aes(x=q3,y=N, color="q3")) + theme_bw() + 
  labs(title="A. Cumulative 25th, 50th and 75th percentiles of correlations \n by df (Sample 1)", 
       x="Correlation",y="degrees of freedom (log10)") +
  theme(plot.title = element_text(size=18), plot.subtitle = element_text(size=14), 
        legend.text = element_text(size = 14), axis.title=element_text(size=14), legend.position = "none") + 
  scale_color_discrete(name = "", labels=c("25th percentile","50th percentile", "75th percentile")) + geom_point(data=chdf_median, aes(x=median, y=dfmedian), color="black", shape=17, size=3)


#Figure 4 A, B and C
#plotting density of sig and nonsig records at different sample sizes for Sample 1
#Figure 4A
d1<-ggplot(excl_non_spec, aes(x=N,y = (..count..)/sum(..count..), fill=factor(Nonsignificant))) + geom_density(alpha=0.5) + scale_x_log10() + labs(title= "A. Probability density distribution of Sample size for significant and nonsignificant correlation values", y="Probability", x = "Sample size (log10)") + scale_fill_discrete(name = "", labels=c("Significant","Nonsignificant")) +
  theme_bw() + theme(plot.title = element_text(size=18), plot.subtitle = element_text(size=14), legend.text = element_text(size = 14), axis.title=element_text(size=14)) 
#Figure 4B
d2<-ggplot(excl_non_spec, aes(x=Correlation,y = (..count..)/sum(..count..), color=factor(Nonsignificant))) + geom_density( size=1 ) + labs(title= "B. Probability density distribution of significant and nonsignificant correlation values", y="Probability", x = "Correlation values (r)") + scale_color_discrete(name = "", labels=c("Significant","Nonsignificant")) +
  theme_bw() + theme(plot.title = element_text(size=18), plot.subtitle = element_text(size=14), legend.text = element_text(size = 14), axis.title=element_text(size=14)) 
#Figure 4C
d3 <-ggplot(excl_non_spec, aes(x=Correlation, color=factor(Nonsignificant))) + stat_ecdf(size=1, geom="line") + labs(title= "C. Cumulative probability distribution of significant and nonsignificant correlation values", y="Cumulative probability", x = "Correlation values (r)") + scale_color_discrete(name = "", labels=c("Significant","Nonsignificant")) +
  theme_bw() + theme(plot.title = element_text(size=18), plot.subtitle = element_text(size=14), legend.text = element_text(size = 14), axis.title=element_text(size=14)) 
#arranging together
grid.arrange(d1,d2,d3)

##SAMPLE2 
#data are uploaded from this csv 
sample2 <- read_csv("CorrResults140820.csv")
sample2<- sample2%>%filter(abs(r)<=1) #ensuring that only realistic values remain in the sample
#creating column with different years
sample2$Year <- ifelse(endsWith(sample2$Dir, "10"),2010, 
                       ifelse(endsWith(sample2$Dir, "11"), 2011,
                              ifelse(endsWith(sample2$Dir, "12"),2012,
                                     ifelse(endsWith(sample2$Dir, "13"),2013,
                                            ifelse(endsWith(sample2$Dir, "14"),2014,
                                                   ifelse(endsWith(sample2$Dir, "15"),2015,
                                                          ifelse(endsWith(sample2$Dir, "16"),2016,
                                                                 ifelse(endsWith(sample2$Dir, "17"),2017,
                                                                        ifelse(endsWith(sample2$Dir, "18"),2018,
                                                                               ifelse(endsWith(sample2$Dir, "19"),2019, 0)))))))))) 



#creating a column specifying subfield for sample 2

targetdev2 <- c("\\ChildDev","\\DevPsy","\\JAppDevPsy","\\JCPP","\\JExpChildPs")
targetsoc2<- c("\\EJP","\\JExpSocPsy","\\JPSP", "\\JResInPer", "\\SPPS")

sample2$subfield <- ifelse(startsWith(sample2$Dir,"\\ChildDev"),"dev",
                           ifelse(startsWith(sample2$Dir,"\\DevPsy"),"dev",
                                  ifelse(startsWith(sample2$Dir,"\\JAppDevPsy"),"dev",
                                         ifelse(startsWith(sample2$Dir,"\\JCPP"),"dev",
                                                ifelse(startsWith(sample2$Dir,"\\JExpChildPs"),"dev", "soc")))))


#histogram of r-values for each year (Figure 5)
plot1<-ggplot(sample2, aes(x=abs(r), y=(..count..)/sum(..count..), fill="Year")) + geom_histogram(bins=100) + facet_wrap(~Year, nrow=2) + theme_bw() + 
  labs(title = "Distribution of r-values across years 2010-2019", x="Correlations", y="Probability", fill="") + 
  theme(legend.position = "none", plot.title = element_text(size=18), plot.subtitle = element_text(size=14), legend.text = element_text(size = 14), axis.title=element_text(size=14))

#calculating the spearman correlation coefficient between df
sample2df <- sample2 %>% filter(df>0)
s210 <- sample2df %>% filter(Year==2010)
s211 <- sample2df %>% filter(Year==2011)
s212 <- sample2df %>% filter(Year==2012)
s213 <- sample2df %>% filter(Year==2013)
s214 <- sample2df %>% filter(Year==2014)
s215 <- sample2df %>% filter(Year==2015)
s216 <- sample2df %>% filter(Year==2016)
s217 <- sample2df %>% filter(Year==2017)
s218 <- sample2df %>% filter(Year==2018)
s219 <- sample2df %>% filter(Year==2019)
cor.test(abs(s210$r), s210$df, method = "spearman")
cor.test(abs(s211$r), s211$df, method = "spearman")
cor.test(abs(s212$r), s212$df, method = "spearman")
cor.test(abs(s213$r), s213$df, method = "spearman")
cor.test(abs(s214$r), s214$df, method = "spearman")
cor.test(abs(s215$r), s215$df, method = "spearman")
cor.test(abs(s216$r), s216$df, method = "spearman")
cor.test(abs(s217$r), s217$df, method = "spearman")
cor.test(abs(s218$r), s218$df, method = "spearman")
cor.test(abs(s219$r), s219$df, method = "spearman")
#calculating CIs for the correlations
spearmanCI(abs(s210$r), s210$df, level=0.95)

#calculating exact p-values (two-sided and one-sided)
sample2df$Exact_P <- mapply(CorToP, r=abs(sample2df$r), n=(sample2df$df+2)) 
sample2df$OSExact_P <- mapply(OSCorToP, r=abs(sample2df$r), n=(sample2df$df+2))

#calculating bootstrapped 95% confidence intervals for the median correlation values in different years for Figure 6 A
install.packages("boot")
library(boot)
s210 <- sample2%>%filter(Year==2010)
Mboot = boot(abs(s210$r),
             function(x,i) median(x[i]), R=10000)
CI10 <- boot.ci(boot.out=Mboot, type="norm")
s211 <- sample2%>%filter(Year==2011)
Mboot = boot(abs(s211$r),
             function(x,i) median(x[i]), R=10000)
CI11 <- boot.ci(boot.out=Mboot, type="norm")
s212 <- sample2%>%filter(Year==2012)
Mboot = boot(abs(s212$r),
             function(x,i) median(x[i]), R=10000)
CI12 <- boot.ci(boot.out=Mboot, type="norm")
s213 <- sample2%>%filter(Year==2013)
Mboot = boot(abs(s213$r),
             function(x,i) median(x[i]), R=10000)
CI13 <- boot.ci(boot.out=Mboot, type="norm")
s214 <- sample2%>%filter(Year==2014)
Mboot = boot(abs(s214$r),
             function(x,i) median(x[i]), R=10000)
CI14 <- boot.ci(boot.out=Mboot, type="norm")
s215 <- sample2%>%filter(Year==2015)
Mboot = boot(abs(s215$r),
             function(x,i) median(x[i]), R=10000)
CI15 <- boot.ci(boot.out=Mboot, type="norm")
s216 <- sample2%>%filter(Year==2016)
Mboot = boot(abs(s216$r),
             function(x,i) median(x[i]), R=10000)
CI16 <- boot.ci(boot.out=Mboot, type="norm")
s217 <- sample2%>%filter(Year==2017)
Mboot = boot(abs(s217$r),
             function(x,i) median(x[i]), R=10000)
CI17 <- boot.ci(boot.out=Mboot, type="norm")
s218 <- sample2%>%filter(Year==2018)
Mboot = boot(abs(s218$r),
             function(x,i) median(x[i]), R=10000)
CI18 <- boot.ci(boot.out=Mboot, type="norm")
s219 <- sample2%>%filter(Year==2019)
Mboot = boot(abs(s219$r),
             function(x,i) median(x[i]), R=10000)
CI19 <- boot.ci(boot.out=Mboot, type="norm")
CI10 <- CI10$normal
CI11 <- CI11$normal
CI12 <- CI12$normal
CI13 <- CI13$normal
CI14 <- CI14$normal
CI15 <- CI15$normal
CI16 <- CI16$normal
CI17 <- CI17$normal
CI18 <- CI18$normal
CI19 <- CI19$normal
CIwhole<-rbind(CI10,CI11,CI12,CI13,CI14,CI15,CI16,CI17,CI18,CI19) #creating one data frame with all confidence intervals

med <- tapply(abs(sample2$r), sample2$Year, median) #calculating medians for each these
med <- data.frame(med)
CIwhole <- data.frame(CIwhole, med, Year) #joining into one data frame
colnames(CIwhole)<-c("Level","CIlow","CIhigh","median", "Year")
#creating Figure 6 A
CIplot <-ggplot(CIwhole, aes(x=Year, y=median, colour="Median")) + geom_line(size=1.5) + scale_x_continuous(breaks = c(2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019)) + geom_line(data=CIwhole, aes(x=Year, y=CIlow, colour="Confidence Interval"), size=1, alpha=0.5) + geom_line(data=CIwhole, aes(x=Year, y=CIhigh, colour="Confidence Interval"), size=1, alpha=0.5) +
  theme_bw() + theme(legend.position = "right",plot.title = element_text(size=18), legend.title = element_blank(), legend.text = element_text(size = 14), axis.title=element_text(size=14)) + labs(title="A. Median with bootstrapped 95% Confidence Intervals ", y="Correlation", legend.title=element_blank()) + expand_limits(y=0)

#percentile plot (Figure 6 B and C) by subfield
#dividing by subfield
dev2 <- sample2 %>% filter(subfield=="dev")
soc2 <- sample2 %>% filter(subfield=="soc")
#basis for the plot of quartiles by subfield across years
#quartiles by year for developmental
meddev <- tapply(abs(dev2$r),dev2$Year, median,na.rm=TRUE)
q1dev <- tapply(abs(dev2$r),dev2$Year,quantile,prob=0.25, na.rm=TRUE)
q2dev <- tapply(abs(dev2$r),dev2$Year,quantile,prob=0.75, na.rm=TRUE)
Year <- c(2010:2019)
devDD <- data.frame(Year,q1dev,meddev,q2dev) #joined into dataframe
#quartiles by year for social
medsoc <- tapply(abs(soc2$r),soc2$Year, median,na.rm=TRUE)
q1soc <- tapply(abs(soc2$r),soc2$Year,quantile,prob=0.25, na.rm=TRUE)
q2soc <- tapply(abs(soc2$r),soc2$Year,quantile,prob=0.75, na.rm=TRUE)
Year <- c(2010:2019)
socDD <- data.frame(Year,q1soc,medsoc,q2soc) #dataframe

#plot of quartiles of r by subfield across years (Figure 6 B)
plotS1 <- ggplot(devDD, aes(x=Year, y=meddev, colour="50th percentile",linetype="developmental")) + geom_line(size=1.5)  + scale_x_continuous(breaks = c(2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019)) + geom_line(data=devDD, aes(x=Year, y=q1dev, colour="25th percentile"), size=1.5) + geom_line(data=devDD, aes(x=Year, y=q2dev, colour="75th percentile"), size=1.5) +
  theme_bw() + theme(legend.position = "right",legend.box = "vertical",plot.title = element_text(size=18), legend.title = element_blank(), legend.text = element_text(size = 14), axis.title=element_text(size=14)) + labs(title="B. 25th, 50th and 75th percentiles of r-values for developmental and social psychology for years 2010-2019", y="Correlation", legend.title=element_blank()) +
  geom_line(data=socDD, aes(x=Year, y=q1soc, colour="25th percentile", linetype="social"), size=1.5) + geom_line(data=socDD, aes(x=Year, y=medsoc, colour="50th percentile", linetype="social"), size=1.5) + geom_line(data=socDD, aes(x=Year, y=q2soc, colour="75th percentile", linetype="social"), size=1.5) + scale_y_continuous(limits=c(0,NA))

#basis for plot of quartiles for degrees of freedom across years
dev3<-sample2df%>%filter(subfield=="dev")
soc3<-sample2df%>%filter(subfield=="soc")
#quartiles for developmental by year
meddev2 <- tapply(dev3$df,dev3$Year, median,na.rm=TRUE)
q1dev2 <- tapply(dev3$df,dev3$Year,quantile,prob=0.25, na.rm=TRUE)
q2dev2 <- tapply(dev3$df,dev3$Year,quantile,prob=0.75, na.rm=TRUE)
Year <- c(2010:2019)
dfdevDD<- data.frame(Year,q1dev2,meddev2,q2dev2)
#quartiles for social by year
medsoc2 <- tapply(soc3$df,soc3$Year, median,na.rm=TRUE)
q1soc2 <- tapply(soc3$df,soc3$Year,quantile,prob=0.25, na.rm=TRUE)
q2soc2 <- tapply(soc3$df,soc3$Year,quantile,prob=0.75, na.rm=TRUE)
Year <- c(2010:2019)
dfsocDD <- data.frame(Year,q1soc2,medsoc2,q2soc2)

#plots for degrees of freedom quartiles across years (Figure 6 C)
plotS2 <-ggplot(dfdevDD, aes(x=Year, y=meddev2, colour="50th percentile", linetype="developmental")) + geom_line(size=1.5) + scale_x_continuous(breaks = c(2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019)) + geom_line(data=dfdevDD, aes(x=Year, y=q1dev2, colour="25th percentile"), size=1.5) + geom_line(data=dfdevDD, aes(x=Year, y=q2dev2, colour="75th percentile"), size=1.5) +
  theme_bw() + scale_y_log10() + theme(legend.position = "right",legend.box = "vertical",plot.title = element_text(size=18), legend.title = element_blank(), legend.text = element_text(size = 14), axis.title=element_text(size=14)) + labs(title="C. 25th, 50th and 75th percentiles of df values for developmental and social psychology for years 2010-2019", y="degrees of freedom (log10)", legend.title=element_blank()) +
  geom_line(data=dfsocDD, aes(x=Year, y=q1soc2, colour="25th percentile", linetype="social"), size=1.5) + geom_line(data=dfsocDD, aes(x=Year, y=medsoc2, colour="50th percentile", linetype="social"), size=1.5) + geom_line(data=dfsocDD, aes(x=Year, y=q2soc2, colour="75th percentile", linetype="social"), size=1.5)

grid.arrange(CIplot, plotS1, plotS2) #arranging together all parts of Figure 6 A, B, C
             
#the Figure 2 B for Sample 2
#Sample2 percentile and cumulative percentile plot
dfS2 <- sample2df%>%select("r","df") #from sample 2 data frame with only records containing degrees of freedom selecting r column and df column
dfS2 <- mutate(dfS2, rabs = abs(r)) #creating column of absolute r value 
chdf3 <- dfS2 %>% arrange(df) #creating dataframe arranged by df
#loop to calculate cumulative percentiles so that for each df, the percentiles were calculated for all r-values with that df and smaller df.
qTable3 <- data.frame()

for (x in unique(chdf3$df)) {
  chdfRow <- chdf3[chdf3$df<=x,c("df","rabs")]
  qRow <- t(data.frame(quantile(chdfRow[["rabs"]], na.rm = TRUE)))
  qRow["df"] <- max(chdfRow["df"])
  qTable3 <- rbind(qTable3,qRow)
}
colnames(qTable3) <- c("q0","q1","q2","q3","q4","N")

chdf3 <- chdf3 %>% mutate(id = row_number()) #row numbers as id of the rows for use in the below operation
#calculating the medians for each 300 values going from the lowest degrees of freedom
chdf3$category <- cut(chdf3$id, seq(0,3292,300), labels=c(1:10))
View(chdf3)
chdf_median <- chdf3%>%
  group_by(category)%>% 
  summarise(median=median(rabs), dfmedian=median(df))
chdf_median <- chdf_median%>%select(dfmedian, median) #data frame with medians for every 300 values

#plots for the cumulative percentiles: first one is a scatter plot showing the exact values
perplot3 <- ggplot(qTable3, aes(x=q2, y=N, color="q2")) + geom_point() + scale_y_log10() + geom_point(data=qTable3, aes(x=q1,y=N, color="q1")) + geom_point(data=qTable3, aes(x=q3,y=N, color="q3")) + theme_bw() + 
  labs(title="B. Cumulative 25th, 50th and 75th percentiles of correlations \n by df (Sample 2)", 
       x="Correlation",y="degrees of freedom (log10)") +
  theme(plot.title = element_text(size=18), plot.subtitle = element_text(size=14), 
        legend.text = element_text(size = 14), axis.title=element_text(size=14), legend.position = "none") + 
  scale_color_discrete(name = "", labels=c("25th percentile","50th percentile", "75th percentile")) + geom_point(data=chdf_median, aes(x=median, y=dfmedian), color="black", shape=17, size=3)

#arranging both plots together
grid.arrange(perplot1, perplot3, ncol=2) #arranging together Figure 2 A and B
             

#data frame of sample 2 records without studies in sample 1
sample_elim <- read_csv("sample2_elim_17092021.csv")
sample_elim<- sample_elim%>%filter(abs(r)<=1) #ensuring that only realistic values remain in the sample
#adding a column with year
sample_elim$Year <- ifelse(endsWith(sample_elim$Dir, "10"),2010, 
                           ifelse(endsWith(sample_elim$Dir, "11"), 2011,
                                  ifelse(endsWith(sample_elim$Dir, "12"),2012,
                                         ifelse(endsWith(sample_elim$Dir, "13"),2013,
                                                ifelse(endsWith(sample_elim$Dir, "14"),2014,
                                                       ifelse(endsWith(sample_elim$Dir, "15"),2015,
                                                              ifelse(endsWith(sample_elim$Dir, "16"),2016,
                                                                     ifelse(endsWith(sample_elim$Dir, "17"),2017,
                                                                            ifelse(endsWith(sample_elim$Dir, "18"),2018,
                                                                                   ifelse(endsWith(sample_elim$Dir, "19"),2019, 0)))))))))) #want to specify those that are nonsig and at the same time have p < 0.05



#creating a column specifying subfield for sample 2
targetdev2 <- c("\\ChildDev","\\DevPsy","\\JAppDevPsy","\\JCPP","\\JExpChildPs")
targetsoc2<- c("\\EJP","\\JExpSocPsy","\\JPSP", "\\JResInPer", "\\SPPS")
sample_elim$subfield <- ifelse(startsWith(sample_elim$Dir,"\\ChildDev"),"dev",
                               ifelse(startsWith(sample_elim$Dir,"\\DevPsy"),"dev",
                                      ifelse(startsWith(sample_elim$Dir,"\\JAppDevPsy"),"dev",
                                             ifelse(startsWith(sample_elim$Dir,"\\JCPP"),"dev",
                                                    ifelse(startsWith(sample_elim$Dir,"\\JExpChildPs"),"dev", "soc")))))



             
#OVERLAP BETWEEN SAMPLES
#creating a df of the overlap in sample 2
s2overlap <- anti_join(sample2, sample_elim, by=c("X1"))
s1text <- df%>%filter(InTable==0)
colnames(s1text)[25]<-"r"
nrow(s1text)
length(unique(s1text$StudyID))
nrow(unique(s1text[c("Journal","FirstAuthor","PageN")]))
#deleting those studies in Sample1 not detected in Sample2
pages<-c(331,267,555,9,172,96,162,69,211,135,528,391,379,85,291,771,311,302,361,111,182,160,274,246,222,280,206)
nrow(s1text %>% filter(PageN != pages))
s1text_filt<-filter(s1text,(PageN %in% pages == FALSE))
s1text_filt <-filter(s1text_filt,(FirstAuthor != "Borkenau"))
#density plot showing the values from overlapping studies for both samples (Figure 1)
ggplot(s1text_filt, aes(x=r, col="Sample 1"))+ geom_density()+ geom_density(data=s2overlap,aes(x=abs(r), col="Sample 2")) +
  labs(title="Density distribution of the values from studies detected within both samples", 
       x="Correlation") +
  theme(plot.title = element_text(size=18), 
        legend.text = element_text(size = 14), axis.title=element_text(size=14)) + theme_bw()
             
             
##VALIDITY CHECK FOR SAMPLE 2
#sampling 20 random papers for validity check
val_check <- sample2[sample(nrow(sample2), 20), ]
val_check <- val_check$File
VC <- sample2 %>% filter(File %in% val_check)
View(VC)






