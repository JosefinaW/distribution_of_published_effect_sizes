# distribution_of_published_effect_sizes
Data and analysis for the paper "Published effect sizes in Social and Developmental Psychology"
Brief summary of dataset contents, contextualized in experimental procedures and results.
Datasets contain manually collected correlation effect sizes (Data260121p1.csv and Data260121p2.csv) and computer extracted correlation effect sizes (CorrResults140820.csv and sample2_elim_17092021.csv) from the fields of social and developmental psychology.  


## Description of the Data and file structure
Csv files Data260121p1.csv and Data260121p2.csv can be merged together using the StudyID variables.   
In Data260121p1.csv each row refers to one study and the files includes following columns:  
StudyID - unique identifier of each study  
OwnData -study used external dataset (0) or own data (1)  
Journal - journal acronym  
Year - year of publication  
Volume - journal volume  
PageN - first page  
FirstAuthor - first author  
topic - topic of study  
MultipleStudies - study was part of a paper including multiple studies  
tpf (total FINAL participants) -total number of participants as given in methods section  
Groups - number of participant groups  
General_notes -any comments re. previous info  
AlternativeHypothesis - at least one hypothesis given as directional (directional), no hypotheses given as directional (explorative)  
NullHypothesis (mentioned 1, not mentioned 0) - null hypothesis mentioned  
CommentsHypotheses - comments  
TypesCorrelation - type of correlation if mentioned  
CommentsCorrelation - comments  
apri_powercalculation - apriori power calculation included?  
power_test - power test used  
Effect_sizes for power calculations - effect size specified in power calc  
Power_comments - comments  
Preregistration - preregistration link within study  
Replication_exp_in_paper - study was replication  
comments_from_paper - comments of note from within the study  

Data260121p2.csv is organised by single correlation per row with associated statistical and other information.   
It contains following rows.  
StudyID - study identifier  
Correlation - correlation as absolute value  
Directionality - positive or negative  
SignificanceLevel<than - level of significance or critical value given  
Nonsignificant - 1 if reported as nonsignificant  
N - sample size associated with given correlation  
Variable1 -construct measured  
Measure1 - measure  
Variable2 - construct 2 measured  
Measure2 - measure 2  
Same_variable - 1 if variable 1 and 2 are the same  
Type_Same - reason given for looking at same variable  
InPaper (as opposed to supplementary material) - correlation detected in study results or supplementary material  
InTable (1 in table, 0 in text) - 1 in table, 0 in text  
type_p(1: p=, 2: p<, 3: p>) - type of p-value given  
MainCor - attempt to collect data on main effect sizes of interest, later left incomplete due to too many uncertainties  
corr_Type(when different) - type of correlation if different from the main type used within study  
ProvisionallySuspicious - correlation seems strange  
Corr_Comments - any further comments  
N as range - sample size given as range  

CorrResults140820 are data collected from 8 journals from issues between the years 2010 and 2019. It contains following variables:  
Dir - path to original directory, given as journal acronym and year  
File - path to source  
Before - text before the target symbol combination  
Target - target value  
After - text right after the target symbol combination  
FullText - Target + After  
r1 - result of one approach to detect r-value, if approach does not work, gives impossible value  
r2 - result of second approach to detect r-value, if approach does not work, gives impossible value  
df - detected degrees of freedom  
pe - detected p = value  
ps- detected p<value  
pl - detected p>value  
r - r detected correctly by either r1 or r2  

sample2_elim_17092021.csv is a subsample of CorrResults140820 and shares the same column names.  


## Sharing/access Information  

Links to other publicly accessible locations of the data:  
https://github.com/JosefinaW/distribution_of_published_effect_sizes.git  

