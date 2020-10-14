library(ggplot2)
library(gridExtra)
library(cowplot)
library(fitdistrplus)
library(readr)
library(tibble)
library(dplyr)
library(ANOVAreplication)


#read in data and split the probands from sibs
setwd("/Users/jon/Research/scripts/de_novo_sv/sfari_denovos/r_plotting/data_files/")
sv_info=read.csv("dnsv_dataframe.csv", header=TRUE, sep=",")
sv_info$status = recode(sv_info$role, "sib"="Unaffected","proband"="Probands", "mother"="Unaffected","father"="Unaffected")

sv_info = sv_info[(sv_info$status == "Probands") | (sv_info$status == "Unaffected") ,]
sv_info = sv_info[sv_info$father_age_birth_years >0,]
sv_info$father_age_birth_years = as.integer(sv_info$father_age_birth_years)
sv_info$mother_age_birth_years = as.integer(sv_info$mother_age_birth_years)
sv_info$has_dnSV = factor(sv_info$all_sv > 0)
sv_info$dnSV = factor(sv_info$all_sv>0)
dad_sv_info = sv_info[(sv_info$paternal == 1) | (sv_info$dnSV == FALSE) ,]
mom_sv_info = sv_info[(sv_info$maternal == 1) | (sv_info$dnSV == FALSE) ,]



#Mann-Whitney-Wilcoxon for count data.
#########################correllation with dnSNVs
dad_sib_snv_test_stat_snvs = wilcox.test(alternative="greater",
  dad_sv_info$dn_snvs[(dad_sv_info$dnSV == TRUE) & (dad_sv_info$status == "Unaffected") ], 
  dad_sv_info$dn_snvs[(dad_sv_info$dnSV == FALSE) & (dad_sv_info$status == "Unaffected")]
)
dad_proband_snv_test_stat_snvs = wilcox.test(alternative="greater",
  dad_sv_info$dn_snvs[(dad_sv_info$dnSV == TRUE) & (dad_sv_info$status == "Probands") ], 
  dad_sv_info$dn_snvs[(dad_sv_info$dnSV == FALSE) & (dad_sv_info$status == "Probands")]
)

mom_sib_snv_test_stat_snvs = wilcox.test(alternative="greater",
  mom_sv_info$dn_snvs[(dad_sv_info$dnSV == TRUE) & (dad_sv_info$status == "Unaffected") ], 
  mom_sv_info$dn_snvs[(dad_sv_info$dnSV == FALSE) & (dad_sv_info$status == "Unaffected")]
)
mom_proband_snv_test_stat_snvs = wilcox.test(alternative="greater",
  mom_sv_info$dn_snvs[(mom_sv_info$dnSV == TRUE) & (mom_sv_info$status == "Probands") ], 
  mom_sv_info$dn_snvs[(mom_sv_info$dnSV == FALSE) & (mom_sv_info$status == "Probands")]
)

#########################correllation with parental age
dad_sib_test_stat_age = wilcox.test(alternative="greater",
  dad_sv_info$father_age_birth_years[(dad_sv_info$dnSV == TRUE) & (dad_sv_info$status == "Unaffected") ], 
  dad_sv_info$father_age_birth_years[(dad_sv_info$dnSV == FALSE) & (dad_sv_info$status == "Unaffected")]
)
dad_proband_test_stat_age = wilcox.test(alternative="greater",
  dad_sv_info$father_age_birth_years[(dad_sv_info$dnSV == TRUE) & (dad_sv_info$status == "Probands") ], 
  dad_sv_info$father_age_birth_years[(dad_sv_info$dnSV == FALSE) & (dad_sv_info$status == "Probands")]
)
mom_sib_test_stat_age = wilcox.test(alternative="greater",
  mom_sv_info$father_age_birth_years[(mom_sv_info$dnSV == TRUE) & (mom_sv_info$status == "Unaffected") ], 
  mom_sv_info$father_age_birth_years[(mom_sv_info$dnSV == FALSE) & (mom_sv_info$status == "Unaffected")]
)
mom_proband_test_stat_age = wilcox.test(alternative="greater",
  mom_sv_info$father_age_birth_years[(mom_sv_info$dnSV == TRUE) & (mom_sv_info$status == "Probands") ], 
  mom_sv_info$father_age_birth_years[(mom_sv_info$dnSV == FALSE) & (mom_sv_info$status == "Probands")]
)


##################dad plot
dad_lengths=c(
  length(dad_sv_info$dnSV[(dad_sv_info$status == "Unaffected") & (dad_sv_info$dnSV == FALSE)]),
  length(dad_sv_info$dnSV[(dad_sv_info$status == "Unaffected") & (dad_sv_info$dnSV == TRUE)]),
  length(dad_sv_info$dnSV[(dad_sv_info$status == "Probands") & (dad_sv_info$dnSV == FALSE)]),
  length(dad_sv_info$dnSV[(dad_sv_info$status == "Probands") & (dad_sv_info$dnSV == TRUE)])
)
p = ggplot(dad_sv_info[dad_sv_info$status != "NA",], aes(
  x=status,
  y=father_age_birth_years, 
  color=dnSV,fill=dnSV
)) + 
  # scale_y_continuous(limits = c(10, 80)) +
  geom_violin(trim=FALSE,width=.75 ) + 
  labs(y= "Father's age at offspring birth", x = "Sample group")+#, title="Comparison of father's ages between samples with/without paternal dnSVs") +
  annotate("text", x=.84, y=10, label = paste("N=",dad_lengths[1], sep="")) +
  annotate("text", x=1.23, y=10, label = paste("N=",dad_lengths[2], sep="")) +
  annotate("text", x=1.84, y=10, label = paste("N=",dad_lengths[3], sep="")) +
  annotate("text", x=2.23, y=10, label = paste("N=",dad_lengths[4], sep="")) +
  annotate("text", x=1, y=68, label = paste("P-value=",round(dad_sib_test_stat_age$p.value,3), sep=""))+
  annotate("text", x=2, y=68, label = paste("P-value=",round(dad_proband_test_stat_age$p.value,3), sep=""))+
  theme_cowplot()

p = p + geom_segment(aes(x=.8, y=65, xend=1.2, yend=65),linetype='dotted',color='black')
p = p + geom_segment(aes(x=.8, y=65, xend=.8, yend=62),linetype='dotted',color='black')
p = p + geom_segment(aes(x=1.2, y=65, xend=1.2, yend=50),linetype='dotted',color='black')

p = p + geom_segment(aes(x=1.8, y=65, xend=2.2, yend=65),linetype='dotted',color='black')
p = p + geom_segment(aes(x=1.8, y=65, xend=1.8, yend=61),linetype='dotted',color='black')
p = p + geom_segment(aes(x=2.2, y=65, xend=2.2, yend=56),linetype='dotted',color='black')

p
# ggsave("/Users/jon/Research/scripts/de_novo_sv/sfari_denovos/r_plotting/plots/violins-dad-phased.jpg", width = 9, height = 4, dpi=200)
ggsave("/Users/jon/Desktop/supp-fig-4A.svg", width = 9, height = 4, dpi=200)



##################mom plot
mom_lengths=c(
  length(mom_sv_info$dnSV[(mom_sv_info$status == "Unaffected") & (mom_sv_info$dnSV == FALSE)]),
  length(mom_sv_info$dnSV[(mom_sv_info$status == "Unaffected") & (mom_sv_info$dnSV == TRUE)]),
  length(mom_sv_info$dnSV[(mom_sv_info$status == "Probands") & (mom_sv_info$dnSV == FALSE)]),
  length(mom_sv_info$dnSV[(mom_sv_info$status == "Probands") & (mom_sv_info$dnSV == TRUE)])
)
p = ggplot(mom_sv_info[mom_sv_info$status != "NA",], aes(
  x=status,
  y=mother_age_birth_years, 
  color=dnSV,fill=dnSV
)) + 
  geom_violin(trim=FALSE,width=.75 ) + 
  labs(y= "Mother's age at offspring birth", x = "Sample group")+
  annotate("text", x=.84, y=10, label = paste("N=",mom_lengths[1], sep="")) +
  annotate("text", x=1.23, y=10, label = paste("N=",mom_lengths[2], sep="")) +
  annotate("text", x=1.84, y=10, label = paste("N=",mom_lengths[3], sep="")) +
  annotate("text", x=2.23, y=10, label = paste("N=",mom_lengths[4], sep="")) +
  annotate("text", x=1, y=58, label = paste("P-value=",round(mom_sib_test_stat_age$p.value,3), sep=""))+
  annotate("text", x=2, y=58, label = paste("P-value=",round(mom_proband_test_stat_age$p.value,3), sep=""))+
  theme_cowplot()

p = p + geom_segment(aes(x=.8, y=55, xend=1.2, yend=55),linetype='dotted',color='black')
p = p + geom_segment(aes(x=.8, y=55, xend=.8, yend=50),linetype='dotted',color='black')
p = p + geom_segment(aes(x=1.2, y=55, xend=1.2, yend=50),linetype='dotted',color='black')

p = p + geom_segment(aes(x=1.8, y=55, xend=2.2, yend=55),linetype='dotted',color='black')
p = p + geom_segment(aes(x=1.8, y=55, xend=1.8, yend=50),linetype='dotted',color='black')
p = p + geom_segment(aes(x=2.2, y=55, xend=2.2, yend=50),linetype='dotted',color='black')

p
ggsave("/Users/jon/Research/scripts/de_novo_sv/sfari_denovos/r_plotting/plots/violins-mom-phased.jpg", width = 9, height = 4, dpi=200)
ggsave("/Users/jon/Desktop/supp-fig-4B.svg", width = 9, height = 4, dpi=200)
######################################################################################################################################################