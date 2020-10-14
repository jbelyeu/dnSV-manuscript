library(ggplot2)
library(gridExtra)
library(cowplot)
library(fitdistrplus)
library(readr)
library(tibble)
library(dplyr)
library(svglite)
library(ANOVAreplication)


#read in data and split the probands from sibs
setwd("/Users/jon/Research/scripts/de_novo_sv/sfari_denovos/r_plotting/data_files/")
sv_info=read.csv("dnsv_dataframe.csv", header=TRUE, sep=",")
sv_info$status = recode(sv_info$role, "sib"="Unaffected","proband"="Probands", "mother"="Unaffected","father"="Unaffected")

# sv_info = sv_info[(sv_info$father != "8440")  ,]
# sv_info = sv_info[(sv_info$sample != "8412")  ,]
# sv_info = sv_info[(sv_info$project == "SFARI")  ,]
# sv_info = sv_info[(sv_info$sample != "8441")  ,]


sv_info = sv_info[(sv_info$status == "Probands") | (sv_info$status == "Unaffected") ,]
sv_info = sv_info[sv_info$father_age_birth_years >0,]
sv_info$father_age_birth_years = as.integer(sv_info$father_age_birth_years)
sv_info$mother_age_birth_years = as.integer(sv_info$mother_age_birth_years)
sv_info$has_dnSV = factor(sv_info$all_sv > 0)
sv_info$dnSV = factor(sv_info$all_sv>0)
#sv_info = sv_info[sv_info$project == "SFARI",]

svtypes = c("MEI","DUP", "DEL")
for (i in 1:length(svtypes)) {
  svtype = svtypes[i]
  if (svtype == "MEI") {
    typed_sv_info = sv_info[(sv_info$MEI > 0) | (sv_info$dnSV == FALSE),]
  } else if (svtype == "DUP") {
    typed_sv_info = sv_info[(sv_info$DUP > 0) | (sv_info$dnSV == FALSE),]
  } else if (svtype == "DEL") {
    typed_sv_info = sv_info[(sv_info$DEL > 0) | (sv_info$dnSV == FALSE),]
  }
  
  sib_test_stat = wilcox.test(alternative="greater",
    typed_sv_info$father_age_birth_years[(typed_sv_info$dnSV == TRUE) & (typed_sv_info$status == "Unaffected") ], 
    typed_sv_info$father_age_birth_years[(typed_sv_info$dnSV == FALSE) & (typed_sv_info$status == "Unaffected")]
  )
  proband_test_stat = wilcox.test(alternative="greater",
    sv_info$father_age_birth_years[(typed_sv_info$dnSV == TRUE) & (typed_sv_info$status == "Probands") ], 
    sv_info$father_age_birth_years[(typed_sv_info$dnSV == FALSE) & (typed_sv_info$status == "Probands")]
  )
  lengths=c(
    length(typed_sv_info$dnSV[(typed_sv_info$status == "Unaffected") & (typed_sv_info$dnSV == FALSE)]),
    length(typed_sv_info$dnSV[(typed_sv_info$status == "Unaffected") & (typed_sv_info$dnSV == TRUE)]),
    length(typed_sv_info$dnSV[(typed_sv_info$status == "Probands") & (typed_sv_info$dnSV == FALSE)]),
    length(typed_sv_info$dnSV[(typed_sv_info$status == "Probands") & (typed_sv_info$dnSV == TRUE)])
  )
  p = ggplot(typed_sv_info[typed_sv_info$status != "NA",], aes(
    x=status,
    y=father_age_birth_years, 
    color=dnSV,fill=dnSV
  )) + 
    # scale_y_continuous(limits = c(10, 80)) +
    geom_violin(trim=FALSE,width=.75 ) + 
    # labs(y= "Father's age at offspring birth", x = "Sample group", title=paste("Comparison of father's ages between samples with/without dn",svtype,"s", sep="")) +
    labs(y= "Father's age at offspring birth", x = "Sample group") +
    annotate("text", x=.84, y=10, label = paste("N=",lengths[1], sep="")) +
    annotate("text", x=1.23, y=10, label = paste("N=",lengths[2], sep="")) +
    annotate("text", x=1.84, y=10, label = paste("N=",lengths[3], sep="")) +
    annotate("text", x=2.23, y=10, label = paste("N=",lengths[4], sep="")) +
    annotate("text", x=1, y=68, label = paste("P-value=",round(sib_test_stat$p.value,3), sep=""))+
    annotate("text", x=2, y=68, label = paste("P-value=",round(proband_test_stat$p.value,3), sep=""))+
    theme_cowplot()
  
  p = p + geom_segment(aes(x=.8, y=65, xend=1.2, yend=65),linetype='dotted',color='black')
  p = p + geom_segment(aes(x=.8, y=65, xend=.8, yend=60),linetype='dotted',color='black')
  p = p + geom_segment(aes(x=1.2, y=65, xend=1.2, yend=60),linetype='dotted',color='black')
  
  p = p + geom_segment(aes(x=1.8, y=65, xend=2.2, yend=65),linetype='dotted',color='black')
  p = p + geom_segment(aes(x=1.8, y=65, xend=1.8, yend=60),linetype='dotted',color='black')
  p = p + geom_segment(aes(x=2.2, y=65, xend=2.2, yend=60),linetype='dotted',color='black')
  
  p = p + scale_fill_manual(values=c("grey", "black"))
  p = p + scale_color_manual(values=c("grey", "black"))
  p
  ggsave(paste("/Users/jon/Research/scripts/de_novo_sv/sfari_denovos/r_plotting/plots/",svtype,"_violins.svg",sep=""), width = 8, height = 4, dpi=200)
  ggsave(paste("/Users/jon/Desktop/",svtype,"_violins.png"), width = 8, height = 4, dpi=200)
}
