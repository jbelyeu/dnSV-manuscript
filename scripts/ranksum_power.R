library(ggplot2)
library(gridExtra)
library(cowplot)
library(fitdistrplus)
library(readr)
library(tibble)
library(dplyr)
library(ANOVAreplication)
library(svglite)

#read in data and split the probands from sibs
sv_info=read.csv("../data/dnsv_dataframe.csv", header=TRUE, sep=",")
sv_info$status = recode(sv_info$role, "sib"="Unaffected","proband"="Probands", "mother"="Unaffected","father"="Unaffected")

sv_info = sv_info[(sv_info$father != ".") & (sv_info$mother != ".") ,]
sv_info = sv_info[(sv_info$father_age_birth_years >0) & (sv_info$mother_age_birth_years >0),]
sv_info$father_age_birth_years = as.integer(sv_info$father_age_birth_years)
sv_info$mother_age_birth_years = as.integer(sv_info$mother_age_birth_years)
sv_info$has_dnSV = factor(sv_info$all_sv > 0)
sv_info$dnSV = factor(sv_info$all_sv>0)

# ####################################################################################################
# #dn SNV comparison
# #Mann-Whitney-Wilcoxon for count data with different group sizes
# snv_sv_info = sv_info[(sv_info$all_sv > 0)  & (!is.na(sv_info$dn_snvs)) & (sv_info$father_age_birth_years > 0),]
# 
# snv_counts_with_dnsv = as.numeric(snv_sv_info$dn_snvs[(snv_sv_info$dnSV == TRUE) & (snv_sv_info$status == "Unaffected") ])
# snv_counts_without_dnsv = as.numeric(snv_sv_info$dn_snvs[(snv_sv_info$dnSV == FALSE) & (snv_sv_info$status == "Unaffected") ])
# 
# sib_snv_test_stat_dnsnv = wilcox.test(
#   snv_counts_with_dnsv,
#   snv_counts_without_dnsv
# )
# proband_snv_test_stat_dnsnv = wilcox.test(
#   snv_sv_info$dn_snvs[(snv_sv_info$dnSV == TRUE) & (snv_sv_info$status == "Probands") ],
#   snv_sv_info$dn_snvs[(snv_sv_info$dnSV == FALSE) & (snv_sv_info$status == "Probands") ]
# )
# 
# 
# lengths=c(
#   length(snv_sv_info$dnSV[(snv_sv_info$status == "Probands") & (snv_sv_info$dnSV == FALSE) ]),
#   length(snv_sv_info$dnSV[(snv_sv_info$status == "Probands") & (snv_sv_info$dnSV == TRUE) ]),
#   length(snv_sv_info$dnSV[(snv_sv_info$status == "Siblings") & (snv_sv_info$dnSV == FALSE) ]),
#   length(snv_sv_info$dnSV[(snv_sv_info$status == "Siblings") & (snv_sv_info$dnSV == TRUE) ])
# )
# 
# p = ggplot(snv_sv_info, aes(
#   x=status,
#   y=dn_snvs,
#   color=dnSV,
#   fill=dnSV
# )) +
#   geom_violin(trim=FALSE,width=.75 ) +
#   labs(y= "dnSNV count", x = "Sample group")+
#   annotate("text", x=.84, y=0, label = paste("N=",lengths[1], sep="")) +
#   annotate("text", x=1.23, y=0, label = paste("N=",lengths[2], sep="")) +
#   annotate("text", x=1.84, y=0, label = paste("N=",lengths[3], sep="")) +
#   annotate("text", x=2.23, y=0, label = paste("N=",lengths[4], sep="")) +
#   annotate("text", x=1, y=260, label = paste("P-value=",round(sib_snv_test_stat_dnsnv$p.value,3), sep=""))+
#   annotate("text", x=2, y=260, label = paste("P-value=",round(proband_snv_test_stat_dnsnv$p.value,3), sep=""))+
#   theme_cowplot()
# 
# p = p + geom_segment(aes(x=.8, y=250, xend=1.2, yend=250),linetype='dotted',color='black')
# p = p + geom_segment(aes(x=.8, y=250, xend=.8, yend=200),linetype='dotted',color='black')
# p = p + geom_segment(aes(x=1.2, y=250, xend=1.2, yend=245),linetype='dotted',color='black')
# 
# p = p + geom_segment(aes(x=1.8, y=250, xend=2.2, yend=250),linetype='dotted',color='black')
# p = p + geom_segment(aes(x=1.8, y=250, xend=1.8, yend=220),linetype='dotted',color='black')
# p = p + geom_segment(aes(x=2.2, y=250, xend=2.2, yend=170),linetype='dotted',color='black')
# 
# p = p + scale_fill_manual(values=c("blue", "green"))
# p = p + scale_color_manual(values=c("blue", "green"))
# 
# p
# ggsave("../plots/violins_dnsnv.svg", width = 8, height = 4, dpi=200)
####################################################################################################################################

#Age comparison
####################################################################################################
#"greater" tests whether X(has dnSV) is right-shifted vs y (doesn't have dnSV)
sib_test_stat_age = wilcox.test(alternative="greater",
  sv_info$father_age_birth_years[(sv_info$dnSV == TRUE) & (sv_info$status == "Unaffected") ], 
  sv_info$father_age_birth_years[(sv_info$dnSV == FALSE) & (sv_info$status == "Unaffected")]
)
#"greater" tests whether X(has dnSV) is right-shifted vs y (doesn't have dnSV)
proband_test_stat_age = wilcox.test(alternative="greater",
  sv_info$father_age_birth_years[(sv_info$dnSV == TRUE) & (sv_info$status == "Probands") ], 
  sv_info$father_age_birth_years[(sv_info$dnSV == FALSE) & (sv_info$status == "Probands")]
)

lengths=c(
  length(sv_info$dnSV[(sv_info$status == "Unaffected") & (sv_info$dnSV == FALSE)]),
  length(sv_info$dnSV[(sv_info$status == "Unaffected") & (sv_info$dnSV == TRUE)]),
  length(sv_info$dnSV[(sv_info$status == "Probands") & (sv_info$dnSV == FALSE)]),
  length(sv_info$dnSV[(sv_info$status == "Probands") & (sv_info$dnSV == TRUE)])
)

#current version of the violins, includes overlay
ggplot(sv_info[sv_info$status != "NA",], aes(
  x=status,
  y=father_age_birth_years, 
  color=dnSV,
  fill=dnSV
)) + 
  theme_cowplot() +
  # geom_violin(trim=FALSE, width=.6, position = position_dodge(1.2)) + 
  # geom_violin(trim=FALSE, width = .3, position = position_identity(), alpha=0, color = NA)+
  # geom_violin(trim=FALSE, width = .3, position = position_identity(), alpha=0)+
  geom_violin(trim=FALSE, width=.5) + 
  geom_violin(trim=FALSE, width = .25, position = position_nudge(0.4), alpha=0)+
  labs(y= "Father's age at offspring birth", x = "Sample group")+
  annotate("text", x=.87, y=8, label = paste("N=",lengths[1], sep=""), size=5) +
  annotate("text", x=1.125, y=8, label = paste("N=",lengths[2], sep=""), size=5) +
  annotate("text", x=1.4, y=10, label = paste("Unaffected\noverlay", sep=""), size=5) +
  annotate("text", x=1.87, y=8, label = paste("N=",lengths[3], sep=""), size=5) +
  annotate("text", x=2.125, y=8, label = paste("N=",lengths[4], sep=""), size=5) +
  annotate("text", x=2.4, y=10, label = paste("Probands\noverlay", sep=""), size=5) +
  annotate("text", x=1, y=70, label = paste("P-value=",round(sib_test_stat_age$p.value,3), sep=""), size=6)+
  annotate("text", x=2, y=70, label = paste("P-value=",round(proband_test_stat_age$p.value,3), sep=""), size=6) +
  
  geom_segment(aes(x=.87, y=67, xend=1.125, yend=67),linetype='dotted',color='black') +
  geom_segment(aes(x=.87, y=67, xend=.87, yend=62),linetype='dotted',color='black') +
  geom_segment(aes(x=1.125, y=67, xend=1.125, yend=62),linetype='dotted',color='black') +
  geom_segment(aes(x=1.87, y=67, xend=2.125, yend=67),linetype='dotted',color='black') +
  geom_segment(aes(x=1.87, y=67, xend=1.87, yend=62),linetype='dotted',color='black') +
  geom_segment(aes(x=2.125, y=67, xend=2.125, yend=62),linetype='dotted',color='black') +
  theme(
    legend.position=c(0.93,0.85),
    legend.key.size = unit(.8, "cm"),
    legend.text = element_text(size=15),
    axis.text =  element_text(size=16),
    axis.title =  element_text(size=18)
  )
  ggsave("../plots/violins_age.png", width = 10, height = 5, dpi=200)
  
#
#
#
#
#
#
#
#
#
#
#
# power analysis
library(pwr)
powers = c()
effect_sizes = c()
for (i in seq(0,1, by=0.01)){
  pwr_test = (pwr.t2n.test(n1=lengths[1],n2=lengths[2],d=i,sig.level=0.05))
  powers = append(powers, pwr_test$power)
  effect_sizes=append(effect_sizes,i)
}

powers = c()
effect_sizes = c()
for (i in seq(0,1, by=0.01)){
  pwr_test = (pwr.t2n.test(n1=lengths[3],n2=lengths[4],d=i,sig.level=0.05))
  powers = append(powers, pwr_test$power)
  effect_sizes=append(effect_sizes,i)
}

sib_power_analysis = data.frame("Power"=powers,"Effect.size"=effect_sizes)
sibling_ef = pwr.t2n.test(n1=lengths[1],n2=lengths[2],power=0.8,sig.level=0.05)$d
proband_power_analysis = data.frame("Power"=powers,"Effect.size"=effect_sizes)
proband_ef = pwr.t2n.test(n1=lengths[3],n2=lengths[4],power=0.8,sig.level=0.05)$d

pooled_sd = pooled.sd(data.frame(ages=sv_info$father_age_birth_years, group=sv_info$has_dnSV))
sibs_ef_years = round(pooled_sd*sibling_ef,3)
probands_ef_years = round(pooled_sd*proband_ef,3)

ggplot(proband_power_analysis, aes(x=Effect.size, y=Power))+
  geom_segment(aes(x=proband_ef, y=0, xend=proband_ef, yend=1),color='#29AAE2',linetype="dotted", size=1.2)+
  geom_segment(aes(x=sibling_ef, y=0, xend=sibling_ef, yend=1),color='#FBB03A',linetype="dotted", size=1.2)+
  geom_line(color="black", size=1) +
  geom_point(color="black", size=2)+
  geom_line(data=sib_power_analysis,color="black",aes(x=Effect.size, y=Power),size=1)+
  geom_point(data=sib_power_analysis,color="black", aes(x=Effect.size, y=Power),size=2) +
  labs(y= "Power", x = "Effect size", fontsize=16)+#, title="Power analysis of father's age effect on dnSV risk")+
  annotate("text", x = .55, y = (.20), size=7, label = paste('Unaffected Samples (',sibs_ef_years,' years, d=',round(sibling_ef,3),')',sep=""), colour="#FBB03A") +
  annotate("text", x = .47, y = (.10), size=7, label = paste('Probands (',probands_ef_years,' years, d=',round(proband_ef,3),')',sep=""), colour="#29AAE2")+
  theme_cowplot() +
  theme(axis.text=element_text(size=18),axis.title=element_text(size=18))
  
ggsave("../plots/pwr.svg", dpi=250, width=10, height=5)
#
#
#
#
#
#
#
#
#
#dnSNV counts vs paternal age
adjust_ar=FALSE
alpha=.4
proband_snv_info = sv_info[(!is.na(sv_info$dn_snvs)) &(sv_info$dn_snvs>0) &(sv_info$dn_snvs<125) & (sv_info$status=='Probands'),]
proband_model = glm(dn_snvs ~ father_age_birth_years, data=proband_snv_info, family = poisson(link="identity"))
proband_pred = predict(proband_model, type='response', se.fit=TRUE)
#CIs
dist=1.96
proband_snv_info$dad_ci_lo = proband_pred$fit - dist * proband_pred$se.fit
proband_snv_info$dad_ci_hi = proband_pred$fit + dist * proband_pred$se.fit

sib_snv_info = sv_info[(!is.na(sv_info$dn_snvs)) & (sv_info$dn_snvs>0) &(sv_info$dn_snvs<125) & (sv_info$status=='Unaffected'),]
sib_model = glm(dn_snvs ~ father_age_birth_years, data=sib_snv_info, family = poisson(link="identity"))
sib_pred = predict(sib_model, type='response', se.fit=TRUE)
#CIs
dist=1.96
sib_snv_info$dad_ci_lo = sib_pred$fit - dist * sib_pred$se.fit
sib_snv_info$dad_ci_hi = sib_pred$fit + dist * sib_pred$se.fit

# get min and max X and Y values for plot limits
min_age = min(min(sib_snv_info$father_age_birth_years), min(proband_snv_info$father_age_birth_years))
max_age = max(max(sib_snv_info$father_age_birth_years), max(proband_snv_info$father_age_birth_years))
min_dnm = min(min(sib_snv_info$dn_snvs),min(proband_snv_info$dn_snvs))
max_dnm = max(max(sib_snv_info$dn_snvs),max(proband_snv_info$dn_snvs))


# set the upper Y limit 
if (max_dnm < 15) {
  max_dnm = max_dnm
} else {
  max_dnm = max_dnm + 15
}

# adjust the aspect ratio if needed.
# these adjustments are specific to plotting either second-generation
# DNMs or gonosomal DNMs, and are for aesthetic purposes only
if (adjust_ar) {
  adjust = (0.075 * min_age/min_dnm)
}else {
  adjust = 2.25
}

proband_color = "blue"
sib_color = "darkorange"

ggplot(proband_snv_info[proband_snv_info$dn_snvs > 0,]) + 
  # plot the raw data - probands
  geom_jitter(data=proband_snv_info, aes(x=father_age_birth_years, y=dn_snvs, alpha=0.8),size=1, pch=21, fill=proband_color, col='white', stroke=0.05, show.legend = FALSE) +
  # plot the predictions from the fitted GLM - probands
  geom_line(data=cbind(proband_snv_info, pred_d=proband_pred$fit), aes(x=father_age_birth_years, y=pred_d), col=proband_color) +
  # plot confidence bands - probands
  geom_ribbon(data=proband_snv_info,aes(x=father_age_birth_years, ymin=dad_ci_lo, ymax=dad_ci_hi), alpha=alpha, fill=proband_color, show.legend = FALSE) +
  # plot the raw data - sibs
  geom_jitter(data=sib_snv_info, aes(x=father_age_birth_years, y=dn_snvs, alpha=0.8), size=1, pch=21, fill=sib_color, col='white', stroke=0.05, show.legend = FALSE) +
  # plot the predictions from the fitted GLM - sibs
  geom_line(data=cbind(sib_snv_info, pred_d=sib_pred$fit), aes(x=father_age_birth_years, y=pred_d), col=sib_color) +
  # plot confidence bands - sibs
  geom_ribbon(data=sib_snv_info, aes(x=father_age_birth_years, ymin=dad_ci_lo, ymax=dad_ci_hi), alpha=alpha, fill=sib_color, show.legend = FALSE) +
  theme_cowplot()+
  annotate("text", x=20.45, y=110, label = paste("Proband coef=",round(coef(proband_model)[2], digits=3), sep=""), color=proband_color) +
  annotate("text", x=20, y=105, label = paste("Sibling coef=",round(coef(sib_model)[2], digits=3), sep=""), color=sib_color) +
  labs(x= "Father's age at offspring birth", y = "Number of dnSNVs")+#, title="Correlation between paternal age and dnSNV count ") 
ggsave("../plots/snv_age.svg", width = 8, height = 4, dpi=200)

