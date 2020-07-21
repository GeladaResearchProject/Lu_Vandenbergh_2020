## Packages needed
library(ggplot2)
library(ggridges)
library(gridExtras)
library(lme4)
library(lmerTest)
library(visreg)
library(dplyr)
library(survival)
setwd("~/Desktop/Data/RData_Vandenbergh_2020")

## FIGURE 1 -- Maturations in days relative to takeover
figure1<-ggplot(maturations, aes(x=TimeDiff)) + geom_histogram(color="black", fill="blue", alpha=0.3, breaks=(seq(-180, 180, 30)), closed="left") + 
  theme_classic() + xlab("Days relative to male takeover") + 
  ylab("Females maturing") + theme(axis.title = element_text(face="bold")) + 
  scale_y_continuous(breaks=seq(0, 12, 2), labels=round(seq(0, 12, 2)), limits=c(0,10))
figure1 + scale_x_continuous(breaks=seq(-180,180, 60))

## ANALYSIS: LIKELIHOOD OF MATURATION
## Survival model (COXPH)
## Likelihood of maturation across female-months from 3.4 years to maturation
## Matured ~ Takeover + Father Presence + Ecological variables (cluster on female identity)
## Time transformation on father presence
surv_object<-Surv(time = survival.data$Start.Age, time2=survival.data$End.Age, event=survival.data$Matured)
fit.coxph1 <- coxph(surv_object ~ Takeover.Range + tt(Father.Presence) + scale(Cumulative.Rainfall) + scale(Avg.Min) + cluster(ID),  data = survival.data)
summary(fit.coxph1)

## FIGURE 2 -- ridgeplot of female ages at maturation
## Fathers absent or present at 3.4 years
## Takeover in the 3 months before maturation
library(ggridges)

maturations$Category[maturations$Takeover.Maturation=="N"]<-"Father absent,\nno takeover"
maturations$Category[maturations$Takeover.Maturation=="Y"]<-"Father absent,\ntakeover"
maturations$Category[maturations$Father.Presence=="Y"]<-"Father present,\nno takeover" 
maturations$Category[maturations$Father.Presence=="Y" & maturations$Takeover.Maturation=="Y"]<-"Father present,\ntakeover"
maturations$Category<-factor(maturations$Category, levels=c("Father absent,\ntakeover", "Father absent,\nno takeover", "Father present,\ntakeover", "Father present,\nno takeover"))

colors<-c("#1F78B4","#A6CEE3","#33A02C","#B2DF8A")

figure2<-ggplot(data=maturations, aes(x = MaturationAge, y = Category, fill=Category)) + scale_fill_manual(values=colors) + 
  stat_density_ridges(alpha=0.7,scale=0.9, jittered_points=T, point_shape="|",point_size=2, position=position_points_jitter(height=0)) + 
  theme_classic() + theme(legend.position = "none") + xlab("Maturation age (years)") + ylab("") + 
  theme(axis.title = element_text(face="bold"), legend.title = element_text(face="bold")) + 
  scale_x_continuous(breaks=seq(3.0, 7, 1))
figure2 <-figure2 +  scale_y_discrete(expand=c(0.04,0)) + coord_cartesian(clip="off")
figure2

## ANALYSIS: MATURATION AGE BY FATHER PRESENCE AND TAKEOVERS
## Linear mixed model (LMM)
## Age at maturation as a function of father presence and takeovers
## Maturation Age ~ Father Presence (y/n) + Takeover Maturation (y/n)
model<-lmer(data=maturations, MaturationAge~ Father.Presence + Takeover.Maturation + (1|Unit))
summary(model)
median(maturations$MaturationAge[maturations$Takeover.Maturation=="Y"])
median(maturations$MaturationAge[maturations$Takeover.Maturation=="N"])
median(maturations$MaturationAge[maturations$Father.Presence=="Y"])
median(maturations$MaturationAge[maturations$Father.Presence=="N"])

## ANALYSIS: TAKEOVER-MATURATIONS BY AGE
## Generalized linear mixed model (GLMM)
## Likelihood of maturing following takeover as a function of age
## Matured (y/n) ~ Age at takeover
binomial.model<-glmer(data=takeovers.binomial, Matured ~ TakeoverAge + (1|ID), family="binomial")
summary(binomial.model)

## ANALYSIS: TAKEOVER-RELATED ESTROGENS
## Linear mixed model (LMM)
## Outcome: logged estrogens
## Fixed effects: Takeover (30 days), Before Maturation (100 days), Female Age
## Rainfall (90 days), Avg min temp (30 days), Wash step, Laboratory
estrogen.model<-lmer(data=imm.gelada.estrogens, log(E2.ng.g) ~ Post.Takeover + Just.Before.Maturation + Age.at.sample + scale(Cumulative.Rainfall) + scale(Avg.Min) + Wash.step + Lab + (1|ID) + (1|Unit), REML=F, control=lmerControl(optimizer="bobyqa"))
estrogen.model.no.takeovers<-lmer(data=imm.gelada.estrogens, log(E2.ng.g) ~ Just.Before.Maturation + Age.at.sample + scale(Cumulative.Rainfall) + scale(Avg.Min) + Wash.step + Lab + (1|ID) + (1|Unit), REML=F, control=lmerControl(optimizer="bobyqa"))
summary(estrogen.model)
summary(estrogen.model.no.takeovers)
imm.gelada.estrogens$Residuals<-residuals(estrogen.model.no.takeovers)

## FIGURE 3a-- INDIVIDUAL ESTROGEN TRAJECTORIES
## FEMALE 1 PLOT
## Use estrogen.examples dataset
plot1<-ggplot(data=estrogen.examples[which(estrogen.examples$ID=="ID_12"),], aes(x=(Age.at.sample), y=round(E2.ng.g, digits = 1))) +geom_line() +theme_classic() + geom_vline(xintercept=estrogen.examples$Age.at.Takeover[which(estrogen.examples$ID=="ID_12")], lty=4) +geom_point() + scale_x_continuous(breaks=seq(3, 4.5, by=.25), limits = c(4.13-.8,4.13+0.5)) + xlab("") + ylab("") + geom_text(label="Female 1", x=3.5, y=20)
plot1<-plot1 + scale_y_continuous(labels = scales::number_format(accuracy = 0.1), limits=c(0,22)) + coord_cartesian(clip="off")
plot1

## FEMALE 2 PLOT
plot2<-ggplot(data=estrogen.examples[which(estrogen.examples$ID=="ID_1"),], aes(x=(Age.at.sample), y=E2.ng.g)) +geom_line() +theme_classic() + geom_vline(xintercept=estrogen.examples$Age.at.Takeover[which(estrogen.examples$ID=="ID_1")], lty=4) +geom_point() + scale_x_continuous(breaks=seq(1.5, 3.0, by=.25), limits = c(2.39-.8,2.39+0.5)) + ylim(0,10) + ylab("Fecal estrogens (ng/g)") + xlab("") + theme(axis.title=element_text(size=11,face="bold")) + geom_text(label="Female 2", x=1.75, y=10) + coord_cartesian(clip="off")
plot2
  
## FEMALE 3 PLOT
plot3<-ggplot(data=estrogen.examples[which(estrogen.examples$ID=="ID_5"),], aes(x=(Age.at.sample), y=E2.ng.g)) +geom_line() +theme_classic() + geom_vline(xintercept=estrogen.examples$Age.at.Takeover[which(estrogen.examples$ID=="ID_5")], lty=4) +geom_point() + scale_x_continuous(breaks=seq(.5, 2, by=.25), limits = c(1.37-.8,1.37+0.5)) + xlab("Age (years)") + ylim(0,10) + theme(axis.title=element_text(size=11,face="bold")) + ylab("") + geom_text(label="Female 3", x=0.75, y=10) + coord_cartesian(clip="off")
plot3

figure3<-grid.arrange(rbind(ggplotGrob(plot1), ggplotGrob(plot2), ggplotGrob(plot3), size = "last"))

## TAKEOVER BY AGE (FIGURE 3b)
## Use all data -- imm.estrogens
std <- function(x) sd(x)/sqrt(length(x))
takeover.summary.age<-imm.gelada.estrogens %>%
  group_by(Post.Takeover, Rounded.age) %>%
  summarise(Mean=mean(Residuals), SE=std(Residuals), count=n())
takeover.summary.age$CI<-1.96*takeover.summary.age$SE

takeover.summary.age<-droplevels(takeover.summary.age)

age<-c("1-yr olds", "2-yr olds", "3-yr olds", "4-yr olds")
names(age) <- c(1,2,3,4)

takeover.summary.age<-takeover.summary.age[grepl("1|2|3|4",takeover.summary.age$Rounded.age),]
takeover.summary.age$Post.Takeover<- recode(takeover.summary.age$Post.Takeover, "Takeover"="Yes", "No Takeover"="No")
takeover.summary.age$Post.Takeover[takeover.summary.age$Post.Takeover=="No Takeover"] <- "No"
takeover.summary.age$Post.Takeover[takeover.summary.age$Post.Takeover=="Takeover"] <- "Yes"

figure3b<-ggplot(data=takeover.summary.age, aes(x=Post.Takeover, y=Mean)) + facet_wrap(.~Rounded.age, labeller=labeller(Rounded.age=age), ncol=2) + geom_errorbar(ymin=takeover.summary.age$Mean-takeover.summary.age$CI, ymax=takeover.summary.age$Mean+takeover.summary.age$CI, width=0)+
  geom_crossbar(ymin=takeover.summary.age$Mean-takeover.summary.age$SE, ymax=takeover.summary.age$Mean+takeover.summary.age$SE, fill=rep(c("#A6CEE3","#1F78B4"), each=4), width=0.4, fatten=1.5) + ylab("Fecal estrogen residuals") + xlab("Takeover?") +
  theme_test() + theme(axis.title = element_text(size=11, face="bold"), axis.text =  element_text(size=9), strip.text = element_text(size=11, face="bold")) + ylim(-0.1, 1.2)                                                                                                                           
figure3b

## ANALYSIS: AGE AT FIRST BIRTH
## Linear model
## Outcome: age at first birth
## Fixed effect: age at maturation
AFR<-maturations[!is.na(maturations$AFR),]
model<-lm(data=AFR, AFR ~ MaturationAge)
summary(model)

## ANALYSIS: TIME TO FIRST POST-MATURATION TAKEOVER
## Linear mixed model (LMM)
## Outcome: Days to first takeover following maturation
## Fixed effect: male-mediated maturation? (y/n)
time.to.takeovers<-lmer(data=first.takeovers, as.numeric(TimeDiff) ~ TakeoverMaturation + (1|Unit))
summary(time.to.takeovers)


## SUPPLEMENTARY PLOTS
## FIGURE S1
maturation.summary<-estrogen.maturations %>%
  group_by(Window) %>%
  summarise(Mean=mean(log(E2.ng.g)), SE=std(log(E2.ng.g)), count=n())
maturation.summary$CI<-1.96*maturation.summary$SE

labels<-c(maturation.summary$count)
figureS1<-ggplot(data=maturation.summary, aes(x=Window, y=Mean)) +geom_errorbar(ymin=maturation.summary$Mean-maturation.summary$CI, ymax=maturation.summary$Mean+maturation.summary$CI, width=15)+ geom_point() + geom_line()  + xlab("Days relative to maturation") + ylab("Fecal estrogens log(ng/g)") +
  theme_classic() + theme(axis.title = element_text(size=11, face="bold"), axis.text =  element_text(size=10)) + geom_vline(xintercept=0, linetype=2, size=.6) +ylim(0, 1.5) + xlim(-160,50)
figureS1<-figureS1+geom_text(aes(y=maturation.summary$Mean-maturation.summary$CI,label = maturation.summary$count), size = 3, vjust=3, fontface="italic")
figureS1

## FIGURE S2 
frame1<-as.data.frame(exp(fit.coxph1$coefficients))
frame2<-as.data.frame(exp(confint(fit.coxph1)))
factors<-c("Takeover\n(prev. 3 months)", "Father presence", "Cumulative rainfall\n(prev. 3 months)", "Avg. min. temp.\n(prev. month)")
estimates<-cbind(factors, frame1,frame2)
names(estimates)[2]<-"Estimate"
names(estimates)[3]<-"LL"
names(estimates)[4]<-"UL"

cols<-c("#1F78B4", "#33A02C", 'black', "black")

figureS2 <- ggplot(data=estimates,
                   aes(x = factors,y = Estimate, ymin = LL, ymax = UL ))+ theme_classic() +
  geom_pointrange(colour=cols, cex=0.4) + geom_hline(yintercept =1, lty="dotted") + 
  xlab('')+ ylab("Hazard Ratio (95% Confidence Interval)")+ scale_y_continuous(breaks=0:6) +
  geom_errorbar(aes(ymin=LL, ymax=UL),width=0.1,cex=1, colour=cols) +   theme(plot.title=element_text(size=11,face="bold"),
                                                                              axis.text.x=element_text(face="bold"),
                                                                              axis.title=element_text(size=11,face="bold"),
                                                                              strip.text = element_text(hjust=0,vjust = 1,angle=180,face="bold"))+  coord_flip() 

figureS2

## FIGURE S3 
## Likelihood of maturing following takeover as a function of age
## Visualization of the model above
figureS3<-visreg(binomial.model, scale="response", gg=T) 
figureS3+theme_classic() + ylab("Probability of maturing") + 
  xlab("Age at takeover (months)") + 
  theme(axis.title = element_text(size=11, face="bold")) + 
  scale_x_continuous(limits=c(3,6))

## FIGURE S4
## Relationship between maturation age and age at first birth
figureS4<-ggplot(data=AFR, aes(x=MaturationAge, y=AFR)) + scale_color_manual(values=cols) + theme_classic() + xlab("Age at maturation (years)") + ylab("Age at first birth (years)") + theme(axis.title = element_text(size=11, face="bold"))
figureS4<-figureS4 + geom_point() + geom_smooth(se=T, method="lm", color="black") + theme(legend.position = "none") + scale_x_continuous(limits = c(3,7))
figureS4

