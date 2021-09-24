rm(list=ls())
#Prep PEDI data 
library(readxl)
library(dplyr)
library(data.table)

dat <- read.csv("harmonized_brain_data.csv")[,-1]

#Diffs in motion
summary(lm(dat1$motion ~ dat1$Diagnosis + dat1$Age +dat1$Study))$coefficients

#Node based analyses
p <- c()
graph.anal <- data.frame()
graph.AN <- data.frame()
roi.fr <- data.frame()
p.perm <- data.frame()
sig <- c()
perm <- data.frame(n=seq(1,10000,1))
graph$subname <- NULL
Effects <- data.frame()

test <- colnames(graph)[3:86]
for(i in 1:length(test)){
  outcome <- test[i]
  mod <- summary(lm(dat1[[outcome]] ~ dat1$Diagnosis + dat1$Age + dat1$motion + dat1$Study))$coefficients
  fr <- as.data.frame(mod)
  fr$outcome <- outcome
  p <- c(p, fr$`Pr(>|t|)`[2])
  print(fr[2,])
  if(fr$`Pr(>|t|)`[2] < 0.05){
    graph.anal <- rbind(graph.anal, fr)
    sig <- c(ROI[i])
  }
  set.seed(12345)
  m <- permuco::lmperm(lm(dat1[[outcome]] ~ dat1$Diagnosis + dat1$Age + dat1$motion + dat1$Study), np = 10000)
  perm <- cbind(perm, m$distribution[,1])
  roi <- c(ROI[i], fr$Estimate[2], 'Efficiency')
  roi <- as.data.frame(roi)
  roi <- t(roi)
  colnames(roi) <- c("connection", "effect", "type")
  roi.fr <- rbind(roi.fr, roi)
  fr <- data.frame(fr)
  fr <- fr[2,]
  colnames(fr) <- c("Estimate",   "Std..Error",     "t.value",     "P",  "outcome")
  Effects <- rbind(Effects, fr)
}
PAF <- data.frame()
perm <- perm[,-1]
pa <- compute_minP(perm, alternative = "two.sided")
paf <- data.frame(p=pa$main$pvalue, name=roi.fr$connection)
paf$name <- paste0("EF_",paf$name)
p.perm <- rbind(p.perm,paf)
paf <- subset(paf, paf$p < 0.05)
PAF <- rbind(PAF, paf)


p <- p.adjust(p, "fdr")
f <- data.frame(p=p, name=test)
fr <- subset(f, f$p < 0.05)
M <- rbind(M,fr)
p <- c()
perm <- data.frame(n=seq(1,10000,1))
graph.anal <- data.frame()
graph.AN <- data.frame()
test <- colnames(graph)[87:170]
for(i in 1:length(test)){
  outcome <- test[i]
  mod <- summary(lm(dat1[[outcome]] ~ dat1$Diagnosis + dat1$Age + dat1$motion + dat1$Study))$coefficients
  fr <- as.data.frame(mod)
  fr$outcome <- outcome
  p <- c(p, fr$`Pr(>|t|)`[2])
  print(fr[2,])
  if(fr$`Pr(>|t|)`[2] < 0.05){
    graph.anal <- rbind(graph.anal, fr)
    sig <- c(ROI[i])
  }
  set.seed(12345)
  m <- permuco::lmperm(lm(dat1[[outcome]] ~ dat1$Diagnosis + dat1$Age + dat1$motion + dat1$Study), np = 10000)
  perm <- cbind(perm, m$distribution[,1])
  roi <- c(ROI[i], fr$Estimate[2], 'Strength')
  roi <- as.data.frame(roi)
  roi <- t(roi)
  colnames(roi) <- c("connection", "effect", "type")
  roi.fr <- rbind(roi.fr, roi)
  fr <- data.frame(fr)
  fr <- fr[2,]
  colnames(fr) <- c("Estimate",   "Std..Error",     "t.value",     "P",  "outcome")
  Effects <- rbind(Effects, fr)
}
perm <- perm[,-1]
pa <- compute_minP(perm, alternative = TRUE)
paf <- data.frame(p=pa$main$pvalue, name=test)
p.perm <- rbind(p.perm,paf)

paf <- subset(paf, paf$p < 0.05)
PAF <- rbind(PAF, paf)



p.perm$permute <- p.perm$p
bu <- Effects
Effects$lci <- Effects$Estimate - 1.96*Effects$Std..Error
Effects$uci <- Effects$Estimate + 1.96*Effects$Std..Error
Effects <- merge(Effects, p.perm, by.x="outcome", by.y="name")
PLOT1 <- Effects

Effects$est <- paste0(round(Effects$Estimate,2), " [", round(Effects$lci,2), ",", round(Effects$uci,2), "]")
Effects <- Effects[,c("outcome", "est", "p", "permute")]
Effects$est <- gsub("0 [0,0]", "0.00 [0.00,0.00]", Effects$est)
Effects <- Effects %>% separate(outcome, c("measure", "region"))
Effects <- Effects %>% pivot_wider(id_cols = "region", names_from="measure", values_from=c(est, p, permute))

write.csv(Effects, "graph_est.csv")

#Diffs in global efficiency
summary(lm(dat1$global_ef ~ dat1$Diagnosis + dat1$Age + dat1$motion + dat1$Study))$coefficients


PLOT1 <- PLOT1[85:252,]
PLOT1$est <- paste0(round(PLOT1$Estimate,2))dat
PLOT1$sig <- ifelse(PLOT1$P < 0.05,"< 0.05","> 0.05")
PLOT1$sig <- ifelse(PLOT1$permute < 0.05,"< 0.05, corrected",PLOT1$sig)
PLOT1 <- PLOT1 %>% separate(outcome, c("measure", "region"))
PLOTe <- PLOT1[1:84,]
PLOTs <- PLOT1[85:168,]


p <- ggplot(PLOTs, aes(x=region, y=Estimate)) + geom_point(aes(colour=sig)) + geom_errorbar(aes(ymin=lci, ymax=uci, colour=sig), width=0.2) +  coord_flip()


cols <- c("#F8766D","#619CFF")
names(cols) <- levels(PLOTe$sig)
colScale <- scale_colour_manual(name = "sig",values = cols)

p <- ggplot(PLOTe, aes(x=region, y=Estimate)) + geom_point(aes(colour=sig)) + geom_errorbar(aes(ymin=lci, ymax=uci, colour=sig), width=0.2) +  coord_flip() 
p + colScale


#Read in results of NBS analyses to calculate component connectivity for each participants
sigs <- read.csv("tracts_stronger_HC.csv")

sig.t <- sigs %>% select(conn)
sig.t <- as.data.frame(sig.t)
dat$ave.comp1sum <- dat %>% select(sig.t$conn) %>% rowSums()
dat$ave.comp1 <- dat %>% select(sig.t$conn) %>% rowMeans()


SIG <- as.character(sig.t$conn)

sigs <- read.csv("tracts_stronger_AN.csv")sig.t <- sigs %>% select(conn)
sig.t <- as.data.frame(sig.t)
dat$ave.comp2 <- dat %>% select(sig.t$conn) %>% rowMeans()
dat$ave.comp2sum <- dat %>% select(sig.t$conn) %>% rowSums()
dat$rel <- dat$ave.comp2/dat$ave.comp1
SIG1 <- as.character(sig.t$conn)
SIG <- c(SIG, SIG1)

#descriptive stats per group
vars <- colnames(dat1)[c(2:272)]
tab1 <- dat1 %>%  group_by(Study, Dx) %>% summarise_at(c(vars), c(mean), na.rm=TRUE)
tab2 <- dat1 %>%  group_by(Study, Dx) %>% summarise_at(c(vars), c(sd), na.rm=TRUE)
tab1 <- data.frame(t(tab1))
tab2 <- data.frame(t(tab2))
tab3 <- tab1
nu <- sapply(tab1[3:270,], as.character)
nu1 <- sapply(tab2[3:270,], as.character)


for(i in 1:268){
  for(n in 1:10){
    x <- tryCatch({  
      nu[i,n] <- paste0(round(as.numeric(nu[i,n]),2), " (",round(as.numeric(nu1[i,n]),2), ")")}, error = function(e) {})
  }
}

nu <- data.frame(nu)
rownames(nu) <- rownames(tab1)[3:270]
write.csv(nu, "btwngrp_descriptives_bystudy.csv")



tab1 <- dat1 %>%  group_by(Dx) %>% summarise_at(c(vars), c(mean), na.rm=TRUE)
tab2 <- dat1 %>%  group_by(Dx) %>% summarise_at(c(vars), c(sd), na.rm=TRUE)
tab1 <- data.frame(t(tab1))
tab2 <- data.frame(t(tab2))

for(i in 2:268){
  for(n in 1:2){
    x <- tryCatch({  
      tab1[i,n] <- paste0(round(as.numeric(tab1[i,n]),2), " (",round(as.numeric(tab2[i,n]),2), ")")}, error = function(e) {})
  }
}



tab1 <- data.frame(tab1)

write.csv(tab1, "btwngrp_descriptives.csv")

#make plots for certain variables

for(i in c(1,2,10)){
  outcome <- vars[i]
  p <- ggplot(dat, aes_string(x="Study", y=outcome)) + geom_boxplot(aes(colour=Study)) + geom_jitter(aes(colour=Study)) + facet_grid(~Diagnosis)
  ggsave(paste0("dist_bystudy", outcome,".png"),p)
}


vars <- c("ave.comp1","ave.comp2","rel")
labs <- c("Subcortical component connectivity (mean streamline count)","Cortical component connectivity (mean streamline count)","Ratio of cortical to subcortical component connectivity (mean streamline count)")



for(i in 1:length(vars)){
  outcome <- vars[i]
  xlabel <- labs[i]
  p <- ggplot(dat, aes_string(x="Study", y=outcome)) + geom_boxplot(aes(colour=Study)) + geom_jitter(aes(colour=Study)) + facet_grid(~forcats::fct_rev(Diagnosis)) + theme(legend.position="none") + ylab(xlabel)
  ggsave(paste0("dist_bystudy", outcome,".png"),p)
}

#Subgroup analyses
MOD <- data.frame()


AN$sub.num <- ifelse(AN$Subtype=="Restrictive",-1,1)
vars <- c("ave.comp1","ave.comp2","rel", as.character(PAF$name))
for(i in 1:length(vars)){
  outcome <- vars[i]
  xlabel=labs[i]
  mod <- (summary(lm(AN[[outcome]] ~ AN$sub.num + AN$motion + AN$Study + AN$Age + AN$centile_scale))$coefficients)
  tab <- data.frame(mod[2,])
  tab <- data.frame(t(tab))
  tab$name <- xlabel
  MOD <- rbind(MOD, tab)
  print(ggplot(AN, aes_string(x="forcats::fct_rev(Subtype)", y=outcome)) + geom_boxplot(aes(colour=Subtype)) + geom_jitter(aes(colour=Subtype)) + ylab(xlabel)) + theme(legend.position="none")
}

MOD$lci <- MOD$Estimate - 1.96*MOD$Std..Error
MOD$uci <- MOD$Estimate + 1.96*MOD$Std..Error

MOD$est <- paste0(round(MOD$Estimate,2), " [", round(MOD$lci,2), ",", round(MOD$uci,2), "]")
MOD <- MOD[,c("name", "est", "Pr...t..")]

write.csv(MOD, "explore_assoc_subtype.csv")



#Associations of structural connectivity metrics with BMI
means <- c("ave.comp1", "ave.comp2", "rel", as.character(PAF$name))
AN[means] <- sapply(AN[means], scale)

labs <- c(rep("Subcortical component connectivity (mean streamline count)",2),rep("Cortical component connectivity (mean streamline count)",2), rep("Ratio of cortical to subcortical component connectivity (mean streamline count)",2), rep("Strength of left parahippocampal gyrus",2),
          rep("Strength of left cerebellem",2), rep("Strength of left hippocampus",2), rep("Strength of right lateral occipital gyrus",2), rep("Strength of right cerebellem",2))


MOD.AN <- data.frame()
PLOT <- data.frame()
labs <- unique(labs)
for(i in 1:length(means)){
  xlabel <- labs[i]
  var <- means[i]  
  mod <- summary(lm(AN$centile_scale ~ AN$Age.x + AN[[var]] + AN$motion.x + AN$Study.x))$coefficients
  fr <- as.data.frame(mod)
  fr$outcome <- "BMI"
  fr$pred <- var
  PLOT <- rbind(PLOT,(fr[3,]))
  if(fr$`Pr(>|t|)`[3] < 0.05){
    MOD.AN <- rbind(MOD.AN, fr)
  }
  p <- ggplot(AN, aes_string(x=var, y="centile_scale")) + geom_point(aes_string(colour="name"), show.legend = FALSE) + geom_smooth(method="lm", aes_string(colour="name"), show.legend = FALSE) + scale_colour_viridis(discrete=TRUE) + scale_fill_viridis(discrete=TRUE) + ylab("BMI scaled score") + xlab(xlabel) + theme_minimal()
  ggsave(paste0("explore_assoc_BMI_",var,".png"),p)
  p <- ggplot(AN, aes_string(x=var, y="Illness.Duration.x")) + geom_point(aes_string(colour="name"), show.legend = FALSE) + geom_smooth(method="lm", aes_string(colour="name"), show.legend = FALSE) + scale_colour_viridis(discrete=TRUE) + scale_fill_viridis(discrete=TRUE) + ylab("Illness Duration (months)") + xlab(xlabel) + theme_minimal()
  ggsave(paste0("explore_assoc_IllDur_",var,".png"),p)
  mod <- summary(lm(AN$Illness.Duration.x ~ AN[[var]] + AN$centile_scale  + AN$motion.x + AN$Study.x))$coefficients
  fr <- as.data.frame(mod)
  fr$outcome <- "Illness Duration"
  fr$pred <- var
  PLOT <- rbind(PLOT,(fr[2,]))
  if(fr$`Pr(>|t|)`[2] < 0.05){
    MOD.AN <- rbind(MOD.AN, fr)
  }
}

PLOT$ucl <- PLOT$Estimate + 1.96*PLOT$`Std. Error`
PLOT$lcl <- PLOT$Estimate - 1.96*PLOT$`Std. Error`

PLOT$est <- paste0(round(PLOT$Estimate,2), " [", round(PLOT$lcl,2), ",", round(PLOT$ucl,2), "]")
PLOT <- PLOT[,c("outcome","pred", "est", "Pr(>|t|)")]
PLOT <- PLOT %>% pivot_wider(id=pred, names_from=outcome, values_from=3:4)
write.csv(PLOT, "explore_assoc.csv")



#Visualise
library(brainconn)
atlas <- load("atlas.rda")
conmat <- read_excel("/Volumes/EDRU/Individual Folders/Caitlin/DTI/Paper1/NBS_fullres_tstat.xlsx", col_names = FALSE)

atlas <- load("atlas.rda")
atlas.own$ROI.Name <- atlas.own$name

atlas.own$x.mni <- as.integer(atlas.own$x.mni)
atlas.own$y.mni <- as.integer(atlas.own$y.mni)
atlas.own$z.mni <- as.integer(atlas.own$z.mni)



sig <- c("rTHAL","rHIPP", "lENT","lPARA","lCer","lTHAL","lCAUD","lPARH","rPARH", "rCer", "lHIPP", "rAMYG")
atlas.own$ROI.Name[!(atlas.own$name %in% sig)] <- ""
atlas.own$col[(atlas.own$name %in% sig)] <- "pink"
atlas.own$col[atlas.own$name=="lCer"] <- "green"
atlas.own$col[atlas.own$name=="rCer"] <- "green"
atlas.own$col[atlas.own$name=="lHIPP"] <- "green"


conmat[conmat < 3.1] <- 0
conmat[conmat > 3.1] <- 1

col <- subset(atlas.own, atlas.own$name %in% sig)
col <- col$col


atlas.own$node.size <- 8
atlas.own$node.size[(atlas.own$name %in% sig)] <- 10

non <- subset(atlas.own, !(atlas.own$name %in% sig))
num <- non$index
for(i in 1:length(num)){
  n <- num[i]
  conmat[,n] <- 0
  conmat[n,] <- 0
}



p <-brainconn(atlas ="atlas.own", conmat=conmat, edge.width = 1, node.size = 5, edge.color="cyan", node.color="cyan", labels=T,label.size = 5, all.nodes = F, view=c("top"))
ggsave("moreHC_static1.png",p)
p <-brainconn(atlas ="atlas.own", conmat=conmat, edge.width = 1, node.size = 5, edge.color="cyan", node.color="cyan", labels=T,label.size = 5, all.nodes = F, view=c("front"))
ggsave("moreHC_static2.png",p)
p <-brainconn(atlas ="atlas.own", conmat=conmat, edge.width = 1, node.size = 5, edge.color="cyan", node.color="cyan", labels=T,label.size = 5, all.nodes = F, view=c("left"))
ggsave("moreHC_static3.png",p)

conmat <- read_excel("NBS_fullres_tstat.xlsx", col_names = FALSE)

atlas <- load("atlas.rda")
atlas.own$ROI.Name <- atlas.own$name

atlas.own$x.mni <- as.integer(atlas.own$x.mni)
atlas.own$y.mni <- as.integer(atlas.own$y.mni)
atlas.own$z.mni <- as.integer(atlas.own$z.mni)

atlas.own$col <- "pink"

sig <- c("lpORB", "lpOPER", "lpreC", "rrACC", "lLOF", "lSTG", "rFP", "rPUT")
atlas.own$ROI.Name[!(atlas.own$name %in% sig)] <- ""

col <- c(rep("cyan",2), "pink", rep("cyan",5))


atlas.own$node.size <- 8
atlas.own$node.size[(atlas.own$name %in% sig)] <- 10


conmat[conmat > -3.1] <- 0
conmat[conmat < -3.1] <- 1

non <- subset(atlas.own, !(atlas.own$name %in% sig))
num <- non$index
for(i in 1:length(num)){
  n <- num[i]
  conmat[,n] <- 0
  conmat[n,] <- 0
}



p <-brainconn(atlas ="atlas.own", conmat=conmat, edge.width = 1, node.size = 5, edge.color="pink", node.color="pink", labels=T,label.size = 5, all.nodes = F, view=c("top"))
ggsave("ANmoreHC_static1.png",p)
p <-brainconn(atlas ="atlas.own", conmat=conmat, edge.width = 1, node.size = 5, edge.color="pink", node.color="pink", labels=T,label.size = 5, all.nodes = F, view=c("front"))
ggsave("ANmoreHC_static2.png",p)
p <-brainconn(atlas ="atlas.own", conmat=conmat, edge.width = 1, node.size = 5, edge.color="pink", node.color="pink", labels=T,label.size = 5, all.nodes = F, view=c("left"))
ggsave("ANmoreHC_static3.png",p)
