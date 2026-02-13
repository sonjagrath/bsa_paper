data_CCRT <- read_excel("phenotype_RIL_IL.xlsx", sheet = "CCRT")

levels <- c("BKK5", "BKK6", "BKK10", "BKK12", 
            "BKK13", "BKK16", "BKK17", "BKK18", 
            "KATH14", "KATH19", "KATH23", 
            "RIL7", "RIL14", "RIL15", "RIL20", "RIL22", "RIL23", "RIL25", "RIL30", 
            "RIL41", "RIL47", "RIL50", "RIL57", "RIL58", "RIL80", "RIL81", "RIL93")
data_CCRT$Line <- factor(data_CCRT$Line, levels = levels, labels = levels)
data_CCRT$Sex <- as.factor(data_CCRT$Sex)

data_CCRT[which(data_CCRT$Line == "BKK5"),]$Root <- "FastBKK"
data_CCRT[which(data_CCRT$Line == "BKK6"),]$Root <- "FastBKK"
data_CCRT[which(data_CCRT$Line == "BKK10"),]$Root <- "FastBKK"
data_CCRT[which(data_CCRT$Line == "BKK12"),]$Root <- "FastBKK"
data_CCRT[which(data_CCRT$Line == "BKK13"),]$Root <- "SlowBKK"
data_CCRT[which(data_CCRT$Line == "BKK16"),]$Root <- "SlowBKK"
data_CCRT[which(data_CCRT$Line == "BKK17"),]$Root <- "SlowBKK"
data_CCRT[which(data_CCRT$Line == "BKK18"),]$Root <- "SlowBKK"
data_CCRT[which(data_CCRT$Line == "KATH14"),]$Root <- "KATH"
data_CCRT[which(data_CCRT$Line == "KATH19"),]$Root <- "KATH"
data_CCRT[which(data_CCRT$Line == "KATH23"),]$Root <- "KATH"

data_CCRT$Root <- as.factor(data_CCRT$Root)

anova <- aov(Time~Sex*Line, data = data_CCRT)
summary(anova)
tukey <- TukeyHSD(anova) 
print(multcompLetters4(anova, tukey))

plot_CCRT <- ggboxplot(data_CCRT, "Line", "Time", fill = "Root", facet.by = "Sex") + scale_fill_viridis_d(begin = 0.2, end = 0.9, alpha = 0.8)
ggpar(plot_CCRT, x.text.angle = 90, ylab = "Chill coma recovery time (min)", yticks.by = 10, legend = "bottom", ggtheme = theme_classic(), legend.title = "Population")

female_CCRT <- data_CCRT[which(data_CCRT$Sex == "Female"),]
female_CCRT$Line <- factor(female_CCRT$Line, levels = levels, labels = levels)

male_CCRT <- data_CCRT[which(data_CCRT$Sex == "Male"),]
male_CCRT$Line <- factor(male_CCRT$Line, levels = levels, labels = levels)

female_crt <- female_CCRT %>%
#male_crt <- male_CCRT %>%
  group_by(Line) %>%
  summarize(sd(Time), mean(Time))

########### LTi50 ###########

#LTi50 was calculated by following code provided by Dirk Metzler:

# import the excel sheet
library(readxl)
LTi50 <- read_excel("phenotype_RIL_IL.xlsx", sheet = "RIL_Mortlity") # for RILs
# LTi50 <- read_excel("phenotype_RIL_IL.xlsx", sheet = "IL_Mortlity") # for iso-female lines (separate due to difference in time points)
View(LTi50)

# choose Line and Sex

data <- LTi50[which(LTi50$Line == "RIL1" & LTi50$Sex == "Male"),]

# generate the model and test with ANOVA
mod <- glm(cbind(Mortality, Total-Mortality) ~ Time, data, family="quasibinomial")
summary(mod)
anova(mod, test="LRT")

# make predictions according to the data

newdata <- data.frame(Time=0:250/10)
predict(mod, newdata, type="response")

# plot

plot(data$Time, data$Mortality/data$Total, xlab = "Time of exposure (h)", ylab = "Proportional mortality")
lines(newdata$Time, predict(mod, newdata, type="response"))

# find the time point 50% mortality
-mod$coef[1]/mod$coef[2]

# the resulting numbers correspond to LTi50 (of each line and sex) and are used in later analysis

#######

data_LTi <- read_excel("phenotype_RIL_IL.xlsx", sheet = "LTi_CS")

data_LTi$Line <- factor(data_LTi$Line, levels = levels, labels = levels)
data_LTi$Sex <- as.factor(data_LTi$Sex)
data_LTi$Root <- as.factor(data_LTi$Root)

anova <- aov(LTi50~Sex*Line, data = data_LTi)
summary(anova)
anova <- aov(LTi50~Line, data = data_LTi)
tukey <- TukeyHSD(anova) 
print(multcompLetters4(anova, tukey))
anova <- aov(LTi50~Sex, data = data_LTi)
tukey <- TukeyHSD(anova) 
print(multcompLetters4(anova, tukey))

female_LTi <- data_LTi[which(data_LTi$Sex == "Female"),]
male_LTi <- data_LTi[which(data_LTi$Sex == "Male"),]

anova <- aov(LTi50~Line, data = female_LTi)
tukey <- TukeyHSD(anova) 
print(multcompLetters4(anova, tukey))

anova <- aov(LTi50~Line, data = male_LTi)
tukey <- TukeyHSD(anova) 
print(multcompLetters4(anova, tukey))

df <- data_LTi[, -4]
df_wide <- df %>%
  pivot_wider(names_from = Sex, values_from = LTi50)
t.test(df_wide$Female, df_wide$Male, paired = TRUE)

########### CS Mortality ###########

data_CS <- read_excel("phenotype_RIL_IL.xlsx", sheet = "CS")

data_CS$Line <- factor(data_CS$Line, levels = levels, labels = levels)
data_CS$Sex <- as.factor(data_CS$Sex)
data_CS$Root <- as.factor(data_CS$Root)

data_CS[which(data_CS$Line == "BKK5"),]$Root <- "FastBKK"
data_CS[which(data_CS$Line == "BKK6"),]$Root <- "FastBKK"
data_CS[which(data_CS$Line == "BKK10"),]$Root <- "FastBKK"
data_CS[which(data_CS$Line == "BKK12"),]$Root <- "FastBKK"
data_CS[which(data_CS$Line == "BKK13"),]$Root <- "SlowBKK"
data_CS[which(data_CS$Line == "BKK16"),]$Root <- "SlowBKK"
data_CS[which(data_CS$Line == "BKK17"),]$Root <- "SlowBKK"
data_CS[which(data_CS$Line == "BKK18"),]$Root <- "SlowBKK"
data_CS[which(data_CS$Line == "KATH14"),]$Root <- "KATH"
data_CS[which(data_CS$Line == "KATH19"),]$Root <- "KATH"
data_CS[which(data_CS$Line == "KATH23"),]$Root <- "KATH"


anova <- aov(Mortality~Sex*Line, data = data_CS)
summary(anova)
tukey <- TukeyHSD(anova) 
print(multcompLetters4(anova, tukey))

data_CS$Percentage <- data_CS$Mortality/10

plot_CS <- ggbarplot(data_CS, "Line", "Percentage", fill = "Root", facet.by = "Sex", merge = T, add = c("mean", "dotplot"), add.params = list(size = 0.5, fill = "Root")) + scale_fill_viridis_d(begin = 0.2, end = 0.9, alpha = 0.5)
ggpar(plot_CS, x.text.angle = 90, ylab = "Percent mortality upon 8h cold shock", yticks.by = 0.1, legend = "bottom", ggtheme = theme_classic(), legend.title = "Population")

female_CS <- data_CS[which(data_CS$Sex == "Female"),]
male_CS <- data_CS[which(data_CS$Sex == "Male"),]

anova <- aov(Mortality~Line, data = female_CS)
tukey <- TukeyHSD(anova) 
print(multcompLetters4(anova, tukey))

anova <- aov(Mortality~Line, data = male_CS)
tukey <- TukeyHSD(anova) 
print(multcompLetters4(anova, tukey))

#female <- female_CS %>%
male <- male_CS %>%
  group_by(Line) %>%
  summarize(sd(Mortality), mean(Mortality))

