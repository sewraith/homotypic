setwd("/Users/sewraith/Box/Research/Homotypic")
getwd()

#packages
require(dplyr)
require(data.table)
require(ggplot2)
require(tidyr)
library(plotly)

#getting data set up
H1N1_cumulative <- read.csv(file="/Users/sewraith/Box/Research/Homotypic/H1N1_cumulative.csv", header=TRUE, sep=",")
H1N1_agestrat_under5 <- read.csv(file="/Users/sewraith/Box/Research/Homotypic/H1N1_agestrat_under5.csv", header=TRUE, sep=",")
H1N1_agestrat_over5 <- read.csv(file="/Users/sewraith/Box/Research/Homotypic/H1N1_agestrat_over5.csv", header=TRUE, sep=",")
H1N1_agestrat_under5lim <- read.csv(file="/Users/sewraith/Box/Research/Homotypic/H1N1_agestrat_under5lim.csv", header=TRUE, sep=",")
H1N1_agestrat_over5lim <- read.csv(file="/Users/sewraith/Box/Research/Homotypic/H1N1_agestrat_over5lim.csv", header=TRUE, sep=",")

H3N2_cumulative <- read.csv(file="/Users/sewraith/Box/Research/Homotypic/H3N2_cumulative.csv", header=TRUE, sep=",")
H3N2_agestrat_under5 <- read.csv(file="/Users/sewraith/Box/Research/Homotypic/H3N2_agestrat_under5.csv", header=TRUE, sep=",")
H3N2_agestrat_over5 <- read.csv(file="/Users/sewraith/Box/Research/Homotypic/H3N2_agestrat_over5.csv", header=TRUE, sep=",")
H3N2_agestrat_under5lim <- read.csv(file="/Users/sewraith/Box/Research/Homotypic/H3N2_agestrat_under5lim.csv", header=TRUE, sep=",")
H3N2_agestrat_over5lim <- read.csv(file="/Users/sewraith/Box/Research/Homotypic/H3N2_agestrat_over5lim.csv", header=TRUE, sep=",")

IB_cumulative <- read.csv(file="/Users/sewraith/Box/Research/Homotypic/IB_cumulative.csv", header=TRUE, sep=",")
IB_agestrat_under5 <- read.csv(file="/Users/sewraith/Box/Research/Homotypic/IB_agestrat_under5.csv", header=TRUE, sep=",")
IB_agestrat_over5 <- read.csv(file="/Users/sewraith/Box/Research/Homotypic/IB_agestrat_over5.csv", header=TRUE, sep=",")

IB_yamagata <- read.csv(file="/Users/sewraith/Box/Research/Homotypic/IB_yamagata.csv", header=TRUE, sep=",")
IB_yamagata_agestrat_under5 <- read.csv(file="/Users/sewraith/Box/Research/Homotypic/IB_yamagata_agestrat_under5.csv", header=TRUE, sep=",")
IB_yamagata_agestrat_over5 <- read.csv(file="/Users/sewraith/Box/Research/Homotypic/IB_yamagata_agestrat_over5.csv", header=TRUE, sep=",")

IB_victoria <- read.csv(file="/Users/sewraith/Box/Research/Homotypic/IB_victoria.csv", header=TRUE, sep=",")
IB_victoria_agestrat_under5 <- read.csv(file="/Users/sewraith/Box/Research/Homotypic/IB_victoria_agestrat_under5.csv", header=TRUE, sep=",")
IB_victoria_agestrat_over5 <- read.csv(file="/Users/sewraith/Box/Research/Homotypic/IB_victoria_agestrat_over5.csv", header=TRUE, sep=",")
IB_victoria_agestrat_under5lim <- read.csv(file="/Users/sewraith/Box/Research/Homotypic/IB_victoria_agestrat_under5lim.csv", header=TRUE, sep=",")
IB_victoria_agestrat_over5lim <- read.csv(file="/Users/sewraith/Box/Research/Homotypic/IB_victoria_agestrat_over5lim.csv", header=TRUE, sep=",")


#H1N1 figures
h1n1all <- H1N1_cumulative
h1n1all$Season <- as.character(h1n1all$Season)
h1n1all$Season <- factor(h1n1all$Season, levels = h1n1all$Season)

figure1 <- ggplot(h1n1all, aes(x = Season, y = OR, ymin = l, ymax = u, color = Season)) +
  geom_pointrange(position = position_dodge(width = 0.50), size=.75 ) +
  scale_y_log10(breaks = c(0,0.2,1,5,20,80)) +
  geom_hline(yintercept = 1) +
  theme_classic() +
  coord_flip() 

#use this one
figure1 <- ggplot(h1n1all, aes(x = OR, y = Season, xmin = l, xmax = u, color = year)) +
  geom_pointrange(position = position_dodge(width = 0.50), size=.75 ) +
  scale_x_log10(breaks = c(-0.2,0,0.01,0.2,1,5,20,80)) +
  scale_y_discrete(limits = rev) +
  geom_vline(xintercept = 1) +
  theme_classic() 

figure1 + theme(text = element_text(size=20), axis.text.x=element_blank(), axis.ticks.x=element_blank())

#use this one
figure1 + theme(text = element_text(size=30), legend.position="none")

h1n1stratlimu <- H1N1_agestrat_under5lim
h1n1stratlimu$Season <- factor(h1n1stratlimu$Season, levels = h1n1stratlimu$Season)

figure1 <- ggplot(h1n1stratlimu, aes(x = OR, y = Season, xmin = l, xmax = u, color = year)) +
  geom_pointrange(position = position_dodge(width = 0.50), size=.75 ) +
  scale_x_log10(breaks = c(-0.2,0,0.01,0.2,1,5,20,80)) +
  geom_vline(xintercept = 1) +
  theme_classic() 

figure2 <- ggplot(h1n1stratlimu, aes(x = Season, y = OR, ymin = l, ymax = u, color = year)) +
  geom_pointrange(position = position_dodge(width = 0.50), size=.75) +
  geom_hline(yintercept = 1) +
  scale_y_log10(breaks = c(0,0.01,0.2,1,5,20,80)) +
  scale_x_discrete(limits = rev) +
  theme_classic() +
  coord_flip()

figure2 + theme(text = element_text(size=30), legend.position="none")

H1N1_agestrat_under5 <- H1N1_agestrat_under5 %>% mutate(age="4 and under")
H1N1_agestrat_under5$Season <- factor(H1N1_agestrat_under5$Season, levels = H1N1_agestrat_under5$Season)
H1N1_agestrat_over5 <- H1N1_agestrat_over5 %>% mutate(age="5 and over")
H1N1_agestrat_over5$Season <- factor(H1N1_agestrat_over5$Season, levels = H1N1_agestrat_over5$Season)


H1N1_agestrat <- rbind(H1N1_agestrat_under5, H1N1_agestrat_over5)


h1n1age <- H1N1_agestrat

#use this one
figure2 <- ggplot(h1n1age, aes(x = Season, y = OR, ymin = l, ymax = u, color = age, shape=age)) +
  geom_pointrange(position = position_dodge(width = 0.50), size=.75) +
  geom_hline(yintercept = 1) +
  scale_y_log10(breaks = c(0,0.2,1,5,20,80)) +
  theme_classic() +
  coord_flip()

figure2 <- ggplot(h1n1age, aes(x = Season, y = OR, ymin = l, ymax = u, color = year, shape=age)) +
  geom_pointrange(position = position_dodge(width = 0.50), size=.75) +
  geom_hline(yintercept = 1) +
  scale_y_log10(breaks = c(0,0.01,0.2,1,5,20,80)) +
  scale_x_discrete(limits = rev) +
  theme_classic() +
  coord_flip()

figure2 <- ggplot(h1n1age, aes(x = OR, y = Season, xmin = l, xmax = u, color = Season)) +
  geom_pointrange(position = position_dodge(width = 0.50), size=.75 ) +
  scale_x_log10(breaks = c(-0.2,0,0.2,1,5,20,80)) +
  geom_vline(xintercept = 1) +
  theme_classic() 

#use this one
figure2 + theme(text = element_text(size=30), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

figure2 + theme(text = element_text(size=30))

#H3N2 figures
h3n2all <- H3N2_cumulative
h3n2all$Season <- factor(h3n2all$Season, levels = h3n2all$Season)

figure3 <- ggplot(h3n2all, aes(x = Season, y = OR, ymin = l, ymax = u, color = Season)) +
  geom_pointrange(position = position_dodge(width = 0.50), size=.75 ) +
  scale_y_log10(breaks = c(0,0.2,1,5,20,80)) +
  geom_hline(yintercept = 1) +
  theme_classic()+
  coord_flip()

#use this one
figure3 <- ggplot(h3n2all, aes(x = OR, y = Season, xmin = l, xmax = u, color = year)) +
  geom_pointrange(position = position_dodge(width = 0.50), size=.75 ) +
  scale_x_log10(breaks = c(-0.2,0,0.01,0.2,1,5,20,80)) +
  scale_y_discrete(limits = rev) +
  geom_vline(xintercept = 1) +
  theme_classic() 

figure3 + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

figure3 + theme(text = element_text(size=30), legend.position="none")

H3N2_agestrat_under5 <- H3N2_agestrat_under5 %>% mutate(age="4 and under")
H3N2_agestrat_under5$Season <- factor(H3N2_agestrat_under5$Season, levels = H3N2_agestrat_under5$Season)
H3N2_agestrat_over5 <- H3N2_agestrat_over5 %>% mutate(age="5 and over")
H3N2_agestrat_over5$Season <- factor(H3N2_agestrat_over5$Season, levels = H3N2_agestrat_over5$Season)


H3N2_agestrat <- rbind(H3N2_agestrat_under5, H3N2_agestrat_over5)

h3n2age <- H3N2_agestrat
figure4 <- ggplot(h3n2age, aes(x = Season, y = OR, ymin = l, ymax = u, color = age, shape=age)) +
  geom_pointrange(position = position_dodge(width = 0.50), size=.75) +
  geom_hline(yintercept = 1) +
  scale_y_log10(breaks = c(0, 0.2,1,5,20,80)) +
  theme_classic()+
  coord_flip()

figure4 <- ggplot(h3n2age, aes(x = Season, y = OR, ymin = l, ymax = u, color = year, shape=age)) +
  geom_pointrange(position = position_dodge(width = 0.50), size=.75) +
  geom_hline(yintercept = 1) +
  scale_y_log10(breaks = c(0,0.01,0.2,1,5,20,80)) +
  scale_x_discrete(limits = rev) +
  theme_classic() +
  coord_flip()

figure4 + theme(text = element_text(size=30), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

figure4 + theme(text = element_text(size=20))

#IB figures
iball <- IB_cumulative

figure5 <- ggplot(iball, aes(x = Season, y = OR, ymin = l, ymax = u, color = Season)) +
  geom_pointrange(position = position_dodge(width = 0.50), size=.75 ) +
  scale_y_log10(breaks = c(0,0.2,1,5,20,80)) +
  geom_hline(yintercept = 1) +
  theme_classic()+
  coord_flip()

figure5 + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

figure5 + theme(text = element_text(size=20), legend.position="none")

IB_agestrat_under5 <- IB_agestrat_under5 %>% mutate(age="4 and under")
IB_agestrat_over5 <- IB_agestrat_over5 %>% mutate(age="5 and over")

IB_agestrat <- rbind(IB_agestrat_under5, IB_agestrat_over5)

ibage <- IB_agestrat
figure6 <- ggplot(ibage, aes(x = Season, y = OR, ymin = l, ymax = u, color = age, shape=age)) +
  geom_pointrange(position = position_dodge(width = 0.50), size=.75) +
  geom_hline(yintercept = 1) +
  scale_y_log10(breaks = c(0, 0.2,1,5,20,80)) +
  theme_classic()+
  coord_flip()

figure6 + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

figure6 + theme(text = element_text(size=20))

ibyam <- IB_yamagata
ibyam$Season <- factor(ibyam$Season, levels = ibyam$Season)

figure7 <- ggplot(ibyam, aes(x = Season, y = OR, ymin = l, ymax = u, color = Season)) +
  geom_pointrange(position = position_dodge(width = 0.50), size=.75 ) +
  scale_y_log10(breaks = c(0,0.2,1,5,20,80)) +
  geom_hline(yintercept = 1) +
  theme_classic()+
  coord_flip()

figure7 <- ggplot(ibyam, aes(x = OR, y = Season, xmin = l, xmax = u, color = year)) +
  geom_pointrange(position = position_dodge(width = 0.50), size=.75 ) +
  scale_x_log10(breaks = c(-0.2,0,0.01,0.2,1,5,20,80)) +
  scale_y_discrete(limits = rev) +
  geom_vline(xintercept = 1) +
  theme_classic() 

figure7 + theme(text = element_text(size=30), legend.position="none")


IB_yamagata_agestrat_under5 <- IB_yamagata_agestrat_under5 %>% mutate(age="4 and under")
IB_yamagata_agestrat_under5$Season <- factor(IB_yamagata_agestrat_under5$Season, levels = IB_yamagata_agestrat_under5$Season)
IB_yamagata_agestrat_over5 <- IB_yamagata_agestrat_over5 %>% mutate(age="5 and over")
IB_yamagata_agestrat_over5$Season <- factor(IB_yamagata_agestrat_over5$Season, levels = IB_yamagata_agestrat_over5$Season)

IByam_agestrat <- rbind(IB_yamagata_agestrat_under5, IB_yamagata_agestrat_over5)

ibyamage <- IByam_agestrat
figure8 <- ggplot(ibyamage, aes(x = Season, y = OR, ymin = l, ymax = u, color = age, shape=age)) +
  geom_pointrange(position = position_dodge(width = 0.50), size=.75) +
  geom_hline(yintercept = 1) +
  scale_y_log10(breaks = c(0, 0.2,1,5,20,80)) +
  theme_classic()+
  coord_flip()

figure8 <- ggplot(ibyamage, aes(x = Season, y = OR, ymin = l, ymax = u, color = year, shape=age)) +
  geom_pointrange(position = position_dodge(width = 0.50), size=.75) +
  geom_hline(yintercept = 1) +
  scale_y_log10(breaks = c(0,0.01,0.2,1,5,20,80)) +
  scale_x_discrete(limits = rev) +
  theme_classic() +
  coord_flip()


figure8 + theme(text = element_text(size=30),axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))


ibvic <- IB_victoria
ibvic$Season <- factor(ibvic$Season, levels = ibvic$Season)

figure8 <- ggplot(ibvic, aes(x = Season, y = OR, ymin = l, ymax = u, color = Season)) +
  geom_pointrange(position = position_dodge(width = 0.50), size=.75 ) +
  scale_y_log10(breaks = c(0,0.2,1,5,20,80)) +
  geom_hline(yintercept = 1) +
  theme_classic()+
  coord_flip()

figure8 <- ggplot(ibvic, aes(x = OR, y = Season, xmin = l, xmax = u, color = year)) +
  geom_pointrange(position = position_dodge(width = 0.50), size=.75 ) +
  scale_x_log10(breaks = c(-0.2,0,0.01,0.2,1,5,20,80)) +
  scale_y_discrete(limits = rev) +
  geom_vline(xintercept = 1) +
  theme_classic() 


figure8 + theme(text = element_text(size=30), legend.position="none")


IB_victoria_agestrat_under5 <- IB_victoria_agestrat_under5 %>% mutate(age="4 and under")
IB_victoria_agestrat_under5$Season <- factor(IB_victoria_agestrat_under5$Season, levels = IB_victoria_agestrat_under5$Season)
IB_victoria_agestrat_over5 <- IB_victoria_agestrat_over5 %>% mutate(age="5 and over")
IB_victoria_agestrat_over5$Season <- factor(IB_victoria_agestrat_over5$Season, levels = IB_victoria_agestrat_over5$Season)

IBvic_agestrat <- rbind(IB_victoria_agestrat_under5, IB_victoria_agestrat_over5)

ibvicage <- IBvic_agestrat
figure8 <- ggplot(ibvicage, aes(x = Season, y = OR, ymin = l, ymax = u, color = age, shape=age)) +
  geom_pointrange(position = position_dodge(width = 0.50), size=.75) +
  geom_hline(yintercept = 1) +
  scale_y_log10(breaks = c(0, 0.2,1,5,20,80)) +
  theme_classic()+
  coord_flip()

figure8 <- ggplot(ibvicage, aes(x = Season, y = OR, ymin = l, ymax = u, color = year, shape=age)) +
  geom_pointrange(position = position_dodge(width = 0.50), size=.75) +
  geom_hline(yintercept = 1) +
  scale_y_log10(breaks = c(0,0.01,0.2,1,5,20,80)) +
  scale_x_discrete(limits = rev) +
  theme_classic() +
  coord_flip()

figure8 + theme(text = element_text(size=30),axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))


