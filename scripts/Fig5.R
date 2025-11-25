rm(list=ls())
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(moments)
library(hrbrthemes)
library(extrafont)
extrafont::loadfonts() 

d_all <- read.csv("~/data_ind.csv")
d_all$Group <- factor(d_all$Group)
levels(d_all$Group)
levels(d_all$Group)[9] <- "Eukaryotic Microalgae"
levels(d_all$Group)[c(6,8,15)] <- "Crustacea"

n_ind <- d_all %>% group_by(Species) %>% dplyr::summarise(n_ind=n())
d_all <- merge(d_all, n_ind, by="Species")
dd <- d_all %>% filter(n_ind>4) 

## Empirical species allometric exponents -------------------------------------------
dd2 <- dd[dd$Group %in% c("Mammals", "Birds", "Amphibians", "Reptiles", "Fishes", "Crustacea"),] 
id <- unique(dd2$Species)
beta <- count <- n <- gr <- M <- B <- r2 <- list()
for(i in 1:length(id)){
  dtemp <- dd2 %>% filter(Species==id[i] )
  beta[[i]] <- coef(lm( log(dtemp$MR_25) ~ log(dtemp$Wet_Mass) ))[2]
  n[[i]] <- unique(dtemp$n_ind)
  gr[[i]] <- unique(dtemp$Group)
  M[[i]] <- mean(dtemp$Wet_Mass)
  B[[i]] <- mean(dtemp$MR_25)
  r2[[i]] <- summary(lm(log( dtemp$MR_25 ) ~ log(dtemp$Wet_Mass) ))$r.sq
}

out <- data.frame(Group=unlist(gr), 
                  Species=id, 
                  n=unlist(n), 
                  Mass=unlist(M), 
                  MR=unlist(B),
                  beta=unlist(beta),  
                  R2= unlist(r2)) %>% arrange(Group)
out
xtable::xtable(out,  display = c("s","s","s", "s","g", "g","g", "g"))

## FIG 5 ------------------------------------------------------------
summary(lm(log10(dd$MR_25) ~ log10(dd$Wet_Mass)))

fig5a <- ggplot(dd2)+
  geom_point(aes(x=log10(Wet_Mass), y=log10(MR_25) , fill=Group), shape=21, size=1.5)+
  geom_smooth(aes(x=log10(Wet_Mass), y=log10(MR_25)), method="lm", se=F, color="black", size=0.5 )+
  theme_ipsum(axis_title_size = 14, axis_title_just = "m" ) +
  theme(panel.grid.major = element_line(color = "white", size = 0.25),
        panel.grid.minor = element_line(color = "white", size = 0.25),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5),
        legend.position = c(0.2, 0.75), legend.title = element_blank(),
        legend.text = element_text(size = 12) )+
  labs(x = "log Mass",
       y = "log Metabolic Rate") +
  geom_text(x = Inf, y = Inf,
            hjust = 1.3, vjust = 17.5,  
            label = expression(paste( beta, " = 0.87")),
            color = 'black', size=5, family = "Arial Narrow")+
  geom_text(x = Inf, y = Inf,
            hjust = 1.3, vjust = 18,  
            label = expression(paste( R^2, " = 0.9")),
            color = 'black', size=5, family = "Arial Narrow")+
  scale_fill_brewer(palette="Spectral")
fig5a

ggsave("~/fig5a.png", width=5, height=5)

library(plyr)
slopes <- ddply(dd2, c("Group"), summarise,
                slope = round(coef(lm(log10(MR_25)~log10(Wet_Mass)))[2], 2),
                r2    = round(summary(lm(log10(MR_25)~log10(Wet_Mass)))$r.squared, 2))

fig5b <- ggplot(dd2)+
  facet_wrap(~Group, scales = "free", ncol=2)+
  geom_smooth(aes(x=log10(Carbon_Mass), y=log10(MR_25), group=Species), color="grey40",
              method="lm", se=F,size=0.35 )+
  geom_smooth(aes(x=log10(Carbon_Mass), y=log10(MR_25)), method="lm", se=F, color="black", size=0.5 )+
  theme_ipsum(axis_title_size = 14, axis_title_just = "m" ) +
  theme(panel.grid.major = element_line(color = "white", size = 0.25),
        panel.grid.minor = element_line(color = "white", size = 0.25),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5),
        panel.spacing = unit(0, "cm"),
        legend.position = "none")+
  scale_color_brewer(palette="Spectral")+
  labs(x = "log Mass",
       y = "log Metabolic Rate") +
  geom_text(data = slopes,
            aes(label = paste("beta ==", slope)),
            parse = TRUE,
            x = Inf, y = Inf,
            hjust = 1.3, vjust = 11,  
            size = 3, color = "black", family = "Arial Narrow") +
  geom_text(data = slopes,
            aes(label = paste("R^2 ==", r2)),
            parse = TRUE,
            x = Inf, y = Inf,
            hjust = 1.3, vjust = 12,  
            size = 3, color = "black", family = "Arial Narrow")
fig5b

ggsave("~/fig5b.png", width=5, height=7)

slopes2 <- ddply(dd2, c("Species", "Group"), summarise,
                beta = round(coef(lm(log10(MR_25)~log10(Wet_Mass)))[2], 2),
                r2    = round(summary(lm(log10(MR_25)~log10(Wet_Mass)))$r.squared, 2),
                n = unique(n_ind))

beta.sp <- slopes2 %>% dplyr::filter( r2 > 0.9)
cols <- brewer.pal(n = 6, name = "Spectral")

fig5c <- ggplot(beta.sp, aes(x = Group, y = beta)) +
  annotate("rect",
           xmin = -Inf, xmax = Inf,
           ymin = 0.5, ymax = 1,
           fill = "gray75", alpha = 0.5) +
  stat_boxplot(geom = "errorbar", width = 0.4, color = "steelblue") +
  geom_boxplot(aes(group = Group), fill = "lightblue", color = NA, outliers = F, width = 0.5) +
  stat_summary(fun = median, geom = "crossbar", width = 0.5, fatten = 1, color = "black") +
  geom_jitter(aes(color = Group), shape = 18, size = 3, stroke = 1, width = 0.15) +
  labs(y = expression(paste("Empirical scaling exponents ", tilde(beta)[S])), x = "")+
  scale_color_manual(values=c("#D53E4F", "#FEE08B", "#E6F598", "#99D594", "#3288BD") )+
 # scale_color_brewer(palette="Spectral")+
  theme_ipsum(axis_title_size = 14, axis_title_just = "m" ) +
  theme(
    panel.grid.major = element_line(color = "white", size = 0.25),
    panel.grid.minor = element_line(color = "white", size = 0.25),
    panel.border = element_rect(colour = "black", fill=NA, linewidth=0.25),
    panel.spacing = unit(0, "cm"),
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1, size= 12, color = "black")) +
  coord_cartesian(clip = "off")
fig5c

ggsave("~/fig5c.png", width=5, height=4)








