load("~/analytical.RData")
library(ggplot2)
library(RColorBrewer)
library(hrbrthemes)
library(extrafont)
extrafont::loadfonts() 


### Plots --------------------------------------------------------

slopes <- d %>%
  group_by(gamma, phi) %>%
  do({
    fit <- lm(log(B) ~ log(mass), data = .)
    coef_df <- tidy(fit) %>% filter(term == "log(mass)")
    r2 <- glance(fit)$r.squared
    coef_df$r2 <- r2
    coef_df
  })

## FIG S1 ---------------------------------------------------------------------

ggplot(slopes)+
  geom_histogram(aes(x=r2), color="grey80", fill="grey10")+
  theme_ipsum(axis_title_size = 18, axis_title_just = "m" ) +
  theme(panel.grid.major = element_line(color = "white", linewidth = 0.25),
        panel.grid.minor = element_line(color = "white", linewidth = 0.25),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5))+
  labs(x = expression(paste("Model fit, ", R^2) )) 

ggsave("~/FigS1.png", height = 5, width = 5)

## FIG 2 ---------------------------------------------------------------------

fig2a <- ggplot(slopes)+
  geom_tile(aes(x=gamma, y=phi, fill=estimate))+
  geom_contour(aes(x = gamma, y = phi, z = estimate),
               breaks = 0.75, color = "black", linewidth = 0.75, linetype="dashed") +
  labs(y = expression(paste("Baseline fluctuation level ", phi)),
       x = expression(paste("Stocastic inefficiency ", gamma)),
       fill = expression(beta) ) +
  theme_ipsum(axis_title_size = 22, axis_title_just = "m" ) +
  scale_fill_distiller(palette="Spectral")+
  theme(panel.grid.major = element_line(color = "white", linewidth = 0.25),
        panel.grid.minor = element_line(color = "white", linewidth = 0.25),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5))
fig2a

fig2b <- ggplot(d)+geom_tile(aes(x=gamma, y=phi, fill=beta))+
  facet_wrap(~animal)+
  labs(
    y = expression(paste("Baseline fluctuation level ", phi)),
    x = expression(paste("Stocastic inefficiency ", gamma)),
    fill = expression(beta[S])
  ) +
  theme_ipsum(axis_title_size = 22, axis_title_just = "m") +
  theme(panel.grid.major = element_line(color = "white", linewidth = 0.25),
        panel.grid.minor = element_line(color = "white", linewidth = 0.25),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5))+
  scale_fill_distiller(palette="Spectral")

fig2b

library(cowplot)
setwd("~/figures/")

grid <- plot_grid( fig2a, fig2b,
                   labels = "auto", 
                   nrow=1, align="h", label_size = 24, 
                   rel_widths = c(1, 1.3),
                   label_fontface = "plain")

save_plot("Fig2.png", grid,
          ncol = 2, 
          nrow = 1.5, 
          base_aspect_ratio = 1.9)

## FIG S2 ---------------------------------------------------------------------

figS2 <- ggplot(d)+
  facet_wrap(~animal)+
  geom_tile(aes(x=gamma, y=phi, fill=sd_lap))+
  labs(y = expression(paste(phi)),
       x = expression(gamma),
       fill = expression(sigma[r]) ) +
  theme_ipsum(axis_title_size = 18, axis_title_just = "m" ) +
  scale_fill_distiller(palette="RdYlBu")+
  theme(panel.grid.major = element_line(color = "white", size = 0.25),
        panel.grid.minor = element_line(color = "white", size = 0.25),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5))

figS2

ggsave("~/FigS2.png", height = 5, width = 7)


## FIG 3 ---------------------------------------------------------------------
mod_raw <- lm(beta ~ poly(sd_lap, 2, raw = TRUE), data = d)
summary(mod_raw)

dd <- d[sample(1:nrow(d), 0.5*nrow(d), replace = F),]
cols <- brewer.pal(n = 6, name = "Spectral")
cols

fig3 <- ggplot(dd, aes(x = sd_lap, y = beta)) +
  geom_point(aes(color=animal), size=3, pch=21) + 
  facet_wrap(~animal)+
  geom_smooth(method = "lm", formula = y ~ poly(x, 2, raw = T), color = "blue", se = FALSE, linewidth=0.75) +
  labs(y = expression(paste("Scaling exponent ", beta[S])),
       x = expression(paste("Metabolic fluctuation ", sigma[r]))) +
  theme_ipsum(axis_title_size = 18, axis_title_just = "m") +
  scale_color_brewer(palette = "Spectral") +
  theme(panel.grid.major = element_line(color = "white", linewidth = 0.25),
        panel.grid.minor = element_line(color = "white", linewidth = 0.25),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
        legend.position = "none",
        legend.title=element_blank() ) 

fig3

ggsave("~/Dropbox/scaling/figures/fig3.png", height = 5, width = 6)


## FIG S4 ---------------------------------------------------------------------
dat2 <- d %>% group_by(animal) %>% summarise(std = mean(sd_lap),
                                             m = mean(mass),
                                             b = mean(B))
dat2
summary(lm(log10(d$sd_lap)~log10(d$mass)))
summary(lm(log10(dat2$std)~log10(dat2$m)))

figS4a <- ggplot(dat2, aes(x=log10(m), y=log10(std) )) +
  geom_smooth(method="lm", color="black", fill="lightcyan", se=F, size=0.5)+
  geom_point(shape=21, size=2, stroke=0.5)+
  theme_ipsum(axis_title_size = 18, axis_title_just = "m" ) +
  theme(panel.grid.major = element_line(color = "white", linewidth = 0.25),
        panel.grid.minor = element_line(color = "white", linewidth = 0.25),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5))+
  labs(x = "log Mass",
       y = expression(paste("log ", sigma[r]))) +
  geom_text(x = Inf, y = Inf,
            hjust = 1.3, vjust = 3,  # adjust distance from edges
            label = expression(paste( nu, " = -0.305 ", "± 0.038" )),
            color = 'black', size=5, family = "Arial Narrow")

figS4a

summary(lm(log(d$sd_lap)~log(d$B)))
summary(lm(log10(dat2$std)~log10(dat2$b)))

figS4b <- ggplot(dat2, aes(x=log10(b), y=log10(std))) +
  geom_smooth(method="lm", color="black", fill="lightcyan", se=F, size=0.5)+
  geom_point(shape=22, size=2, stroke=0.5)+
  theme_ipsum(axis_title_size = 18, axis_title_just = "m" ) +
  theme(panel.grid.major = element_line(color = "white", linewidth = 0.25),
        panel.grid.minor = element_line(color = "white", linewidth = 0.25),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5))+
  labs(x = "log Metabolic rate",
       y = expression(paste("log ", sigma[r]))) +
  geom_text(x = Inf, y = Inf,
            hjust = 1.3, vjust = 3,  # adjust distance from edges
            label = expression(paste( delta, " = -0.386 ± 0.025")),
            color = 'black', size=5, family = "Arial Narrow")

figS4b

library(cowplot)
setwd("~/Dropbox/scaling/figures/")


grid <- plot_grid( fig6a, fig6b,
                   labels = "auto", 
                   nrow=1,  align="hv", label_size = 20, 
                   label_fontface = "plain")

grid

save_plot("FigS4.png", grid,
          ncol = 2, 
          nrow = 1.1,
          base_aspect_ratio = 1.1 )

### BOOTSTRAP CONF INTERVALS -------------------------------------------
nboot <- 5000  
slopes <- numeric(nboot)
set.seed(123)  

for (b in 1:nboot) {
  boot_sample <- d %>%
    group_by(animal) %>%
    slice_sample(n = 1, replace = TRUE) %>%   
    ungroup()
  fit <- lm(log(B) ~ log(mass), data = boot_sample)
  slopes[b]     <- coef(fit)[2]
}

mean_slope <- mean(slopes)
ci_slope   <- quantile(slopes, c(0.025, 0.975))

list(slope_mean = mean_slope,slope_CI   = ci_slope)

## FIG S3 ---------------------------------------------------------------------
load( "~/sims.RData")

ggplot(df_results)+
  facet_wrap(~animal)+
  geom_point(aes(x=gamma, y=jf, color=best_fit))+
  scale_color_manual(values = c("red", "midnightblue")) +
  labs(y = expression(paste("Baseline fluctuation level ", phi)),
       x = expression(paste("Stocastic inefficiency ", gamma)) ) +
  theme_ipsum(axis_title_size = 16, axis_title_just = "m" ) +
  theme(panel.grid.major = element_line(color = "white", linewidth = 0.25),
        panel.grid.minor = element_line(color = "white", linewidth = 0.25),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
        legend.text=element_text(size=12),
        legend.box.background = element_rect(colour = "black"),
        legend.title = element_blank())

ggsave("~/FigS3.png", height = 5, width = 7.5)







