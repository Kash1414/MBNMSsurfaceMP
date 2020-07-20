# Load pacakges and data
library(tidyverse)
#library(lme4)
library(MASS)
library(MuMIn)


SE = function(x){sd(x, na.rm = TRUE)/sqrt(sum(!is.na(x)))} #function to calculate standard error

col_summ <- function(tbl, f) {
  rbind(
    tbl,
    purrr::map(tbl, ~ if(is.numeric(.x)) f(.x) else NA)
  )
}

LK_data <- read.csv("Sample Particle Counts - Sheet2.csv") %>% #creating new column that removes the underscore and picture number from each ID for grouping
  drop_na(Count)

LK_data 

# create summary tables
LK_data_summ <- LK_data %>% 
  group_by(sampleid) %>% 
  summarise(type = first(Type),
            location = first(Location),
            total_count = sum(Count)*4,
            total_fibers = sum(Fibers)*4,
            total_particles = sum(Particles)*4)
LK_data_summ

groups_sum <- LK_data_summ %>% 
  group_by(type) %>% 
  summarise(mean_count = mean(total_count),
            SE_count = SE(total_count))
groups_sum

#adjust for blanks/water volume
LK_data_summ <- LK_data_summ %>% 
  filter(type %in% c("nearshore", "offshore")) %>% 
  mutate(blank_adj_count = total_count - 22,
         blank_adj_count_zero_corr = ifelse(blank_adj_count < 0, 0, blank_adj_count),
         blank_adj_count_m3 = blank_adj_count/120,
         blank_adj_count_m3_zero_corr = ifelse(blank_adj_count_m3 < 0, 0, blank_adj_count_m3),
         prop_fibers = ifelse(blank_adj_count_m3_zero_corr > 0,total_fibers/total_count, NA),
         prop_frag = ifelse(blank_adj_count_m3_zero_corr > 0, total_particles/total_count, NA))


LK_adj_summ<-LK_data_summ %>% 
  filter(type %in% c("nearshore", "offshore")) %>% 
  #filter(location != "Boardwalk") %>%
  group_by(type) %>% 
  summarise(mean_count = mean(blank_adj_count_m3_zero_corr),
            SE_count = SE(blank_adj_count_m3_zero_corr),
            med_count = median(blank_adj_count_m3_zero_corr))

LK_adj_summ %>% summarise(tot_samp = n_distinct(sampleid))

LK_total_avg<-LK_data_summ %>% 
  summarise(mean_count = mean(blank_adj_count_m3_zero_corr),
            SE_count = SE( blank_adj_count_m3_zero_corr))
LK_total_avg


#Stats here ----

# GLM testing for nearshore v. offshore difference
d_mod <- LK_data_summ %>%  
  filter(sampleid != "B211") %>%  # Removes the one Boardwalk outlier
  filter(type %in% c("nearshore", "offshore")) %>% 
  col_summ(sum)



m_nb_type <- glm.nb(blank_adj_count_m3_zero_corr ~ type, 
                    data = d_mod)

summary(m_nb_type)


## Check for over/underdispersion in the model
E2 <- resid(m_nb_type, type = "pearson")
N  <- nrow(d_mod)
p  <- length(coef(m_nb_type))   
sum(E2^2) / (N - p)


m_p <- glm(blank_adj_count_m3_zero_corr ~ type, 
           family = "poisson",
           data = d_mod)

summary(m_p)

model.sel(m_p, m_nb_location)

# Negative binomial model much better! 



# GLM testing for Sampling Location difference

m_nb_location <- glm.nb(blank_adj_count_m3_zero_corr ~ location -1, 
           data = d_mod)

summary(m_nb_location)


# GLM testing for Sampling Location difference

m_pt_type <- glm(cbind(total_fibers, total_count-total_fibers) ~ type,
                 data = d_mod,
                     family = binomial)
summary(m_pt_type)
r.squaredGLMM(m_pt_type)




#Figure 3 ----
LK_adj_summ_for_plot<-LK_data_summ %>% 
  filter(type %in% c("nearshore", "offshore")) %>% 
  group_by(type,location) %>% 
  summarise(mean_count = mean(blank_adj_count_m3_zero_corr),
            SE_count = SE(blank_adj_count_m3_zero_corr),
            med_count = median(blank_adj_count_m3_zero_corr),
            avg_fiber = mean(prop_fibers, na.rm = TRUE),
            SE_fiber = SE(prop_fibers),
            avg_frag = mean(prop_frag, na.rm = TRUE),
            SE_frag = SE(prop_frag))

pal <- c("Boardwalk" = "olivedrab3","Marina" = "olivedrab4",  "Davidson" = "steelblue4",  "SurRidge" = "steelblue3")


p <- ggplot(LK_adj_summ_for_plot, aes(fill=location, x=type, y=mean_count))
p + geom_bar(position=position_dodge(width=0.9), stat="identity") + 
  geom_errorbar(aes(ymin=mean_count-SE_count, ymax=mean_count+SE_count),
                position=position_dodge(width=0.9), width=0.25) +
  theme_bw() +
  theme(plot.title = element_text(face="bold", vjust=14, size=24)) +
  ylab(bquote('Particles per'~m^3)) +
  scale_fill_manual(values = pal,
                    labels=c("Boardwalk", "Davidson", "Marina", "Sur Ridge")) +
  theme(axis.title.y = element_text(size=16)) +
  theme(axis.text.y  = element_text(vjust=0.5, size=14)) +
  theme(axis.title.x = element_blank()) +
  theme(axis.text.x  = element_text(vjust=0.5, size=12)) +
  theme(legend.text=element_text(size=11)) + 
  theme(legend.title =element_text(size=11)) +
  guides(fill = guide_legend(title = "Sampling\nlocation"))

dev.copy2pdf(file="Particles by site.pdf", width=6, height=4)





# #figure 3: boxplot
# LK_boxplot<-LK_data_summ %>% 
#   filter(type %in% c("nearshore", "offshore")) %>% 
#   ggplot() +
#   geom_boxplot(aes(type, blank_adj_count_m3_zero_corr), outlier.shape = NA) +
#   geom_jitter(aes(type, blank_adj_count_m3_zero_corr, color = location)) +
#   scale_y_log10(labels = scales::comma) +
#   theme_classic() +
#   labs(x="Type", 
#        y=bquote("Number of MP"~(m~3)), 
#        title="Microplastic Particle Distribution")
# 
# LK_boxplot


#Figure 4: particle type by location----
LK_adj_summ_long <- LK_adj_summ_for_plot %>% 
  pivot_longer(cols = c("avg_fiber","avg_frag"), values_to = "avg_part_prop", names_to = "part_type") %>%
  mutate(location = ifelse(as.character(location) == "SurRidge", "Sur Ridge", as.character(location)))

pal2 <- c("avg_fiber" = "darkorchid3", "avg_frag" = "goldenrod3")


p2 <- ggplot(LK_adj_summ_long, 
             aes(fill=part_type, 
                 x=fct_relevel(location, "Boardwalk", "Marina", "Davidson"),
                 y=avg_part_prop))
p2 + geom_bar(position=position_dodge(width=0.9), stat="identity") + 
  geom_errorbar(aes(ymin=avg_part_prop-SE_frag, ymax=avg_part_prop+SE_frag),
                position=position_dodge(width=0.9), width=0.25) +
  facet_grid(.~type, scales = "free") +
  theme_bw() +
  theme(plot.title = element_text(face="bold", vjust=14, size=24)) +
  labs(x = "Sampling location",
       y = "Particle type proportions",
       fill = "Particle\ntype") +
  ylim(0,1) +
  scale_fill_manual(values = pal2,
                    labels=c("fiber", "non-fiber")) +
  #scale_x_discrete(labels=ifelse("SurRidge", "Sur Ridge", location)) +
  theme(axis.title.y = element_text(size=16)) +
  theme(axis.text.y  = element_text(vjust=0.5, size=14)) +
  theme(axis.title.x = element_text(size=16)) +
  theme(axis.text.x  = element_text(vjust=0.5, size=12)) +
  theme(legend.text=element_text(size=11),
        strip.text = element_text(size = 12)) + 
  theme(legend.title =element_text(size=11))


dev.copy2pdf(file="Particle type proportions by site.pdf", width=6, height=4)



# particle_df=mpdata %>%
#   group_by(type) %>%
#   summarise(total=sum(total_count), fiber=sum(total_fibers), particles=sum(total_particles)) %>%
#   mutate(percentfibers=(fiber/total)*100, percentparticles=(particles/total)*100) %>%
#   filter(type %in% c("nearshore","offshore")) %>%
#   pivot_longer(cols = -type, names_to = "plastic", values_to = "stat") %>%
#   filter(plastic %in% c("percentfibers","percentparticles")) %>%
#   mutate(Plastic=fct_recode(plastic, Fibers="percentfibers", Particles="percentparticles")) %>%
#   rename(Percent_Total = "stat") 
# 
# 
# 
# ggplot() +
#   geom_col(aes(type,Percent_Total, fill=Plastic, main="Particle Composition"))
# 
# particle_df
# percent_bar

