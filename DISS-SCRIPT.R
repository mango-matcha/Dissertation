# install packages
install.packages("ggplot2")
install.packages("readr")
install.packages("dplyr")
install.packages("DHARMa")
install.packages("remotes")
install.packages("glmmTMB")
install.packages("vegan")
install.packages("tidyverse")
install.packages("tibble")
install.packages("emmeans")
install.packages("ggeffects")
install.packages("gridExtra")
install.packages("lmerTest")
install.packages("goeveg")
install.packages("devtools")
devtools::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")


#libraries
library(ggplot2) # graphs 
library(readr) #reads in data set
library(DHARMa)    #check assumptions of GLMER
library(vegan) #shannon diveristy
library(tidyverse)#data cleaning, manipuating data
library(emmeans)# negbinomial interpet 
library(glmmTMB) # gamma models
library(pairwiseAdonis) # post hoc for permanova 
library(tibble) #part of NMDS
library(goeveg)#part of NMDS
library(ggeffects) #makes graph from model
library(gridExtra) #arrange graphs next to eachother

#load in data 
final_data <- read_csv("Final-Diss.csv")

#rename data 
final_data<-final_data%>%
  mutate(Species = recode(Species,
                          "CB" = "Mesoporos",
                          "Space" = "Ceratium_other",
                          "spike" = "Ceratium_fusus",
                          "pouch" = "Dinophysis",
                          "star" = "Thalassionema",
                          "2ball" = "Coscinodiscus",
                          "ladder" = "Paralia",
                          "ex" = "Dictyocha",
                          "fluffy" = "Acantharia",
                          "sally" = "copepod_B",
                          "copopod" = "copepod_A",
                          "naupilis"= "copepod_nauplii",
                          "foramifera" = "foraminifera"))%>%
  mutate(Zoo_Phyto= case_when(Species== "Acanthareans"~ "zooplankton",
                              TRUE ~ Zoo_Phyto))


#-----------------------filtering-----------------------------------------------

#scale data to be standardized for sampling methods 
final_data_scaled <- final_data%>%
  mutate(volume = 0.049* depth)%>% #calculates volume for each depth in m^3.Radius of net = 125mm = 0.125m. equation: 0.125^2 *depth(in metres)* pi =0.049*Depth
  mutate(volume = volume * 1000)%>% # converts to litres instead of m^3
  mutate(Abundance_per_ml = if_else(Zoo_Phyto== "phytoplankton",
                                    Abundance/0.03,
                                    Abundance/0.2))%>% #calculates abundance per ml to standardise phyto and zoo
  mutate(Abundance_per_sample= Abundance_per_ml*10*2)%>% #calculates amount per jar bc took 1ml from 20ml sample, 10 to scale to 10 ml of sample *2 bc 20 ml is only half sample half preseration
  mutate(Density = Abundance_per_sample/volume) #calculates total volume(depth towed) (units organisms per litre)



#change coverage, zoovphyto and site to factors to help model
final_data_scaled$coverage <- as.factor(final_data_scaled$coverage)
final_data_scaled$Zoo_Phyto <- as.factor(final_data_scaled$Zoo_Phyto)
final_data_scaled$site <- as.factor(final_data_scaled$site)

#Filter for zoo and phyto separately-> for abudnace model
final_scale <- final_data_scaled %>%
  group_by(site, ID, coverage, Zoo_Phyto,depth) %>%
  summarise(
    Density = sum(Density, na.rm = TRUE),
    .groups = "drop")

#filter for species data set 
species_scale <- final_data_scaled %>%
  group_by(site, ID, coverage, Species) %>%
  summarise(
    Density = sum(Density, na.rm = TRUE),
    .groups = "drop")


#----------------------------Graphing-------------------------------------------
#data set only Zoo -> for 
zoo_scale <- final_scale %>%
  filter(Zoo_Phyto == "zooplankton") %>%
  group_by(site, ID, coverage, depth) %>%
  summarise(
    Density = sum(Density, na.rm = TRUE),
    .groups = "drop")


#data set only phyto 
phyto_scale <- final_scale %>%
  filter(Zoo_Phyto == "phytoplankton") %>%
  group_by(site, ID, coverage, depth) %>%
  summarise(
    Density = sum(Density, na.rm = TRUE),
    .groups = "drop")

# change the order of coverage so goes from low to medium to high
zoo_scale$coverage <- factor(zoo_scale$coverage,
                             levels = c("low", "medium", "high"))

phyto_scale$coverage <- factor(phyto_scale$coverage,
                               levels = c("low", "medium", "high"))

#Boxplot of Zooplankton vs kelp coverage
(zooplankton <- ggplot(zoo_scale, aes(x = coverage, y = Density, fill = coverage)) +
    geom_boxplot(fill = "#b31529") +
    geom_point(fill = "#b31529") +
    labs(y = "Zooplankton Density per Litre (Organisms Per Litre)",
         x = "Kelp Coverage") +
    ylim(0, 1336.883) +
    scale_y_continuous(labels = scales::comma, limits = c(0, 1336.883)) +
    theme_classic())


#Boxplot of Phytoplankton vs Kelp Coverage
(phytoplankton <- ggplot(phyto_scale, aes(x = coverage, y = Density, fill = coverage)) +
    geom_boxplot(fill = "#8ec3de") +
    geom_point(fill = "#8ec3de") +
    labs(y = "Phytoplankton Density (Organisms Per Litre) ",
         x = "Kelp Coverage") +
    ylim(0, 1336.883) + # uses a mutual range to put both plots on same scale
    scale_y_continuous(labels = scales::comma, limits = c(0, 1336.883)) + # this changes the y axis to Real numbers 
    theme_classic())

#arrange the plots together
raw_plot<- grid.arrange(zooplankton, phytoplankton, ncol=2)

ggsave("raw_plot.png", plot = raw_plot, width = 12, height = 8, dpi = 600) #saves modelled plot as png 

# check the full range across both datasets to find mutual range to plot on same scale
range(c(zoo_scale$Density, phyto_scale$Density))

#getting descriptive statistics for phytoplankton
mean_phyto<- phyto_scale %>%
  group_by(coverage)%>%
  summarise(
    mean = mean(Density),
    median = median(Density),
    sd = sd(Density),
    min = min(Density),
    max = max(Density),
    Q1 = quantile(Density, 0.25),
    Q3 = quantile(Density, 0.75),
    n = n()
  )
#getting descriptive statistics for zooplankton
mean_zoo<- zoo_scale %>%
  group_by(coverage)%>%
  summarise(
    mean = mean(Density),
    median = median(Density),
    sd = sd(Density),
    min = min(Density),
    max = max(Density),
    Q1 = quantile(Density, 0.25),
    Q3 = quantile(Density, 0.75),
    n = n()
  )


#species density 
species_density_each <- final_data_scaled %>%
  group_by(Species) %>%
  summarise(
    total_density = sum(Density, na.rm = TRUE),
    mean_density = mean(Density, na.rm = TRUE),
    n = n()) 

print(species_density_each)

#---------------------Modeling Abundance----------------------------------------

#tweedie link: becuase has zeros, non integers, continous (due to scaling), not normal, data overdispersed
#AIC better then doing a gamma log link and changing zeros* Tweedie is the best option*

#combination of both zoo and phyto
zoovphyto <- glmmTMB(Density ~ coverage * Zoo_Phyto + (1|site),
                     family = tweedie(link = "log"), data=final_scale)

#graph residuals 
simulateResiduals(fittedModel = zoovphyto, plot = T) 

#overall statistics
summary(zoovphyto) 

#confidence intervals
confint(zoovphyto)

#gives values already exp
emmeans(zoovphyto, ~ coverage * Zoo_Phyto, type = "response")

#emmeans tests
#overall comparison
joint_tests(zoovphyto)

#each comparison indivudally
em <- emmeans(zoovphyto, ~ coverage * Zoo_Phyto)
pairs(em, adjust = "tukey", type = "response")

#each comparison with confidence intervals
pairs(em, adjust = "tukey") %>% 
  confint(type = "response")


#graph modeled results
final_scale$coverage <- factor(final_scale$coverage,
                               levels = c("low", "medium", "high"))

predicted <- ggpredict(zoovphyto, terms = c("coverage", "Zoo_Phyto"))

density<- (ggplot(predicted, aes(x = x, y = predicted))+
             geom_jitter(data = final_scale, #plots raw data 
                         aes(x = coverage,
                             y = Density,
                             color = Zoo_Phyto),  
                         width = 0.2,
                         alpha = 0.3) +
             
             geom_point(data = predicted, # plots predicted data
                        aes(x = x,
                            y = predicted,
                            color = group),       
                        size = 4,
                        position = position_dodge(width = 0.3)) +
             
             geom_errorbar(data = predicted, # predicted error bars
                           aes(x = x,
                               ymin = conf.low,
                               ymax = conf.high,
                               color = group), 
                           width = 0.2,
                           position = position_dodge(width = 0.3)) + #this seperates the two lines so zoo and phyto r diff position
             scale_color_manual(values = c("phytoplankton" = "#8ec3de",
                                           "zooplankton" = "#b31529")) +
             
             labs(x = "Kelp Coverage", # names axis 
                  y = "Predicted Plankton Density per Litre",
                  color = "Plankton Type") +
             theme_classic())

ggsave("density_plot.png", plot = density, width = 12, height = 8, dpi = 600) #saves modelled plot as png 

#
gamma_scale <- final_data%>%
  mutate(volume = 0.049* depth)%>% #calculates volume for each depth in m^3.Radius of net = 125mm = 0.125m. equation: 0.125^2 *depth(in metres)* pi =0.049*Depth
  mutate(volume = volume * 1000)%>% # converts to litres instead of m^3
  mutate(Abundance_per_ml = if_else(Zoo_Phyto== "phytoplankton",
                                    Abundance/0.03,
                                    Abundance/0.2))%>% #calculates abundance per ml to standardise phyto and zoo
  mutate(Abundance_per_sample= Abundance_per_ml*10*2)%>% #calculates amount per jar bc took 1ml from 20ml sample, 10 to scale to 10 ml of sample *2 bc 20 ml is only half sample half preseration
  mutate(Density = Abundance_per_sample/volume)  #calculates total volume(depth towed) (units organisms per litre)


#--------------------GAMMA Density Model -----------------------------------------------------------
#MODEL NOT CHOSEN
#final_gam <- final_data_scaled %>%
#  group_by(site, ID, coverage, Zoo_Phyto,depth) %>%
#  summarise(
#    Density = sum(Density, na.rm = TRUE),
#    .groups = "drop")%>%
#  mutate(Density = if_else(Density == 0, 1e-6, Density)) #adjust zeros 

#modela <- glmmTMB(Density ~ coverage * Zoo_Phyto + (1|site),
#                  family = Gamma(link = "log"),
#                  data = final_gam)
#graph residuals 
#simulateResiduals(fittedModel = modela, plot = T) 

#overall statistics
#summary(modela) 

#confidence intervals
#confint(modela)


#--------------------overall -----------------------------------------------------------
#Filter all together
final_scale_sum <- final_data_scaled %>%
  group_by(site, ID, coverage, depth) %>%
  summarise(
    Density = sum(Density, na.rm = TRUE),
    .groups = "drop")

final_scale_sum$coverage <- factor(final_scale_sum$coverage,
                                   levels = c("low", "medium", "high"))


#Model of overall abundance
overall <- glmmTMB(Density ~ coverage + (1|site),
                   family = Gamma(link = "log"), data =final_scale_sum)

#plot residuals
simulateResiduals(fittedModel = overall, plot = T) 

#overall summary of statisics 
summary(overall) 

#emmeans overall comparison 
emmeans(overall, ~ coverage, type = "response")

#emmeans for each indivudal coverage comparison
pairs(emmeans(overall, ~ coverage, type = "response"))


#graph results
predicted_2 <- ggpredict(overall, terms = "coverage")

overall_plot<- (ggplot(predicted_2, aes(x = x, y = predicted))+
                  geom_jitter(data = final_scale_sum, #plots raw data 
                              aes(x = coverage,
                                  y = Density),  
                              width = 0.2,
                              color ="#7CC4FF",
                              alpha = 0.3) +
                  
                  geom_point(data = predicted_2, # plots predicted data
                             aes(x = x,
                                 y = predicted),       
                             size = 4,
                             color  ="#0B3C5D") +
                  
                  geom_errorbar(data = predicted_2, # predicted error bars
                                aes(x = x,
                                    ymin = conf.low,
                                    ymax = conf.high), 
                                width = 0.2,
                                color  ="#0B3C5D",
                                position = position_dodge(width = 0.3)) + #this seperates the two lines so zoo and phyto r diff position
                  
                  
                  labs(x = "Kelp Coverage", # names axis 
                       y = "Predicted Plankton Density (Organsims Per Litre) ") +
                  theme_classic())

overall_plot


ggsave("overall_plot.png", plot = overall_plot, width = 12, height = 8, dpi = 600) #saves modelled plot as png 

#------------------diversity model---------------------------------------------- 
#filter to make diversity set

diversity_data <- final_data_scaled %>%
  dplyr::select(ID, site, coverage, Species, Density)%>%
  group_by(site, coverage)%>%
  summarise(Diversity = vegan::diversity(Density , index = "shannon"),
            .groups = "drop")


#plot diversity
diversity_data$coverage <- factor(diversity_data$coverage,
                                  levels = c("low", "medium", "high"))

(div <- ggplot(diversity_data, aes(x = coverage, y = Diversity)) +
    geom_boxplot(fill = "#8ec3de") +
    labs(y = "Plankton Diversity (Shannon Diveristy Indices)",
         x = "Kelp Coverage") +
    theme_classic())
ggsave("div_plot.png", plot = div, width = 12, height = 8, dpi = 600) #saves modelled plot as png 



mean_div<- diversity_data %>%
  group_by(coverage)%>%
  summarise(
    mean = mean(Diversity),
    median = median(Diversity),
    sd = sd(Diversity),
    min = min(Diversity),
    max = max(Diversity),
    Q1 = quantile(Diversity, 0.25),
    Q3 = quantile(Diversity, 0.75),
    n = n()
  )


#model of diversity to test signficance
model_diversity <- glm(Diversity ~ coverage, family=Gamma(link = "log"), data = diversity_data2)

#plot residuals 
simulateResiduals(fittedModel = model_diversity, plot = T) 

#overall summary of stats
summary(model_diversity)

#emmeans
emmeans(model_diversity, ~ coverage, type = "response")

#-------------------NMDS--------------------------------------------------------
#simper analysis 

com_matrix <- final_data_scaled %>% 
  dplyr::select(ID, Species, Density)%>%
  # Turn our Site, Species, and Forest columns into factors
  group_by(ID, Species)%>%
  mutate_at(vars(ID), as.factor) %>%  
  pivot_wider(names_from = Species, # pivot to wide format
              values_from = Density) %>% 
  column_to_rownames(var = "ID") # change our column "site" to our rownames


nmds <-
  metaMDS(com_matrix,
          distance = "bray",
          k = 2)

print(nmds)

stressplot(nmds)

dimcheck_out <- 
  dimcheckMDS(com_matrix,
              distance = "bray",
              k = 6)

print(dimcheck_out)
plot(nmds)

coverage_type <- final_data %>% 
  distinct(ID, coverage)


# Extract NMDS scores for sites 
nmds_SiteScores <-
  as.data.frame(scores(nmds, display ="sites")) %>%  # get nmds scores 
  rownames_to_column("ID")%>% # change rownames (site) to a column 
  left_join(coverage_type, by="ID")   # join our habitat type (grouping variable) to each site 


# Extract NMDS scores for species 
nmds_SpeciesScores <-  
  as.data.frame(scores(nmds, "species"))

# create a column of species, from the rownames of species.scores
nmds_SpeciesScores$species <- rownames(nmds_SpeciesScores) 

#plot NMDS
ggplot() + 
  geom_point(data = nmds_SiteScores, 
             aes(x=NMDS1, y=NMDS2, colour = coverage), 
             size = 3) + 
  
  geom_text(data = nmds_SpeciesScores, 
            aes(x=NMDS1, y=NMDS2, label = species),
            size = 3) +
  theme_classic()


coverage_Centroid <- 
  nmds_SiteScores %>% 
  group_by(coverage) %>% 
  summarise(axis1 = mean(NMDS1),
            axis2 = mean(NMDS2)) %>% 
  ungroup()

# extract convex hull
coverage.hull <- 
  nmds_SiteScores %>% 
  group_by(coverage) %>%
  slice(chull(NMDS1, NMDS2))

nmds_stress <- nmds$stress


#final plot 
# use ggplot to plot 
nmds_plot <-(ggplot() + 
               
               # add site scores
               geom_point(data = nmds_SiteScores, 
                          aes(x=NMDS1, y=NMDS2, colour = coverage), size = 2) + 
               
               # add species scores 
               geom_text(data = nmds_SpeciesScores, 
                         aes(x=NMDS1, y=NMDS2, label = species)) +
               
               # add centroid 
               geom_point(data = coverage_Centroid, 
                          aes(x = axis1, y = axis2, color = coverage), 
                          size = 5, shape = 17) +
               
               # add convex hull
               geom_polygon(data = coverage.hull, 
                            aes(x = NMDS1, y = NMDS2, fill = coverage, group = coverage), 
                            alpha = 0.30) +
               
               # add stress value
               annotate("text", x = 0.75, y = 0.65, 
                        label = paste("2d stress =", round(nmds_stress, 3))) +
               
               scale_color_manual(values = c("low" = "#EE7733",
                                             "high" = "#66CCEE",
                                             "medium" = "#EE3377")) +
               scale_fill_manual(values = c("low" = "#EE7733",
                                            "high" = "#66CCEE",
                                            "medium" = "#EE3377"))+
               
               # edit theme
               labs(x = "NMDS1", y = "NMDS2") + 
               theme(axis.text = element_blank(), axis.ticks = element_blank(),
                     panel.background = element_rect(fill = "white"),
                     panel.border = element_rect(color = "black", 
                                                 fill = NA, linewidth = .5),
                     axis.line = element_line(color = "black"),
                     plot.title = element_text(hjust = 0.5),
                     legend.key.size = unit(.25, "cm")))

ggsave("nmds.png", plot = nmds_plot, width = 12, height = 8, dpi = 600) #saves modelled plot as png 

#permanova
adonis2(com_matrix ~ coverage_type$coverage, data = coverage_type)

#more indepth results 
pairwise.adonis2(com_matrix ~ coverage, data = coverage_type)



#-------------- Which species most abudnant in coverage---------------
#Filter for zoo and phyto separately-> for comparing species
final_scale_2 <- final_data_scaled %>%
  group_by(site, ID, coverage, Zoo_Phyto, Species) %>%
  summarise(
    Density = sum(Density, na.rm = TRUE),
    .groups = "drop")

final_scale_2$coverage <- factor(final_scale_2$coverage,
                             levels = c("low", "medium", "high"))

#mean of all species 
species_summary <- final_scale_2 %>%
  group_by(coverage, Species) %>%
  summarise(mean_density = mean(Density),
            total_density = sum(Density)) %>%
  arrange(coverage, desc(mean_density))

#get top 5 species
top_species <- species_summary %>%
  group_by(coverage) %>%
  slice_max(mean_density, n = 5)

#plot top 5 species
top<-(ggplot(top_species, aes(x = reorder(Species, mean_density), 
                              y = mean_density, 
                              fill = Species)) +
        geom_bar(stat = "identity") +
        scale_fill_manual(values = c("Mesoporos" = "#8ec3de",
                                     "Ceratium_other"="#8ec3de",
                                     "Ceratium_fusus"= "#8ec3de",
                                     "Dinophysis"= "#8ec3de",
                                     "copepod_A" = "#b31529")) +
        facet_wrap(~ coverage, scales = "free_y") +
        coord_flip() +
        labs(x = "Species", 
             y = "Mean Density",
             title = "Top 5 Most Abundant Species per Kelp Coverage") +
        theme_classic() +
        theme(legend.position = "none"))


ggsave("top.png", plot = top, width = 12, height = 8, dpi = 600) #saves modelled plot as png 


#simper analysis 
simper_coverage <- simper(com_matrix, coverage_type$coverage)
summary(simper_coverage)

