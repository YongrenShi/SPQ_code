## Author: Yongren Shi
## Paper: Lynn, Freda B., Yongren Shi, and Kevin Kiley. "Intersectional Group Agreement on the Occupational Order.” Social Psychology Quarterly


library(haven)
library(ggplot2)
library(tidyverse)
library(labelled)
library(cluster)
library(hrbrthemes)
library(parallel)
library(ggpubr)
library(igraph)
library(RColorBrewer)
library(ForceAtlas2) #devtools::install_github("analyxcompany/ForceAtlas2")

setwd("set the directory to folder")

split_tibble <- function(tibble, col = 'col') tibble %>% split(., .[, col])
colMax <- function(data) sapply(data, max, na.rm = TRUE)
colMin <- function(data) sapply(data, min, na.rm = TRUE)
color_gradient <- colorRampPalette(brewer.pal(9, "Blues"))


resampling <- function(gss_demog){ # bootstrap resampling
  resampled_gss_demog <- gss_demog %>% mutate(across(Accountant:Welder, ~ sample(.x, replace=TRUE)))
  return(resampled_gss_demog)
}

pairwise_manhattan_distance <- function(df_pair){ # across all occupations
  i = df_pair$Var1; j = df_pair$Var2
  dis_vec = abs(g_mat[i,] - g_mat[j,])
  n_occ = sum(!is.na(dis_vec))
  dis = ( sum(dis_vec, na.rm=TRUE) ) / (n_occ)
  return(data.frame(i, j, dis))
}

pairwise_spearman_similarity <- function(df_pair){ # across all occupations
  i = df_pair$Var1; j = df_pair$Var2
  df_ = data.frame("i"=g_mat[i,], "j"=g_mat[j,])
  df_ = df_ %>% filter(!is.na(i) & !is.na(j))
  sim = NA
  if (nrow(df_)>0) sim = cor(g_mat[i,], g_mat[j,], method = "spearman", use="pairwise.complete.obs")
  return(data.frame(i, j, sim))
}

Consensus_point_estimation <- function(g_df, pairwise_method){
  #g_df = gss_demog %>% filter(sex=="Male" & race=="White") #(educ == "Graduate" & race=="White")
  df_pairs <- expand.grid(unique(g_df$new_id), unique(g_df$new_id)) %>% filter(Var1!=Var2) %>% rename("i"="Var1", "j"="Var2")
  df_pairs = pairwise_result_df %>% inner_join(df_pairs, by=c("i", "j"))
  if (pairwise_method == "spearman") {consensus =  mean(df_pairs$sim, na.rm=TRUE)}
  if (pairwise_method == "manhattan") {
    E = 2.962
    consensus = (E - mean(df_pairs$dis, na.rm=TRUE)) / E
  }
  return(consensus)
}

# finding stretch of a hierarchy
group_Hierarchy_estimation <- function(g_df){
  return(g_df %>% select(Accountant:Welder) %>% colMeans(na.rm=TRUE) %>% sd())
}  

# intersections homogeneity
calc_intersect_H <- function(gss_demog){
  g_homogeneity = data.frame()
  for (degree_ in c("No Diploma" , "High School", "Junior College", "Bachelor's", "Graduate"))
    for (sex_ in c("Male", "Female"))
      for (race_ in c("White", "Black", "Other")) {# for (class_ in c(1,2,3,4))
        g_df = gss_demog %>% filter(degree==degree_ & sex==sex_ & race==race_)
        if (nrow(g_df)>1){
          homogeneity = Consensus_point_estimation(g_df, pairwise_method)
          hierarchy <- group_Hierarchy_estimation(g_df)
          g_homogeneity <- rbind(g_homogeneity, c(degree_, sex_, race_, homogeneity, hierarchy, nrow(g_df)))
        }
      }
  colnames(g_homogeneity) = c("Education", "Sex","Race",  "homogeneity", "hierarchy" ,"size")
  g_homogeneity <- g_homogeneity %>% mutate("intersects" = rownames(g_homogeneity)) 
  g_homogeneity <- g_homogeneity %>% mutate(grouping = paste( Race, Sex, Education, sep=", "))
  return(g_homogeneity)
}

# base groups homogeneity
calc_base_H <- function(gss_demog){
  g_base_homogeneity = data.frame()
  for (base in c("age", "degree", "sex", "race", "religion", "polparty", "region", "socialclass")){
    if (base=="age") for(age_ in c("18-35","36-65", "65+")){
      g_df = gss_demog %>% filter(age == age_)
      homogeneity = Consensus_point_estimation(g_df, pairwise_method)
      hierarchy <- group_Hierarchy_estimation(g_df)
      g_base_homogeneity <- rbind(g_base_homogeneity, c(age_, homogeneity, hierarchy, nrow(g_df)))
    }
    if (base=="degree") for(degree_ in c("No Diploma" , "High School", "Junior College", "Bachelor's", "Graduate")){
      g_df = gss_demog %>% filter(degree == degree_)
      homogeneity = Consensus_point_estimation(g_df, pairwise_method)
      hierarchy <- group_Hierarchy_estimation(g_df)
      g_base_homogeneity <- rbind(g_base_homogeneity, c(degree_, homogeneity, hierarchy, nrow(g_df)))
    }
    if (base=="sex") for(sex_ in c("Male", "Female")){
      g_df = gss_demog %>% filter(sex == sex_)
      homogeneity = Consensus_point_estimation(g_df, pairwise_method)
      hierarchy <- group_Hierarchy_estimation(g_df)
      g_base_homogeneity <- rbind(g_base_homogeneity, c(sex_, homogeneity, hierarchy,  nrow(g_df)))
    }
    if (base=="race") for(race_ in c("White", "Black", "Other")){
      g_df = gss_demog %>% filter(race == race_)
      homogeneity = Consensus_point_estimation(g_df, pairwise_method)
      hierarchy <- group_Hierarchy_estimation(g_df)
      g_base_homogeneity <- rbind(g_base_homogeneity, c(race_, homogeneity, hierarchy, nrow(g_df)))
    }
    if (base=="religion") for(relig_ in c("Protestant", "Catholic", "Atheist", "Other Religion", "Jewish")){
      g_df = gss_demog %>% filter(relig == relig_)
      homogeneity = Consensus_point_estimation(g_df, pairwise_method)
      hierarchy <- group_Hierarchy_estimation(g_df)
      g_base_homogeneity <- rbind(g_base_homogeneity, c(relig_, homogeneity, hierarchy, nrow(g_df)))
    }
    if (base=="polparty") for(politics_ in c("Democrat", "Republican", "Independent")){
      g_df = gss_demog %>% filter(polparty == politics_)
      homogeneity = Consensus_point_estimation(g_df, pairwise_method)
      hierarchy <- group_Hierarchy_estimation(g_df)
      g_base_homogeneity <- rbind(g_base_homogeneity, c(politics_, homogeneity, hierarchy, nrow(g_df)))
    }
    if (base=="region") for(region_ in c("South", "Not South")){
      g_df = gss_demog %>% filter(region == region_)
      homogeneity = Consensus_point_estimation(g_df, pairwise_method)
      hierarchy <- group_Hierarchy_estimation(g_df)
      g_base_homogeneity <- rbind(g_base_homogeneity, c(region_, homogeneity, hierarchy, nrow(g_df)))
    }
    if (base=="socialclass") for(socialclass_ in c("Lower Class", "Working Class", "Middle Class", "Upper Class")){
      g_df = gss_demog %>% filter(socialclass == socialclass_)
      homogeneity = Consensus_point_estimation(g_df, pairwise_method)
      hierarchy <- group_Hierarchy_estimation(g_df)
      g_base_homogeneity <- rbind(g_base_homogeneity, c(socialclass_, homogeneity, hierarchy, nrow(g_df)))
    }
  }
  colnames(g_base_homogeneity) = c("grouping", "homogeneity", "hierarchy", "size")
  return(g_base_homogeneity)
}

##############################################################################
## Start of the main analysis
pairwise_method = "manhattan" # manhattan spearman
year = 1989 #2012 

if (year==1989){
  occ_labels <- read.csv("orgprestg_data/occupation_labels.csv")
  occupres_df  <- read_dta("orgprestg_data/ICPSR_09593.dta")
  occupres <- occupres_df %>% select(akk:wel)
}
if (year==2012){
  occ_labels <- read.csv("orgprestg_data/prestg2012_crosswalk.csv") %>% filter(Core==1) %>% rename(label=GSS.occupation.titles, abbre=Code)%>%
    mutate(abbre = tolower(abbre), label = gsub("\\s+$", "", label))
  occupres_df  <- read.csv("orgprestg_data/GSS_prestg_2012.csv") %>% rename(gssid = id) 
  colnames(occupres_df) = tolower(colnames(occupres_df))
  occupres = occupres_df %>% select(occ_labels$abbre)
}
occ_abrv = occ_labels$abbre
colnames(occupres) <- occ_labels$label
occupres[occupres>10] = NA
occupres["n_ratings"] <- (!is.na(occupres)) %>% rowSums()
occupres_df <- cbind(occupres, id=occupres_df$gssid)
occupres_df <- occupres_df %>% filter(n_ratings>0)
# summarizing the average occpuations that the respondents rated.
summary(occupres_df$n_ratings)
sd(occupres_df$n_ratings)
mean(occupres_df$n_ratings)

# finding the mean rating for the core occupations
row_means <- occupres_df %>% select(Accountant:Welder) %>% rowwise() %>%
  mutate(mean_of_row = mean(c_across(everything()), na.rm = TRUE)) %>%
  pull(mean_of_row)
mean(row_means)

occ_character_df  <- read_dta("orgprestg_data/occ_characteristics.dta") %>% inner_join(data.frame(threelettercode=occ_abrv), by="threelettercode") %>% 
  arrange(threelettercode) %>% select(threelettercode, TRAIN) %>% mutate(occupation=occ_labels$label)

if (year==1989) {
  gss_data = read_dta("gssdata/GSS1989.dta")
  gss_data = gss_data %>% rename(income_used = income86) %>% mutate(income_used = ifelse(income_used > 18, 1, 0)) # >$60000 (11.9)
}
if (year==2012) { 
  gss_data = read_dta("orgprestg_data/GSS2012merged_r10 - stata.dta")
  gss_data = gss_data %>% rename(income_used = income06) %>% mutate(income_used = ifelse(income_used > 22, 1, 0))#$ >$90000 ~ 15.6%
  
}
# recode classification systems
gss_demog <- gss_data %>% select(id, sex, race, educ, attend, region, partyid, age, class, prestg10, relig, degree, income_used) %>%
  mutate(sex=as.factor(sex), 
         race=as.factor(race), 
         degree = ifelse(degree==0, "No Diploma", 
                         ifelse(degree==1, "High School",
                                ifelse(degree==2, "Junior College", 
                                       ifelse(degree==3, "Bachelor's",
                                              ifelse(degree==4, "Graduate", NA))))), 
         educ = ifelse(educ<12, "No Diploma", 
                       ifelse(educ==12, "High School", 
                              ifelse(educ>12 & educ<16, "Associate Degree", 
                                     ifelse(educ==16, "Bachelor's", "Graduate")))),
         region_original = region, 
         region = ifelse(region %in% c(5, 6, 7), "South", "Not South"), 
         attend = ifelse(attend %in% c(0), "Never Attend", 
                         ifelse(attend %in% c(1,2,3,4,5), "Occasionally Attend", 
                                ifelse(attend %in% c(6,7,8), "Frequently Attend", NA))),
         relig = ifelse(relig==1, "Protestant", 
                        ifelse(relig==2, "Catholic",
                               ifelse(relig==3, "Jewish", 
                                      ifelse(relig==4, "Atheist",
                                             ifelse(relig==5, "Other Religion", NA))))), 
         age = ifelse(age <= 35, "18-35", 
                      ifelse(age>35 & age<=65, "36-65", 
                             ifelse(age >65, "65+", NA))),
         prestg10_num = prestg10,
         
         prestg10 = ifelse(prestg10<40, "Occ_0_39", 
                           ifelse(prestg10<60, "Occ_40_59", "Occ_60_100")),
         polparty = ifelse(partyid %in% c(0, 1), "Democrat", 
                           ifelse(partyid %in% c(5, 6), "Republican", 
                                  ifelse(partyid %in% c(2, 3, 4, 7), "Independent", NA))), 
         socialclass = ifelse(class == 1, "Lower Class", 
                              ifelse(class == 2, "Working Class", 
                                     ifelse(class == 3, "Middle Class", 
                                            ifelse(class==4, "Upper Class", NA)))),
  )
gss_demog$sex <-  recode(gss_demog$sex, "1" = "Male", "2"="Female", .default = NA_character_)
gss_demog$race <-  recode(gss_demog$race, "1" = "White", "2"="Black", "3"="Other", .default = NA_character_)

grp_size = gss_demog %>% group_by(race, sex, degree) %>% summarise(gss_size = n()) %>% ungroup()
grp_size = grp_size %>% mutate(g1 = paste(race, sex, degree, sep = ", ")) %>% select(g1, gss_size)

gss_demog = occupres_df %>% inner_join(gss_demog, by="id")
gss_demog <- gss_demog %>% mutate(new_id = 1:nrow(gss_demog))



################ calculate within-group consensus in their ratings across all occupations ##############
g_mat = gss_demog %>% select(Accountant:Welder) %>% as.matrix()
df_pairs <- expand.grid(1:nrow(g_mat), 1:nrow(g_mat)) %>% filter(Var1!=Var2) %>% mutate(idx = 1:n())
pair_list <- df_pairs %>% split_tibble(c("idx"))
if (pairwise_method == "spearman"){pairwise_result_list <- mclapply(pair_list, pairwise_spearman_similarity, mc.cores = 32)}
if (pairwise_method == "manhattan"){pairwise_result_list <- mclapply(pair_list, pairwise_manhattan_distance, mc.cores = 32)}
pairwise_result_df = as.data.frame(data.table::rbindlist(pairwise_result_list))

g_homogeneity <- calc_intersect_H(gss_demog)
g_base_homogeneity <- calc_base_H(gss_demog)
g_homogeneity <- rbind(g_homogeneity %>% select("grouping", "homogeneity", "hierarchy", "size"), g_base_homogeneity)
g_homogeneity <- g_homogeneity %>% mutate(size = as.numeric(size), homogeneity = as.numeric(homogeneity),  hierarchy = as.numeric(hierarchy))


#### bootstrapping sampling distributions
n_bs = 1000
df_pairs <- expand.grid(1:nrow(g_mat), 1:nrow(g_mat)) %>% filter(Var1!=Var2) %>% mutate(idx = 1:n())
pair_list <- df_pairs %>% split_tibble(c("idx"))
for (counter in 1:n_bs){
  print(counter)
  gss_demog_resampled = resampling(gss_demog)
  g_mat = gss_demog_resampled %>% select(Accountant:Welder) %>% as.matrix()
  if (pairwise_method == "spearman"){pairwise_result_list <- mclapply(pair_list, pairwise_spearman_similarity, mc.cores = 32)}
  if (pairwise_method == "manhattan"){pairwise_result_list <- mclapply(pair_list, pairwise_manhattan_distance, mc.cores = 32)}
  pairwise_result_df = as.data.frame(data.table::rbindlist(pairwise_result_list))
  
  g_bs_homogeneity <- calc_intersect_H(gss_demog_resampled)
  g_bs_base_homogeneity <- calc_base_H(gss_demog_resampled)
  g_bs_homogeneity <- rbind(g_bs_homogeneity %>% select("grouping", "homogeneity", "hierarchy","size"), g_bs_base_homogeneity)
  g_homogeneity[paste0("hm_var", counter)] = g_bs_homogeneity$homogeneity
  g_homogeneity[paste0("hr_var", counter)] = g_bs_homogeneity$hierarchy
}
write.csv(g_homogeneity, paste0("bs_g_consensus_MS", year, "_", pairwise_method,".csv"))




g_homogeneity = read.csv(paste0("bs_g_consensus_MS", year, "_", pairwise_method,".csv"))
g_homogeneity_outcome = g_homogeneity %>% select(grouping, homogeneity, size, starts_with("hm_")) 

g_homogeneity_outcome = g_homogeneity_outcome %>% mutate(across(homogeneity:hm_var100, ~ as.numeric(.x))) %>%
  mutate(se = apply(g_homogeneity_outcome %>% select(hm_var1:hm_var100), 1, sd)) %>% 
  select(grouping, homogeneity, size, se) 


g_homogeneity_outcome <- g_homogeneity_outcome %>% mutate(upper = homogeneity + 1.96*se, lower = homogeneity - 1.96*se)
max_ = max(g_homogeneity_outcome[g_homogeneity_outcome$size>=5,]$upper, na.rm=TRUE)
min_ = min(g_homogeneity_outcome[g_homogeneity_outcome$size>=5,]$lower, na.rm=TRUE)

g_homogeneity_outcome <- g_homogeneity_outcome %>% arrange(homogeneity)
g_homogeneity_outcome <- g_homogeneity_outcome %>% mutate(homogeneity = as.numeric(homogeneity), 
                                                          size = as.numeric(size), 
                                                          grouping = factor(grouping, levels=g_homogeneity_outcome$grouping))
g_homogeneity_outcome %>% ggplot(aes(x=grouping, y=homogeneity)) + geom_point() +  ggpubr::theme_pubclean() +
  geom_errorbar( aes(x=grouping, ymin=lower, ymax=upper))+coord_flip()

# Figure 2 Within-Group Consensus. 
g_educ_intersect = data.frame()
for (education_ in c("No Diploma" , "High School", "Junior College", "Bachelor's", "Graduate"))
  for (sex_ in c("Male", "Female"))
    for (race_ in c("White", "Black")) {# for (class_ in c(1,2,3,4))
      g_educ_intersect <- rbind(g_educ_intersect, c(education_, sex_, race_))
    }

colnames(g_educ_intersect) <- c("Education", "Sex","Race")
g_educ_intersect <- g_educ_intersect %>% mutate(old_grouping = paste(Race, Sex, Education, sep=", "),
                                                print_grouping = paste(Race, Sex, sep=", ")) 

g_educ_intersect <- g_educ_intersect %>% left_join(g_homogeneity_outcome, by=c("old_grouping"="grouping"))

g_educ_intersect <- g_educ_intersect %>% mutate(Education = factor(Education, levels=c("Graduate", "Bachelor's", "Junior College", "High School", "No Diploma")))
g_educ_intersect <- g_educ_intersect %>%  arrange(Education, (homogeneity))

g_educ_intersect <- g_educ_intersect %>% mutate(old_grouping = factor(old_grouping, levels=g_educ_intersect$old_grouping), 
                                                print_grouping = paste0(print_grouping, " (n=", size, ")"))

g_vertical = g_homogeneity_outcome %>% filter(grouping %in% c("No Diploma" , "High School", "Junior College", "Bachelor's", "Graduate")) %>%
  rename("Education" = "grouping")

g_educ_intersect = g_educ_intersect %>% filter(size>=5)
min_c = min(g_educ_intersect$homogeneity); max_c = max(g_educ_intersect$homogeneity)
g_educ_intersect = g_educ_intersect %>% mutate(normalized_values = (g_educ_intersect$homogeneity - min_c) / (max_c - min_c))
index_ = findInterval(g_educ_intersect$normalized_values, seq(0, 1, length.out = length(g_educ_intersect$normalized_values) ))
index_[which.min(index_)] = index_[which.min(index_)] 
index_[which.max(index_)] = index_[which.max(index_)] 
attribute_colors <- color_gradient(nrow(g_educ_intersect)+1)[index_]


g_educ_intersect %>%
  ggplot(aes(x=old_grouping, y=homogeneity)) +
  geom_point(aes(fill=factor(old_grouping)), shape=24, stroke=0.2, size=4, color="black") +  # Set stroke color
  scale_fill_manual(values = attribute_colors) +  # Use scale_fill_manual for fill colors
  ggpubr::theme_pubclean() +
  geom_errorbar(aes(x=old_grouping, ymin=lower, ymax=upper), width = 0.25) + coord_flip() +
  ggforce::facet_col(vars(Education), scales = "free_y", space = "free") +
  geom_hline(data = g_vertical, mapping = aes(yintercept = homogeneity), linetype = "dashed") +
  theme(axis.title.y=element_blank(), axis.text.y = element_text(size = 10, color = "black"), legend.position = "none") +
  ylab("Within-Group Consensus") + 
  scale_x_discrete(breaks=g_educ_intersect$old_grouping, labels = g_educ_intersect$print_grouping) +
  ylim(min_, 0.6) # switch max_ to 0.6 only for manhattan
ggsave(paste0("Panel2_MS_NoOthers_", year, "_", pairwise_method,".pdf"),   width = 5, height=6, units="in", dpi= 300)

################################################################################################
## S3 in supplement Within-Group Consensus for Single-Characteristic Groups
g_homogeneity = read.csv(paste0("bs_g_consensus_MS", year, "_", pairwise_method,".csv"))

g_homogeneity_outcome = g_homogeneity %>% select(grouping, homogeneity, size, starts_with("hm_")) 
g_homogeneity_outcome = g_homogeneity_outcome %>% mutate(across(homogeneity:hm_var100, ~ as.numeric(.x))) %>%
  mutate(se = apply(g_homogeneity_outcome %>% select(hm_var1:hm_var100), 1, sd)) %>% 
  select(grouping, homogeneity, size, se) 

g_homogeneity_outcome <- g_homogeneity_outcome %>% mutate(upper = homogeneity + 1.96*se, lower = homogeneity - 1.96*se)
max_ = max(g_homogeneity_outcome[g_homogeneity_outcome$size>=5,]$upper, na.rm=TRUE)
min_ = min(g_homogeneity_outcome[g_homogeneity_outcome$size>=5,]$lower, na.rm=TRUE)
g_homogeneity_outcome <- g_homogeneity_outcome %>% arrange(homogeneity)
g_homogeneity_outcome <- g_homogeneity_outcome %>% mutate(homogeneity = as.numeric(homogeneity), 
                                                          size = as.numeric(size), 
                                                          grouping = factor(grouping, levels=g_homogeneity_outcome$grouping))

df_age <- data.frame("main_cat"="Age", "sub_cat"=c("18-35","36-65", "65+"))
df_educ <- data.frame("main_cat" = "Education", "sub_cat"=c("No Diploma" , "High School", "Junior College", "Bachelor's", "Graduate"))
df_race <- data.frame("main_cat"="Race", "sub_cat"=c("White", "Black"))
df_gender <- data.frame("main_cat"="Gender", "sub_cat"=c("Male", "Female"))
df_relig <- data.frame("main_cat"="Religion", "sub_cat"=c("Protestant", "Catholic", "Atheist", "Other Religion", "Jewish"))
#df_attend <- data.frame("main_cat"="Religious Attendance", "sub_cat"=c("Never Attend", "Occasionally Attend", "Frequently Attend"))
df_party <- data.frame("main_cat"="Party Identification", "sub_cat"=c("Democrat", "Republican", "Independent"))
df_class <- data.frame("main_cat"="Social Class", "sub_cat"=c("Lower Class", "Working Class", "Middle Class", "Upper Class"))
df_region <- data.frame("main_cat"="Region", "sub_cat"=c("South", "Not South"))

df_demog <- rbind(df_educ, df_race, df_gender, df_age, df_relig, df_party, df_class, df_region) 
df_demog <- df_demog %>% left_join(g_homogeneity_outcome, by=c("sub_cat" = "grouping")) %>% arrange((homogeneity))
df_demog <- df_demog %>% mutate(sub_cat = paste0(sub_cat, " (n=", size, ")") ) #paste0(sub_cat, " ", size)) paste0(g1, " (n=", g1_size, ")")
df_demog <- df_demog %>% mutate(sub_cat = factor(sub_cat, levels=df_demog$sub_cat), 
                                main_cat = factor(main_cat, levels=c("Education", "Race", "Gender", "Age", "Region", "Party Identification", "Religion", "Social Class")))
min_c = min(df_demog$homogeneity); max_c = max(df_demog$homogeneity)
df_demog = df_demog %>% mutate(normalized_values = (df_demog$homogeneity - min_c) / (max_c - min_c))
index_ = findInterval(df_demog$normalized_values, seq(0, 1, length.out = length(df_demog$normalized_values) ))
index_[which.min(index_)] = index_[which.min(index_)] 
index_[which.max(index_)] = index_[which.max(index_)] 
attribute_colors <- color_gradient(nrow(df_demog)+1)[index_]

df_demog %>% ggplot(aes(x=sub_cat, y=homogeneity)) + 
  geom_point(aes(fill=factor(sub_cat)), shape=24, stroke=0.2, size=3.5, color="black") +  # Set stroke color
  scale_fill_manual(values = attribute_colors) +  # Use scale_fill_manual for fill colors
  ggpubr::theme_pubclean() +
  geom_errorbar( aes(x=sub_cat, ymin=lower, ymax=upper), width = 0.25)+coord_flip()+
  ggforce::facet_col(vars(main_cat), scales = "free_y", space = "free")+
  #facet_wrap(~main_cat, scales="free_y", ncol = 1) + 
  theme(axis.title.y=element_blank(), axis.text.y = element_text(size = 10, color = "black"), legend.position = "none")+
  ylab("Within-Group Consensus") +
  ylim(min_, 0.6)# switch max_ to 0.6 only for manhattan

ggsave(paste0("Panel1_supplement_", year, "_", pairwise_method,".png"),   width = 5, height=8, units="in", dpi= 300)

######################################################################
## Creating the similarity matrix and building the networks
btw_group_similarity_estimation <- function(g1_df, g2_df, pairwise_method){
  # it computes the pairwise similarity between individuals in two groups, and takes an average.
  df_pairs <- expand.grid(unique(g1_df$new_id), unique(g2_df$new_id)) %>% filter(Var1!=Var2) %>% rename("i"="Var1", "j"="Var2")
  df_pairs = pairwise_result_df %>% inner_join(df_pairs, by=c("i", "j"))
  if (pairwise_method == "spearman") {consensus =  mean(df_pairs$sim, na.rm=TRUE)}
  if (pairwise_method == "manhattan") {
    E = 2.962
    consensus = (E - mean(df_pairs$dis, na.rm=TRUE)) / E
  }
  return(consensus)
}

g_mat = gss_demog %>% select(Accountant:Welder) %>% as.matrix()
df_pairs <- expand.grid(1:nrow(g_mat), 1:nrow(g_mat)) %>% filter(Var1!=Var2) %>% mutate(idx = 1:n())
pair_list <- df_pairs %>% split_tibble(c("idx"))
if (pairwise_method == "spearman"){
  pairwise_result_list <- mclapply(pair_list, pairwise_spearman_similarity, mc.cores = 32)
}
if (pairwise_method == "manhattan"){
  pairwise_result_list <- mclapply(pair_list, pairwise_manhattan_distance, mc.cores = 32)
}
pairwise_result_df = as.data.frame(data.table::rbindlist(pairwise_result_list))

btw_group_similarity_df = data.frame()
for(education_1 in c("Graduate", "Bachelor's", "Junior College", "High School", "No Diploma")) for (race_1 in c("White", "Black") ) for (sex_1 in c("Male", "Female")) {
  g1_df = gss_demog %>% filter(sex==sex_1 & degree==education_1 & race==race_1) 
  for(education_2 in c("Graduate", "Bachelor's", "Junior College", "High School", "No Diploma")) for (race_2 in c("White", "Black") ) for (sex_2 in c("Male", "Female")) {
    g2_df = gss_demog %>% filter(sex==sex_2 & degree==education_2 & race==race_2) 
    sim = btw_group_similarity_estimation(g1_df, g2_df, pairwise_method)
    grouping_1 = paste(race_1, sex_1, education_1, sep=", ")
    grouping_2 = paste(race_2, sex_2, education_2, sep=", ")
    g1_size = nrow(g1_df);    g2_size = nrow(g2_df)
    btw_group_similarity_df <- rbind(btw_group_similarity_df, data.frame("g1"=grouping_1, "g2"= grouping_2, "sim"= sim, "g1_size"=g1_size, "g2_size"=g2_size))
  }
}

plot_df <- btw_group_similarity_df %>% filter(g1_size>5 & g2_size>5) 
plot_df <- plot_df %>% mutate(g1_label = paste0(g1, " (n=", g1_size, ")"), g2_label = paste0(g2, " (n=", g2_size, ")")) %>% 
  rename("Similarity" = "sim")
order_df <- plot_df %>% group_by(g1_label) %>% summarise(m_hom = mean(Similarity)) %>% arrange(desc(m_hom))
plot_df <- plot_df %>% mutate(g1_label = factor(g1_label, levels=unique(order_df$g1_label)), g2_label = factor(g2_label, levels=rev(unique(order_df$g1_label))))


plot_df %>% filter(g1 == g2) %>% select(Similarity) %>% summary()
plot_df %>% filter(g1 == g2)  %>% select(Similarity) %>% summarise(sd(Similarity))
plot_df %>% filter(g1 != g2)  %>% select(Similarity) %>% summary() 
plot_df %>% filter(g1 != g2)  %>% select(Similarity) %>% summarise(sd(Similarity))


binarized=TRUE
tenth_percentile <- quantile(plot_df$Similarity, probs = 0.10)
edge_df <- plot_df %>% mutate(weight = Similarity)
edge_df <- edge_df %>% filter(weight > tenth_percentile) 
min_ = min(edge_df$weight); max_ = max(edge_df$weight)
edge_df <- edge_df %>% mutate(weight = 2+10*(weight-min_ )/ (max_ - min_))
edge_df <- edge_df %>% mutate(width = weight^3) # the tie width does not affect layout

vertex_df <- edge_df %>% select(g1, g1_size) %>% distinct() %>% mutate(g1_size = log10(g1_size))
groupings = c(); props = c(); sei10s = c(); print_groupings = c(); income_useds = c(); income_aboves = c()
for(education_ in c("Graduate", "Bachelor's", "Junior College", "High School", "No Diploma")) for (race_ in c("White", "Black", "Other") ) for (sex_ in c("Male", "Female")) {
  g_df = gss_demog %>% filter(sex==sex_ & degree==education_ & race==race_) 
  prop = g_df %>% filter(socialclass == "Middle Class" | socialclass == "Upper Class") %>% nrow() / nrow(g_df)
  sei10 = mean(g_df$prestg10_num, na.rm=TRUE)
  income_used = mean(g_df$income_used, na.rm=TRUE)
  income_above = mean(g_df$income_above, na.rm=TRUE)
  
  grouping = paste(race_, sex_, education_, sep=", ")
  shortened_education_ = ifelse(education_=="Graduate", "MA+", 
                                ifelse(education_=="Bachelor's", "BA", 
                                       ifelse(education_=="Junior College", "JC",
                                              ifelse(education_=="High School", "HS", "LTHS"))))
  print_grouping = paste0(race_, "\n", sex_, "\n", shortened_education_)
  groupings <- c(groupings, grouping)
  print_groupings <- c(print_groupings, print_grouping)
  props <- c(props, prop); sei10s <- c(sei10s, sei10)
  income_aboves <- c(income_aboves, income_above); income_useds <- c(income_useds, income_used)
}
df_attr <- data.frame(g1=groupings, prop = props, sei10=sei10s, income_used=income_useds, income_above=income_aboves,  print_grouping=print_groupings)
df_attr <- df_attr %>% mutate(sei10 = sei10 -30)
df_attr <- df_attr %>% inner_join(grp_size, by = "g1")
vertex_df <- vertex_df %>% inner_join(df_attr, by="g1")
g_homogeneity_outcome = read.csv(paste0("bs_g_consensus_MS", year, "_", pairwise_method,".csv")) %>% select(grouping, homogeneity, size) 
vertex_df <- vertex_df %>% inner_join(g_homogeneity_outcome, by=c("g1"="grouping"))

g <- graph_from_data_frame(edge_df, directed = FALSE, vertices = vertex_df)
g <- igraph::simplify(g, remove.loops=TRUE,  remove.multiple = TRUE, edge.attr.comb = "mean")

# setting the color for binarized and continuous versions
if (binarized == FALSE){
  vertex.frame.color = "lightgray"
  min_ = min(V(g)$homogeneity); max_ = max(V(g)$homogeneity)
  normalized_values <- (V(g)$homogeneity - min_) / (max_ - min_)
  index_ = findInterval(normalized_values, seq(0, 1, length.out = length(normalized_values) + 1))
  index_[which.min(index_)] = index_[which.min(index_)] 
  index_[which.max(index_)] = index_[which.max(index_)] 
  attribute_colors <- color_gradient(length(normalized_values)+1)[index_]
  
  normalized_values <- (E(g)$Similarity - min_) / (max_ - min_)
  index_ = findInterval(normalized_values, seq(0, 1, length.out = length(normalized_values) + 1))
  index_[which.min(index_)] = index_[which.min(index_)] 
  index_[which.max(index_)] = index_[which.max(index_)] 
  edge_colors <- color_gradient(length(normalized_values))[index_]
}
if (binarized == TRUE){
  vertex.frame.color = "white"
  thd <- quantile(plot_df$Similarity, probs = 0.80)
  
  attribute_colors <- ifelse(V(g)$homogeneity > thd, "#0070C0", "#DEEBF7")
  edge_colors <- ifelse(E(g)$Similarity > thd, "#0070C0", "#DEEBF7")
}

E(g)$edge_color = edge_colors
attributes_to_delete = c("Similarity", "g1_size", "g2_size", "g1_label", "g2_label", "width", "edge_color")
g_layout = g
for (attr_name in attributes_to_delete) {g_layout <- delete_edge_attr(g_layout, name = attr_name)}
V(g)$label_color="black"


if (year == 2012 | year == 1989) {
  if (pairwise_method == "manhattan") {set.seed(12345)} # the seed is for preventing label overlap in the net vis.
  if (pairwise_method == "spearman") {set.seed(2)} # the seed is for preventing label overlap in the net vis.
  layout = layout.forceatlas2(g_layout, gravity=100, iterations=1000, plotstep=2000, delta=2., k=10, ks=0.01)
  plot(g, layout=layout,edge.width = 0.02*E(g)$width, vertex.color= attribute_colors, 
       vertex.size=2*sqrt(V(g)$size), vertex.frame.color=vertex.frame.color, vertex.frame.width = 3, 
       vertex.label.cex=1., vertex.label = V(g)$print_grouping, 
       vertex.label.color=V(g)$label_color, label.font = 2)
  position = "topleft"
  if (year==1989 & pairwise_method=="spearman") {position = "topright"}
  if (year==2012 & pairwise_method=="spearman") {position = "bottomleft"}
  if (year==1989 & pairwise_method=="manhattan") {position = "topright"}
  if (year==1989 & pairwise_method=="spearman") {position = "topright"}
  
  # Figure 3a binarized weight of 2012; without label text
  if (binarized==FALSE){
    labeled = FALSE; 
    if (labeled == TRUE){vertex.label = V(g)$print_grouping} else {vertex.label=NA}
    png(paste0("figures/net_",pairwise_method,"_equalsize_", pairwise_method, "_", binarized, "_", labeled, "_", year,".png"), width = 5, height = 5, units = "in", res=300)
    par(mar=c(0,0,0,0))
    plot(g, layout=layout,edge.width = 2, vertex.color= attribute_colors, 
         vertex.size=10, vertex.frame.color=vertex.frame.color, vertex.frame.width = 2, 
         vertex.label.cex=.6, vertex.label = vertex.label, 
         vertex.label.color=V(g)$label_color, vertex.label.font = 2, edge.color=edge_colors)
    filtered_edges <- which(E(g)$Similarity > thd) # 
    g2 <- subgraph.edges(g, filtered_edges, delete.vertices = FALSE)
    plot(g2, layout=layout,edge.width = 2, vertex.color= attribute_colors, 
         vertex.size=10, vertex.frame.color=vertex.frame.color, vertex.frame.width = 2, 
         vertex.label.cex=.6, vertex.label = vertex.label, 
         vertex.label.color=V(g2)$label_color, vertex.label.font = 2, edge.color=E(g2)$edge_color, add=TRUE)
    dev.off()
  }
  
  # Figure 3b continuous weight of 2012; without label text
  if (binarized==TRUE){
    labeled = FALSE
    if (labeled == TRUE){vertex.label = V(g)$print_grouping} else {vertex.label=NA}
    png(paste0("figures/net_",pairwise_method,"_equalsize_", pairwise_method, "_", binarized, "_", labeled, "_", year,".png"), width = 5, height = 5, units = "in", res=300)
    par(mar=c(0,0,0,0))
    plot(g, layout=layout,edge.width = 2, vertex.color= attribute_colors, 
         vertex.size=10, vertex.frame.color=vertex.frame.color, vertex.frame.width = 2, 
         vertex.label.cex=.6, vertex.label = vertex.label, 
         vertex.label.color=V(g)$label_color, vertex.label.font = 2, edge.color=edge_colors)
    filtered_edges <- which(E(g)$Similarity > thd) # 
    g2 <- subgraph.edges(g, filtered_edges, delete.vertices = FALSE)
    plot(g2, layout=layout,edge.width = 2, vertex.color= attribute_colors, 
         vertex.size=10, vertex.frame.color=vertex.frame.color, vertex.frame.width = 2, 
         vertex.label.cex=.6, vertex.label = vertex.label, 
         vertex.label.color=V(g2)$label_color, vertex.label.font = 2, edge.color=E(g2)$edge_color, add=TRUE)
    dev.off()
  }
  
  # Figure 4a: Size
  if (binarized==TRUE) {
    labeled = TRUE
    if (labeled == TRUE){vertex.label = V(g)$print_grouping} else {vertex.label=NA}
    png(paste0("figures/net_",pairwise_method,"_size_", binarized, "_", labeled, "_", year,".png"), width = 5, height = 5, units = "in", res=300)
    par(mar=c(0,0,0,0))
    plot(g, layout=layout,edge.width = 2, vertex.color= attribute_colors, 
         vertex.size=0.05*(V(g)$size), vertex.frame.color=vertex.frame.color, vertex.frame.width = 1.5, 
         vertex.label.cex=.6, vertex.label = vertex.label, 
         vertex.label.color=V(g)$label_color, vertex.label.font = 2, edge.color=edge_colors)
    filtered_edges <- which(E(g)$Similarity > thd)
    g2 <- subgraph.edges(g, filtered_edges, delete.vertices = FALSE)
    
    plot(g2, layout=layout,edge.width = 2, vertex.color= attribute_colors, 
         vertex.size=0.25*(V(g)$size), vertex.frame.color=vertex.frame.color, vertex.frame.width = 1.5, 
         vertex.label.cex=0.6, vertex.label = vertex.label, 
         vertex.label.color=V(g2)$label_color, vertex.label.font = 2, edge.color=E(g2)$edge_color, add=TRUE)
    dev.off()
  }
  # Figure 4b: Income
  if (binarized==TRUE) {
    labeled = TRUE
    if (labeled == TRUE){vertex.label = V(g)$print_grouping} else {vertex.label=NA}
    png(paste0("figures/net_",pairwise_method,"_income_", binarized, "_", labeled, "_", year,".png"), width = 5, height = 5, units = "in", res=300)
    par(mar=c(0,0,0,0))
    plot(g, layout=layout,edge.width = 2, vertex.color= attribute_colors, 
         vertex.size=110*(V(g)$income_used+0.1)^2, vertex.frame.color=vertex.frame.color, vertex.frame.width = 1.5, 
         vertex.label.cex=.6, vertex.label = vertex.label, 
         vertex.label.color=V(g)$label_color, vertex.label.font = 2, edge.color=edge_colors)
    filtered_edges <- which(E(g)$Similarity > thd)
    g2 <- subgraph.edges(g, filtered_edges, delete.vertices = FALSE)
    
    plot(g2, layout=layout,edge.width = 2, vertex.color= attribute_colors, 
         vertex.size=110*(V(g)$income_used+0.1)^2, vertex.frame.color=vertex.frame.color, vertex.frame.width = 1.5, 
         vertex.label.cex=0.6, vertex.label = vertex.label, 
         vertex.label.color=V(g2)$label_color, vertex.label.font = 2, edge.color=E(g2)$edge_color, add=TRUE)
    dev.off()
  }
  # Figure 4c: SEI
  if (binarized==TRUE) {
    labeled = TRUE
    if (labeled == TRUE){vertex.label = V(g)$print_grouping} else {vertex.label=NA}
    png(paste0("figures/net_",pairwise_method,"_SEI_", binarized, "_", labeled, "_", year,".png"), width = 5, height = 5, units = "in", res=300)
    par(mar=c(0,0,0,0))
    plot(g, layout=layout,edge.width = 2, vertex.color= attribute_colors, 
         vertex.size=0.04*(V(g)$sei10+1)^2, vertex.frame.color=vertex.frame.color, vertex.frame.width = 1.5, 
         vertex.label.cex=.6, vertex.label = vertex.label, 
         vertex.label.color=V(g)$label_color, vertex.label.font = 2, edge.color=edge_colors)
    filtered_edges <- which(E(g)$Similarity > thd)
    g2 <- subgraph.edges(g, filtered_edges, delete.vertices = FALSE)
    
    plot(g2, layout=layout,edge.width = 2, vertex.color= attribute_colors, 
         vertex.size=0.04*(V(g)$sei10+1)^2, vertex.frame.color=vertex.frame.color, vertex.frame.width = 1.5, 
         vertex.label.cex=0.6, vertex.label = vertex.label, 
         vertex.label.color=V(g2)$label_color, vertex.label.font = 2, edge.color=E(g2)$edge_color, add=TRUE)
    dev.off()
  }
  # Figure 4d: upper-middle class
  if (binarized==TRUE) {
    labeled = TRUE
    if (labeled == TRUE){vertex.label = V(g)$print_grouping} else {vertex.label=NA}
    png(paste0("figures/net_",pairwise_method,"_prop_", binarized, "_", labeled, "_", year,".png"), width = 5, height = 5, units = "in", res=300)
    par(mar=c(0,0,0,0))
    plot(g, layout=layout,edge.width = 2, vertex.color= attribute_colors, 
         vertex.size=30*(V(g)$prop), vertex.frame.color=vertex.frame.color, vertex.frame.width = 1.5, 
         vertex.label.cex=.6, vertex.label = vertex.label, 
         vertex.label.color=V(g)$label_color, vertex.label.font = 2, edge.color=edge_colors)
    filtered_edges <- which(E(g)$Similarity > thd)
    g2 <- subgraph.edges(g, filtered_edges, delete.vertices = FALSE)
    
    plot(g2, layout=layout,edge.width = 2, vertex.color= attribute_colors, 
         vertex.size=30*(V(g)$prop), vertex.frame.color=vertex.frame.color, vertex.frame.width = 1.5, 
         vertex.label.cex=0.6, vertex.label = vertex.label, 
         vertex.label.color=V(g2)$label_color, vertex.label.font = 2, edge.color=E(g2)$edge_color, add=TRUE)
    dev.off()
  }
}

