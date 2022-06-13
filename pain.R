# Read dependencies
source("_dependencies/dependencies.R")

# Set sheet url
url <- "https://docs.google.com/spreadsheets/d/15pNDO6b-sun56CY3SUcRl-Sfhk4FsU9biyb9J4WfYDw/"

# Prepare dataframe to sent information up to the sheet
messenger <- NULL

################################################################################
#                                   Pain conditions                            #
################################################################################

# Get the data from the gsheet
pc.all <- read_sheet(
  url, sheet = "Pain <> SWB (Conditions)"
)

# Get chronic pain conditions
pc.chronic <- pc.all %>% filter(`condition (general)` == "chronic") %>% 
  filter(!is.na(d))
# Prepare for meta
pc.chronic <- pc.chronic %>% rowwise() %>% mutate(
  n1 = `sample size`/2,
  n2 = `sample size`/2,
  d_se = getDSE(d, n1=n1, n2=n2)
) %>% ungroup()
# meta
pc.chronic.mod <- rma.mv(yi = d, V = d_se^2,
                slab = citation,
                random =  ~1 | citation,
                test = "t", method="REML",
                data = pc.chronic)

# prepare message for the sheet
messenger <- rbind(messenger,
                   data.frame(value=coef(pc.chronic.mod)[[1]], 
                              ci=paste0(round.c(pc.chronic.mod$ci.lb,2), ", ", 
                                        round.c(pc.chronic.mod$ci.ub,2)), 
                              nstudies=uniqueN(pc.chronic$citation), 
                              neffects=nrow(pc.chronic),
                              nobs=sum(pc.chronic$`sample size`), 
                              note="average effect on swb (SD) of having (non-cancer) chronic pain")
)

# Get cancer pain conditions
pc.cancer <- pc.all %>% filter(`condition (general)` == "terminal/cancer")
# Prepare for meta
pc.cancer <- pc.cancer %>% rowwise() %>% mutate(
  n1 = `sample size`/2,
  n2 = `sample size`/2,
  d_se = getDSE(d, n1=n1, n2=n2)
) %>% ungroup()
# meta
pc.cancer.mod <- rma.mv(yi = d, V = d_se^2,
                         slab = citation,
                         random =  ~1 | citation,
                         test = "t", method="REML",
                         data = pc.cancer)

# prepare message for the sheet
messenger <- rbind(messenger,
                   data.frame(value=coef(pc.cancer.mod)[[1]], 
                              ci=paste0(round.c(pc.cancer.mod$ci.lb,2), ", ", 
                                        round.c(pc.cancer.mod$ci.ub,2)), 
                              nstudies=uniqueN(pc.cancer$citation), 
                              neffects=nrow(pc.cancer),
                              nobs=sum(pc.cancer$`sample size`), 
                              note="average effect on swb (SD) of having cancer pain")
)

# Get migraine conditions
pc.migraine <- pc.all %>% filter(grepl("migraine", tolower(`condition (specific)`)))

# Prepare for meta
pc.migraine <- pc.migraine %>% rowwise() %>% mutate(
  n1 = `sample size`/2,
  n2 = `sample size`/2,
  d_se = getDSE(d, n1=n1, n2=n2)
) %>% ungroup()
# meta
pc.migraine.mod <- rma.mv(yi = d, V = d_se^2,
                        slab = citation,
                        random =  ~1 | citation,
                        test = "t", method="REML",
                        data = pc.migraine)

# prepare message for the sheet
messenger <- rbind(messenger,
                   data.frame(value=coef(pc.migraine.mod)[[1]], 
                              ci=paste0(round.c(pc.migraine.mod$ci.lb,2), ", ", 
                                        round.c(pc.migraine.mod$ci.ub,2)), 
                              nstudies=uniqueN(pc.migraine$citation), 
                              neffects=nrow(pc.migraine),
                              nobs=sum(pc.migraine$`sample size`), 
                              note="average effect on swb (SD) of having migraine pain")
)

################################################################################
#                                    Treatments                                #
################################################################################

# Get the data from the gsheet
treat.all <- read_sheet(
  url, sheet = "Pain <> SWB (Treatments)"
)

###################
# Psychotherapies #
###################

# Get the treatments
treat.psyc <- treat.all %>% filter(grepl("psych", tolower(treatment)))
# filter
treat.psyc <- treat.psyc %>% filter(!is.na(d)) # only if we have effect size
treat.psyc <- treat.psyc %>% filter(!is.na(d_se)) # need d_se, 
# for now this excludes the data we extracted from McCracken et al. 2022

# All sorts of FU, comparisons (active, TAU, etc.), samples (adults, children),
# are included for now.

# Convert the effect sizes
treat.psyc <- treat.psyc %>% mutate(d = ifelse(isHighScoreGood, d*-1, d))

# Note whether the outcome is pain or swb
treat.psyc <- treat.psyc %>% mutate(
  measure_general = ifelse(grepl("pain", measure), "pain", "swb"),
  measure_id = paste0(citation, "_", measure)
)
treat.psyc.pain <- treat.psyc %>% filter(measure_general == "pain")
treat.psyc.swb  <- treat.psyc %>% filter(measure_general == "swb")

# meta
treat.psyc.pain.mod <- rma.mv(yi = d, V = d_se^2,
                              slab = citation,
                              random =  ~1 | citation,
                              test = "t", method="REML",
                              data = treat.psyc.pain)

treat.psyc.swb.mod <- rma.mv(yi = d, V = d_se^2,
                             slab = citation,
                             random =  ~1 | citation,
                             test = "t", method="REML",
                             data = treat.psyc.swb)

# prepare message for the sheet
messenger <- rbind(messenger,
                   data.frame(value=coef(treat.psyc.pain.mod)[[1]], 
                              ci=paste0(round.c(treat.psyc.pain.mod$ci.lb,2), ", ", 
                                        round.c(treat.psyc.pain.mod$ci.ub,2)), 
                              nstudies=uniqueN(treat.psyc.pain$citation), 
                              neffects=nrow(treat.psyc.pain),
                              nobs=sum(treat.psyc.pain$`sample size`), 
                              note="average effect on pain (SD) of psychology-based treatments for chronic pain")
)

messenger <- rbind(messenger,
                   data.frame(value=coef(treat.psyc.swb.mod)[[1]], 
                              ci=paste0(round.c(treat.psyc.swb.mod$ci.lb,2), ", ", 
                                        round.c(treat.psyc.swb.mod$ci.ub,2)), 
                              nstudies=uniqueN(treat.psyc.swb$citation), 
                              neffects=nrow(treat.psyc.swb),
                              nobs=sum(treat.psyc.swb$`sample size`), 
                              note="average effect on swb (SD) of psychology-based treatments for chronic pain")
)

# Relationship between pain and swb
# x <- treat.psyc %>% pivot_wider(names_from = measure_general, values_from = c(`sample size`, d, d_se))

################################################################################
#                     Send information back up to the sheet                    #
################################################################################

write_sheet(messenger, url, "Values from R")
