# Ran in R version 4.2.0 (2022-04-22), on aarch64-apple-darwin20 (64-bit) platform under macOS Monterey 12.6.

# I used the following versions of packages employed: dplyr_1.0.9, tidyverse_1.3.1, brms_2.17.0,
# tidybayes_3.0.2, ggplot2_3.3.6 and patchwork_1.1.1.

# set working directory (works only in RStudio)
setwd( dirname(rstudioapi::getSourceEditorContext()$path) )

# list required packages into a character object
pkgs <- c(
  "dplyr", "tidyverse", # for data wrangling
  "brms", "tidybayes", # for Bayesian analyses 
  "ggplot2", "ggpubr" # for plotting
)

# load or install packages as needed
for ( i in pkgs ) {
  if ( i %in% rownames( installed.packages() ) == F ) install.packages(i) # install if it ain't installed yet
  if ( i %in% names( sessionInfo()$otherPkgs ) == F ) library( i , character.only = T ) # load if it ain't loaded yet
}

# set some values for later
options( mc.cores = parallel::detectCores() ) # use all parallel CPU cores
ch = 4 # number of chains
it = 1500 # iterations per chain
wu = 500 # warm-up iterations, to be discarded
ad = .90 # adapt_delta parameter
s = 87542 # seed for reproducibility

# Note that although I set a seed for all models, the results are only exactly
# reproducible on the same operating system with the same C++ compiler and version.

# create folders "models", "figures", "tables" and "sessions" to store results and sessions info in
# prints TRUE and creates the folder if it was not present, prints NULL if the folder was already present.
sapply( c("models", "figures", "tables", "sessions"), function(i) if( !dir.exists(i) ) dir.create(i) )

# set ggplot theme
theme_set( theme_classic(base_size = 14) )

# prepare colors to use in graphs (a colorblind-friendly palette)
cbPal <- c( "#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7" )

# read data
d0 <- read.csv( "data/20220207_dbs_iadl_data_full.csv" , sep = "," ) # the data set
d.par <- read.table( "data/20220224_stim_pars.txt" , header = T ) # stimulation parameters
d.nms <- read.csv( "data/item_nms.csv" , sep = "," , row.names = 1 ) # item names


# ----------- heuristic causal model (DAG)  -----------

# will be created via DiagrammeR
# for now, I know I need to adjust for preop. DRS2, BDI-II and LEDD for total effect
# and adjust for postop. DRS-2, BDI-II and LEDD for direct effect


# ----------- pre-processing  -----------

# keep only those variables I will need (PDAQ, ID, DRS-2, BDI-II, LEDD, and demographics)
d1 <- d0 %>% select( id, assessment_type, levodopa_equivalent, drsii_total, bdi, # predictors
                     time_after_surg, age_neuropsy, sex, edu_type, years_edu, # demographics
                     starts_with("pdaq") # outcome
                     ) %>%
  # rename the variables so that I don't go crazy
  rename( ass_type = assessment_type, led = levodopa_equivalent, drs = drsii_total ) %>%
  # change "post_1y" to "post" because they are all one-year after DBS operation
  mutate( ass_type = ifelse( ass_type == "post_1y", "1_post", "0_pre" ) %>% as.factor() )

# check there ain't no NAs in the data set
sapply( names(d1) , function(i) sum( is.na(d1[[i]]) ) ) # seems ok
all( sapply( names(d1) , function(i) sum( is.na(d1[[i]]) ) ) == 0 ) # indeed, it is ok, no NAs

# extract the number of patients
N = length( unique(d1$id) ) # N = 32

# find a reasonable scaling of predictors
# first get an idea about in-sample variability
sapply( levels(d1$ass_type) , function(i) # loop through pre/post assessments
  sapply( c("drs","bdi","led") , function(j) # loop through all predictors
    round( sd( d1[[j]][ d1$ass_type == i ] , na.rm = T ) , 2 )
    )
  )

# will center all predictors of interest on their pre-surgery mean
# extract the means from data and save them to a table
scl <- data.frame(
  # first create a column of means to use for scaling
  M = c(
    mean( d1$drs[d1$ass_type == "0_pre"] , na.rm = T ),
    mean( d1$bdi[d1$ass_type == "0_pre"] , na.rm = T ),
    mean( d1$led[d1$ass_type == "0_pre"] , na.rm = T )
  ),
  # then add user(me)-specified scaling factors based on theoretical consideration and data variability
  SD = c(
    3, # 3 point for DRS-2, half-way between max (144/144) and likely MCI (138/144), close to SDs pre and post
    7, # 7 points for BDI/II, equivalent to feeling better/worse in a third of the items (7/21), close to SDs pre and post
    500 # 500 mg for LEDD, equivalent to 2 Isicom 250 or 5 Isicom 100 pills, close to SD post
  ),
  # add names to rows so that I know which row belongs to which variable
  row.names = c( "drs" , "bdi" , "led" )
)

# transform the data from wide to long format over PDAQ items
d2 <- d1 %>% pivot_longer( cols = paste0( "pdaq_" , 1:15 ) , names_to = "item" , values_to = "resp" ) %>%
  mutate(
  # change names of PDAQ items to integers
  # the gsub part removes everything before the underscore (text)
  # so that only integer remains and is easy to transform via as.integer()
  # could have also gone via as.integer( as.factor() )
  item = as.integer( gsub( "^.*\\_" , "" , item ) )
)

# calculate number of PDAQ items (K = 15)
length( unique(d2$item) )
max( d2$item ) # A-OK

# check that patient/item/assessment type each combination is there exactly once
table( d2$id , d2$item , d2$ass_type ) # A-OK

# scale DRS and LEDD such that the model samples more efficiently
# make the outcome ordinal such that brms can take it in
df <- d2 %>%
  mutate( drs = ( drs - scl["drs","M"] ) / scl["drs","SD"],
          bdi = ( bdi - scl["bdi","M"] ) / scl["bdi","SD"],
          led = ( led - scl["led","M"] ) / scl["led","SD"],
          resp = as.ordered( resp ), item = paste0("i",item)
          )

# check the contrast for pre/post factor
contrasts( df$ass_type ) # ok, it's a dummy coding as I wanted


# ---- models set-up ---- 

# setp-up linear models
f <- list(
  orig = bf( resp ~ ass_type * led + ass_type * drs + ass_type * bdi + ( 1 | id ) + ( 1 | item ) ) + cumulative("logit"),
  prob = bf( resp ~ ass_type * led + ass_type * drs + ass_type * bdi + ( 1 | id ) + ( 1 | item ) ) + cumulative("probit"),
  ledd = bf( resp ~ ass_type * led + ( 1 | id ) + ( 1 | item ) ) + cumulative("probit")#,
  #flex = bf( resp | thres(gr = item) ~ ass_type * led + ( 1 | id ) + ( 1 | item ) ) + cumulative("probit")
)

# set-up priors (the same for all models)
p <- c( prior( normal( 0 , 0.5 ) , class = b ),
        prior( student_t( 3 , 0 , 2.5 ) , class = Intercept ),
        prior( student_t( 3 , 0 , 2.5 ) , class = sd , group = id ),
        prior( student_t( 3 , 0 , 2.5 ) , class = sd , group = item )
        )


# ---- models fitting ----

# prepare a list for single models
m <- list()

# conduct model fitting
for ( i in names(f) ) m[[i]] <- brm( formula = f[[i]], prior = p,
                                     data = df, sample_prior = T, seed = s, chains = ch,
                                     iter = it, warmup = wu, control = list( adapt_delta = ad ),
                                     file = paste0( "models/",i,".rds" ),
                                     save_model = paste0("models/",i,".stan")
                                     )
