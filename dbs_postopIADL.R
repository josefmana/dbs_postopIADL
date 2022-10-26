# Ran in R version 4.2.0 (2022-04-22), on aarch64-apple-darwin20 (64-bit) platform under macOS Monterey 12.6.

# I used the following versions of packages employed: dplyr_1.0.9, tidyverse_1.3.1, brms_2.17.0,
# tidybayes_3.0.2, ggplot2_3.3.6 and patchwork_1.1.1.

# set working directory (works only in RStudio)
setwd( dirname(rstudioapi::getSourceEditorContext()$path) )

# list required packages into a character object
pkgs <- c(
  "dplyr", "tidyverse", # for data wrangling
  "brms", "tidybayes", "bayestestR", # for Bayesian analyses 
  "ggplot2", "patchwork" # for plotting
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
theme_set( theme_classic(base_size = 14 ) )

# prepare colors to use in graphs (a colorblind-friendly palette)
cbPal <- c( "#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7" )

# read data
d0 <- read.csv( "data/20220207_dbs_iadl_data_full.csv" , sep = "," ) # the data set
d.par <- read.table( "data/20220224_stim_pars.txt" , header = T ) # stimulation parameters
d.nms <- read.csv( "data/item_nms.csv" , sep = "," , row.names = 1 ) # item names


# ----------- heuristic causal model (DAG)  -----------

# will be created via ggdag
# for now, I know I need to adjust for preop. DRS2, BDI-II and LEDD for total effect
# and adjust for postop. DRS-2, BDI-II and LEDD for direct effect


# ----------- pre-processing  -----------

# keep only those variables I will need (PDAQ, ID, DRS-2, BDI-II, LEDD, and demographics)
d1 <- d0 %>% select( id, assessment_type, levodopa_equivalent, drsii_total, bdi, # predictors
                     time_after_surg, age_neuropsy, sex, edu_type, years_edu, # demographics
                     starts_with("pdaq") # outcome
                     ) %>%
  # rename the variables so that I don't go crazy
  rename( post = assessment_type, led = levodopa_equivalent, drs = drsii_total ) %>%
  mutate( post = ifelse( post == "post_1y", 1, 0 ) %>% as.factor() ) # re-code post-surgery indicator variable

# check there ain't no NAs in the data set
sapply( names(d1) , function(i) sum( is.na(d1[[i]]) ) ) # seems ok
all( sapply( names(d1) , function(i) sum( is.na(d1[[i]]) ) ) == 0 ) # indeed, it is ok, no NAs

# find a reasonable scaling of predictors
# first get an idea about in-sample variability
sapply( levels(d1$post) , function(i) # loop through pre/post assessments
  sapply( c("drs","bdi","led") , function(j) # loop through all predictors
    round( sd( d1[[j]][ d1$post == i ] , na.rm = T ) , 2 )
    )
  )

# will center all predictors of interest on their pre-surgery mean
# extract the means from data and save them to a table
scl <- data.frame(
  # first create a column of means to use for scaling
  M = c(
    mean( d1$drs[d1$post == 0] , na.rm = T ),
    mean( d1$bdi[d1$post == 0] , na.rm = T ),
    mean( d1$led[d1$post == 0] , na.rm = T )
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
table( d2$id , d2$item , d2$post ) # A-OK

# scale DRS and LEDD such that the model samples more efficiently
# make the outcome ordinal such that brms can take it in
df <- d2 %>%
  mutate( drs = ( drs - scl["drs","M"] ) / scl["drs","SD"],
          bdi = ( bdi - scl["bdi","M"] ) / scl["bdi","SD"],
          led = ( led - scl["led","M"] ) / scl["led","SD"],
          resp = as.ordered( resp )
          )

# check the number of patients
length( unique(df$id) ) # N = 32

# check contrasts for the session
contrasts( df$post ) # indicator variable coding, A-OK


# ---- models set-up ---- 

# set-up linear models
# m0 = total effect of time ("main effect") (the masked model)
# m1 = total effect of LEDD ("total effect") (the primary model)
# m2 = direct effect of LEDD ("direct effect") (the control model)
# m3 = direct effect of LEDD ("original model") (the original logit model presented in the article)
f <- list(
  m0_meff = bf( resp ~ 1 + post + ( 1 | id ) + ( 1 | item ) ) + cumulative("probit"),
  m1_teff = bf( resp ~ 1 + post * led + ( 1 | id ) + ( 1 | item ) ) + cumulative("probit"),
  m2_deff = bf( resp ~ 1 + post * led + post * drs + post * bdi + ( 1 | id ) + ( 1 | item ) ) + cumulative("probit"),
  m3_orig = bf( resp ~ 1 + post * led + post * drs + post * bdi + ( 1 | id ) + ( 1 | item ) ) + cumulative("logit")
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


# ----------- soft model checking -----------

# check the highest Rhat for chains convergence and Pareto-k for influential outliers
cbind.data.frame(
  R_hat = sapply( names(m) , function(i) max( rhat(m[[i]]), na.rm = T ) ) %>% round(3),
  Pareto_k = sapply( names(m) , function(i) max( loo(m[[i]])$diagnostics$pareto_k, na.rm = T ) ) %>% round(3)
)

# check Pareto-k visually as well
par( mfrow = c(2,2) )
for ( i in names(m) ) plot( loo(m[[i]]), main = i )
par( mfrow = c(1,1) )

# plot posterior predictive check stratified by item
ppc <- list()

# plotting proper
for ( i in names(m) ) {
  
  # set distinct color palettes for different models
  if( i == "m0_meff") bayesplot::color_scheme_set( scheme = "red" )
  else if ( i == "m1_teff" ) bayesplot::color_scheme_set( scheme = "orange" )
  else if ( i == "m2_deff" ) bayesplot::color_scheme_set( scheme = "yellow" )
  else if ( i == "m3_orig" ) bayesplot::color_scheme_set( scheme = "green" )
  
  # plot the per item posterior predictive check
  ppc[[i]] <- pp_check( m[[i]], type = "bars_grouped", ndraws = NULL, group = "item" ) +
    theme( legend.position = "none" )
  
}

# plot them together
( ( ppc$m0_meff + labs(title = "Ordered-probit (main effect)") ) | ( ppc$m1_teff + labs( title = "Ordered-probit (total effect)" ) ) ) /
  ( ( ppc$m2_deff + labs( title = "Ordered-probit (direct effect)" ) ) | ( ppc$m3_orig + labs( title = "Ordered-logit (direct effect)" ) ) )

# save it
ggsave( "figures/per item posterior predictive checks.jpg", width = 1.5*9.62, height = 1.2*11.7, dpi = 600 )


# ----------- posterior predictive check (item/session) -----------

# in this section I extract ppc for m1 only
# extract posterior retrodictions in the form of predicted response for each data row in each MCMC iteration
ppred <- predict( m$m1_teff, summary = F )

# prepare a table containing each item/session combination
d_seq <- expand.grid( unique(df$item), unique(df$post) ) %>% `colnames<-`( c("item","post") ) %>%
  # add labels of both items and sessions
  mutate( Item = factor( rep( d.nms$name , 2 ), levels = d.nms$name , ordered = T ),
          Time = factor( ifelse( post == 0 , "pre-surgery", "post-surgery"),
                         levels = paste0(c("pre","post"),"-surgery"), ordered = T )
          )

# add observed frequencies
for ( i in 0:4 ) {
  # prepare a column for each response type
  d_seq[[paste0("resp_",i)]] <- NA
  # loop through all item/session combinations to fill-in observed frequencies of each response 0:4
  for( j in 1:nrow(d_seq) ) d_seq[ j, paste0("resp_",i) ] <- ( with( df, ( resp[ which( item == d_seq$item[j] & post == d_seq$post[j] ) ] ) ) == i ) %>% sum()
}

# change the table to a long format
d_seq <- d_seq %>% pivot_longer(
  cols = contains("resp_"), values_to = "obs", names_to = "resp",
  names_transform = function(x) sub("^[^_]*_","",x) %>% as.numeric()
) %>%
  # add columns for median prediction, and lower and upper predictive interval boundaries
  mutate( pred = NA, ci.low = NA, ci.upp = NA ) %>%
  as.data.frame()

# add predictions for each item/session/response combination
for ( i in 1:nrow(d_seq) ) d_seq[ i , c("pred","ci.low","ci.upp") ] <- c(
  sapply( 1:nrow(ppred), function(j) ( ppred[ j , with( df, which( item == d_seq$item[i] & post == d_seq$post[i] ) ) ] == d_seq$resp[i]+1 ) %>% sum() ) %>% median(),
  sapply( 1:nrow(ppred), function(j) ( ppred[ j , with( df, which( item == d_seq$item[i] & post == d_seq$post[i] ) ) ] == d_seq$resp[i]+1 ) %>% sum() ) %>% ci(.95,"ETI") %>% as.data.frame() %>% select(CI_low) %>% as.numeric(),
  sapply( 1:nrow(ppred), function(j) ( ppred[ j , with( df, which( item == d_seq$item[i] & post == d_seq$post[i] ) ) ] == d_seq$resp[i]+1 ) %>% sum() ) %>% ci(.95,"ETI") %>% as.data.frame() %>% select(CI_high) %>% as.numeric()
)

# plot it
d_seq %>%
  # shift such that zero-response frequencies are still visualized
  mutate( obs = obs+1 ) %>%
  # plotting proper
  ggplot( aes( x = resp , y = obs , fill = Time ) ) +
  # plot observed frequencies as bar plots
  geom_bar( stat = "identity", position = position_dodge( width = .8 ), width = .7, alpha = .4 ) +
  # plot predicted frequencies as points with whiskers
  geom_pointrange(
    with( d_seq, aes( x = resp, y = pred+1, ymin = ci.low+1, ymax = ci.upp+1, color = Time, fill = Time ) ),
    shape = 21, size = 1.2, fatten = 1.5, position = position_dodge( width = .8 )
    ) +
  # wrap it
  facet_wrap( ~ Item  , scales = "free" , nrow = 5 , ncol = 3 ) +
  labs( x = "Response", y = "Frequency", fill = "", color = "" ) +
  scale_fill_manual( values = cbPal[c(1,6)] ) +
  scale_color_manual( values = cbPal[c(1,6)] ) +
  scale_y_continuous( limits = c(0, 32), breaks = seq(1,31,5) , labels = seq(0,30,5) ) +
  theme( legend.position = "bottom" )

# save it
ggsave( "figures/m1_teff_ppc_item_session.jpg", dpi = 600 )


# ----------- posterior predictive check (composite score per session) -----------

# prepare a table containing each item/session combination
d_seq <- data.frame( post = 0:1 ) %>% mutate(
  Time = factor( ifelse( post == 0 , "pre-surgery", "post-surgery"), levels = paste0(c("pre","post"),"-surgery"), ordered = T ),
  `0` = c( ( df$resp[ which(df$post == 0) ] == 0 ) %>% sum(), ( df$resp[ which(df$post == 1) ] == 0 ) %>% sum() ),
  `1` = c( ( df$resp[ which(df$post == 0) ] == 1 ) %>% sum(), ( df$resp[ which(df$post == 1) ] == 1 ) %>% sum() ),
  `2` = c( ( df$resp[ which(df$post == 0) ] == 2 ) %>% sum(), ( df$resp[ which(df$post == 1) ] == 2 ) %>% sum() ),
  `3` = c( ( df$resp[ which(df$post == 0) ] == 3 ) %>% sum(), ( df$resp[ which(df$post == 1) ] == 3 ) %>% sum() ),
  `4` = c( ( df$resp[ which(df$post == 0) ] == 4 ) %>% sum(), ( df$resp[ which(df$post == 1) ] == 4 ) %>% sum() )
) %>%
  pivot_longer( cols = 3:7, values_to = "obs", names_to = "resp" ) %>%
  # add columns for median prediction, and lower and upper predictive interval boundaries
  mutate( resp = as.numeric(resp), pred = NA, ci.low = NA, ci.upp = NA ) %>% as.data.frame()
                                               
# add predictions for each item/session/response combination
for ( i in 1:nrow(d_seq) ) d_seq[ i , c("pred","ci.low","ci.upp") ] <- c(
  sapply( 1:nrow(ppred), function(j) ( ppred[ j , with( df, which( post == d_seq$post[i] ) ) ] == d_seq$resp[i]+1 ) %>% sum() ) %>% median(),
  sapply( 1:nrow(ppred), function(j) ( ppred[ j , with( df, which( post == d_seq$post[i] ) ) ] == d_seq$resp[i]+1 ) %>% sum() ) %>% ci(.95,"ETI") %>% as.data.frame() %>% select(CI_low) %>% as.numeric(),
  sapply( 1:nrow(ppred), function(j) ( ppred[ j , with( df, which( post == d_seq$post[i] ) ) ] == d_seq$resp[i]+1 ) %>% sum() ) %>% ci(.95,"ETI") %>% as.data.frame() %>% select(CI_high) %>% as.numeric()
)

# plot it
d_seq %>% mutate( obs = obs+1 ) %>%
  ggplot( aes( x = resp , y = obs , fill = Time ) ) +
  geom_bar( stat = "identity", position = position_dodge( width = .8 ), width = .7, alpha = .4 ) +
  geom_pointrange(
    with( d_seq, aes( x = resp, y = pred+1, ymin = ci.low+1, ymax = ci.upp+1, color = Time, fill = Time ) ),
    shape = 21, size = 2, fatten = 2.5, position = position_dodge( width = .8 )
  ) +
  facet_wrap( ~ Time  , scales = "free" , nrow = 2 , ncol = 1 ) +
  labs( x = "Response", y = "Frequency", fill = "", color = "" ) +
  scale_fill_manual( values = cbPal[c(1,6)] ) +
  scale_color_manual( values = cbPal[c(1,6)] ) +
  scale_y_continuous( limits = c(0, 322), breaks = seq(1,301,100) , labels = seq(0,300,100) ) +
  theme( legend.position = "none" )

# save it
ggsave( "figures/m1_teff_ppc_comp_session.jpg", width = .75 * 9.65, height = .75 * 11.7, dpi = 600 )


# ---- parameters description ----

# plot the population-level parameters of the total effect of LEDD (i.e., m1_teff)
m$m1_teff %>%
  # extract parameter estimates
  spread_draws( b_post1, b_led, `b_post1:led` ) %>%
  pivot_longer( cols = contains("b_"), values_to = "est", names_to = "par" ) %>%
  mutate( par = factor( par, levels = rev( paste0("b_",c("post1","led","post1:led") ) ), ordered = T ) ) %>%
  # plot the values
  ggplot( aes( y = par, x = est, fill = stat(abs(x) > .1) ) ) +
  stat_halfeye( .width = .95, fatten_point = 3 ) +
  geom_vline( xintercept = 0, linetype = "dashed" ) +
  labs( x = "Effect size", y = "Parameter", title = "Regression of IADL on time and LEDD") +
  scale_y_discrete( labels = c("Post x LEDD", "LEDD", "Post-surgery") ) +
  scale_x_continuous( breaks = seq(-.2,1,.2), labels = seq(-.2,1,.2) ) +
  scale_fill_manual( values = c("gray80","skyblue") ) +
  theme( legend.position = "none", plot.title = element_text(hjust=.5) )

# save it
ggsave( "figures/m1_teff_par_estimates.jpg", dpi = 600 )

# extract posteriors of the b_post parameter from all three types of models and plot them
sapply( names(m)[1:3], function(i) as_draws_df( m[[i]] ) %>% select(b_post1) ) %>%
  cbind.data.frame() %>%
  rename( "Vanilla" = "m0_meff.b_post1", "LEDD adjusted" = "m1_teff.b_post1", "Direct effect" = "m2_deff.b_post1") %>%
  pivot_longer( everything(), names_to = "Model:", values_to = "Effect size" ) %>%
  mutate( `Model:` = factor( `Model:`, levels = c("Vanilla","LEDD adjusted","Direct effect", ordered = T) ) ) %>%
  # plotting proper
  ggplot( aes(x = `Effect size`, color = `Model:`, fill = `Model:` ) ) +
  stat_slab( geom = "slab", linetype = "solid", size = 2 ) +
  geom_vline( xintercept = 0, linetype = "dashed" ) +
  labs( x = "Effect size", y = "Density", title = "Post-minus-pre-surgery IADL (standardized effect)" ) +
  scale_y_continuous( breaks = seq(0,1,.2), labels = seq(0,1,.2) %>% sprintf("%.1f",.) ) +
  scale_x_continuous( breaks = seq(-.2,1.2,.2), labels = seq(-.2,1.2,.2) %>% sprintf("%.1f",.) ) +
  scale_color_manual( values = alpha( cbPal[c(5,3,4)], 1 ) ) +
  scale_fill_manual( values = alpha( cbPal[c(5,3,4)], .1 ) ) +
  theme( legend.position = "bottom", plot.title = element_text(hjust=.5) )

# save it
ggsave( "figures/b_post1_estimates.jpg", dpi = 600 )


# ---- conditional effects ----

# prepare a list holding the conditional effects tabs
ce <- list()

# extract main effects
for ( i in c("post","led") ) ce[[i]] <- conditional_effects( m$m1_teff , categorical = T , effects = i )

# prepare conditions for interaction extraction¨
cond <- make_conditions( m$m1_teff, vars = "post" )

# add interaction to the list
ce$`post:led` <- conditional_effects( m$m1_teff, categorical = T, effects = "led", conditions = cond )

# add response levels names
for ( i in names(ce) ) {
  ce[[i]][[paste0(sub(".*:","",i),":cats__")]] <- ce[[i]][[paste0(sub(".*:","",i),":cats__")]] %>%
    mutate( `Response:` = recode(cats__ , "4" = '"none"', "3" = '"a little"' , "2" = '"somewhat"' , "1" = '"a lot"' , "0" = '"cannot do"' ) )
}

# prepare the pre-vs-post figure
ce$A <- ce$post[["post:cats__"]] %>%
  ggplot( aes(x = post, y = estimate__, ymin = lower__, ymax = upper__, color = `Response:`) ) +
  geom_point( position = position_dodge(.66), size = 5 ) +
  geom_errorbar( position = position_dodge(.66), size = 1 , width = .5) +
  scale_x_discrete( name = "Time of assessment", labels = c("Pre-surgery","Post-surgery") ) +
  scale_y_continuous( name = "Probability of response (%)",
                      limits = c(0,1) , breaks = seq(0,1,.1), labels = seq(0,100,10)
                      ) +
  theme_classic( base_size = 14 ) + theme( legend.position = "none" )

# prepare the LEDD (main effect) figure
ce$B <- ce$led[["led:cats__"]] %>%
  ggplot( aes(x = led, y = estimate__, ymin = lower__, ymax = upper__, color = `Response:`, fill = `Response:`) ) +
  geom_line( size = 1.5 ) +
  geom_ribbon( alpha = .1 , linetype = 0 ) +
  scale_x_continuous( name = "LEDD (mg)", labels = seq(500,4000,500),
                      breaks = (seq(500,4000,500)-scl["led","M"]) / scl["led","SD"]
                      ) +
  scale_y_continuous( name = "Probability of response (%)",
                      limits = c(0,1) , breaks = seq(0,1,.1) , labels = seq(0,100,10)
                      ) +
  theme_classic( base_size = 14 ) + theme( legend.position = "none", axis.title.y = element_blank() )

# prepare the LEDD:time interaction figure
ce$C <- ce$`post:led`[["led:cats__"]] %>%
  mutate( post = factor( ifelse(post == 0, "pre-surgery", "post-surgery"),
                         levels = paste0(c("pre","post"),"-surgery"), ordered = T )
          ) %>%
  ggplot( aes(x = led, y = estimate__, ymin = lower__, ymax = upper__, color = `Response:`, fill = `Response:`) ) +
  geom_line( size = 1.5 ) +
  geom_ribbon( alpha = .1 , linetype = 0 ) +
  scale_x_continuous( name = "LEDD (mg)", labels = seq(500,4000,500),
                      breaks = (seq(500,4000,500)-scl["led","M"]) / scl["led","SD"]
                      ) +
  scale_y_continuous( name = "Probability of response (%)",
                      limits = c(0,1) , breaks = seq(0,1,.1) , labels = seq(0,100,10)
                      ) +
  facet_wrap( ~ post ) +
  theme_classic( base_size = 14 ) + theme( legend.position = "bottom" )

# collect the plots into a single figure
( ce$A | ce$B ) / ce$C + plot_annotation( tag_levels = "A" )

# save it
ggsave( "figures/m1_teff_conditional_effects.jpg", width = 9.64, height = 1.5 * 6, dpi = 600 )


# ---- save information needed for prediction calculations ----

# prepare a folder
if( !dir.exists("shiny") ) dir.create("shiny")

# remove the original data from m1_teff model such that it's shareable
shiny.m <- m$m1_teff
shiny.m$data <- data.frame(resp = factor(NA, levels = c(0,1)), post = NA, led = NA, id = NA, item = NA)

# save the information needed
saveRDS( list(m = shiny.m, scl = scl), "shiny/stat_model.RDS" )


# ----------- session info -----------

# write the sessionInfo() into a .txt file
capture.output( sessionInfo(), file = "sessions/dbs_postopIADL_modelling.txt" )
