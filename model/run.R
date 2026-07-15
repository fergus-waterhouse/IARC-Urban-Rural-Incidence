"Run mortality model.

Usage:
		run.R <inc> <out_dir> <predictor> <cancer> <sex> <num_iter> <num_burn> [--sigma=<sigma>] [--knots=<knots>] [--num_chains=<num_chains>]
		run.R (-h | --help)

Options:
		-h --help                       Show this help message and exit.
		--num_chains=<num_chains>       Number of parallel chains [default: 4].
		--knots=<knots>                 Number of knots for splines [default: 2].
		--sigma=<sigma>                 The standard deviation of the country-level variation normal prior [default: 1.0].

Arguments:
    <inc>           Incidence dataset.
    <out_dir>       Output directory.
    <predictor>     Name of predictor variable.
    <cancer>        Cancer site of the run.
    <sex>           Sex of the run (1 for male, 2 for female, 0 for both).
    <num_iter>      Number of iterations for the MCMC.
    <num_burn>      Number of burn in iterations for the MCMC.
" -> doc

args <- docopt::docopt(doc)


# ----------------- SETUP -----------------
suppressPackageStartupMessages({
    library(dplyr)
    library(nimble)
    library(tidyverse)
    library(lme4)
    library(glue)
    library(parallel)
    require(coda)
    require(dplyr)
    require(tibble)
})

source("model/utils.R")
source("model/modelspline.R")
source("model/modelspline_sex.R")

dir.create(args$out_dir)

# ----------------- SELECT MODEL LOGIC -----------------

if (args$sex != 0) {
    selected_model_function <- run_MCMC_spline_sex
} else {
    selected_model_function <- run_MCMC_spline
}

if (args$num_chains < 2) {
    stop("[ERROR] INVALID NUMBER OF CHAINS (MUST BE > 1)")
}


# -- STEP 1 --------------- LOAD INCIDENCE DATA -----------------
inc <- read.csv(args$inc) 

if (args$sex != 0) {
    inc <- filter(inc, sex == args$sex)
}

inc <- inc %>% 
    filter(cancer_lab == args$cancer) %>%
    mutate(cov = .data[[args$predictor]]) %>%
    mutate(
        # Tier 1: Continent
        hier1.id = as.integer(as.factor(continent)),
        # Tier 2: Region
        hier2.id = as.integer(as.factor(region)),
        # Tier 3: Country
        hier3.id = as.integer(as.factor(country)),
    )


# ----------------- PRINT INFORMATION -----------------
cat(glue("

    ###### URBANRURALINC ############################################
    
    ====== RUNNING URBANRURALINC MODEL ======
    Cancer Site: {args$cancer}
    Sex:         {args$sex}
    Predictor:   {args$predictor}
    Sigma:       {args$sigma}
    
    --- DATA SUMMARY ---
    Total Rows:          {format(nrow(inc), big.mark=',')}
    Total Cases:         {format(sum(inc$cases), big.mark=',')}
    Total Person-Years:  {format(sum(inc$py), big.mark=',')}
    Age Groups:          {length(unique(inc$agr))}
    Hier 1 (Continents): {length(unique(inc$hier1.id))}
    Hier 2 (Regions):    {length(unique(inc$hier2.id))}
    Hier 3 (Countries):  {length(unique(inc$hier3.id))}
    Registries:          {length(unique(inc$registry))}
    
    --- MCMC SETTINGS ---
    Chains:       {args$num_chains}
    Iterations:   {args$num_iter}
    Burn-in:      {args$num_burn}
    Thin (Mort):  100
    Thin (Param): 100
    
    --- DIRECTORY INFO ---
    Loading Incidence: {args$inc}
    Saving Output To:  {args$out_dir}
    
    --- EQUATION ---
    ln(rate) = beta0 + beta1*age + bs3 %*% Z + [gamma_sex] + deltas3 * cov
"))


meta_info <- data.frame(
    "Meta Info" = c(
        "Incidence Data", 
        "Sex (0: Both, 1: Men, 2: Women)", 
        "Cancer", "Predictor", 
        "Num. of knots", 
        "Sigma", 
        "Num. Iterations", 
        "Num. Burn-In", 
        "Num. Chains"
        ),
    "Value" = c(
        args$inc, 
        args$sex, 
        args$cancer, 
        args$predictor, 
        args$knots, 
        args$sigma, 
        args$num_iter, 
        args$num_burn, 
        args$num_chains
        )
)

write.csv(
    meta_info, 
    file = str_c(args$out_dir, "/meta_info.csv")
    )

# -- STEP 2 --------------- DEFINE FIXED HIER2 & HIER3 TERMS -----------------

fixedhier2 <- detect_fixed_nodes(
    data = inc, 
    group_col = "hier1.id", 
    target_col = "hier2.id", 
    return_col = "hier2.id"
)

fixedhier3 <- detect_fixed_nodes(
    data = inc, 
    group_col = "hier3.id", 
    target_col = "cov", 
    return_col = "hier3.id"
)


# -- STEP 3 --------------- CALCULATE INITIAL VALUES -----------------

cat("\n\n\n--------- INITIATING MODEL ---------\n")
initial <- init_model(inc)
cat(" Model Successfully Initiated")


# -- STEP 4 --------------- RUN MODEL PREPROCESSING -----------------
model_inputs            <- grid_lookup(inc)
model_inputs$var.hier2  <- as.numeric(which(!fixedhier2))
model_inputs$fix.hier2  <- as.numeric(which(fixedhier2))
model_inputs$var.hier3  <- as.numeric(which(!fixedhier3))
model_inputs$fix.hier3  <- as.numeric(which(fixedhier3))
model_inputs$spline.k   <- as.integer(args$knots)
model_inputs$spline.z   <- thin_plate(1:18, as.integer(args$knots))
model_inputs$sigma      <- as.numeric(args$sigma)


# -- STEP 5 --------------- SET UP CLUSTER AND RUN -----------------
cat("\n\n\n--------- RUNNING MCMC -------------\n")

this_cluster <- makeCluster(as.numeric(args$num_chains))
invisible(clusterEvalQ(this_cluster, {
    library(nimble)
}))

cat(paste0(c(" [1] Running MCMC Chains... (Started: ", format(Sys.time(), "%H:%M:%S"), ")")))
chain_output <- parLapply(
    cl         = this_cluster,
    X          = 1:as.numeric(args$num_chains),
    fun        = selected_model_function,
    model_name = "incidence_model",
    incidence  = inc,
    init_vals  = initial,
    n_iter     = as.numeric(args$num_iter),
    n_burn     = as.numeric(args$num_burn),
    inputs     = model_inputs
)
stopCluster(this_cluster)
cat(paste0(c(" (Ended: ", format(Sys.time(), "%H:%M:%S"), ")")))

# -- STEP 6 --------------- SAVE DATA -----------------

cat(sprintf("\n [2] Saving Data to %s\n", args$out_dir))

saveRDS(
    chain_output,
    file = str_c(args$out_dir, "/mcmc_out.rds")
)


# -- STEP 7 --------------- GENERATE REPORT -----------------
cat("\n\n--------- GENERATING REPORTS -------------")

report_path <- str_c(args$out_dir, "/mcmc_summary_report.csv")
waic_path   <- str_c(args$out_dir, "/mcmc_waic_report.csv")

generate_mcmc_report(
    mcmc_raw = chain_output, 
    out_file = report_path,
    waic_file = waic_path
)

cat("\n\n###### URBANRURALINC RUN COMPLETE ###################################\n\n")
