# PPC function with separate p-values for each data source
ppcOcc_intMsPGOcc <- function(object, fit.stat = "chi-squared",  thin.by = 1) {
  
  # Input validation
  if (missing(object)) {
    stop("error: object must be specified")
  }
  if (class(object) != "intMsPGOcc") {
    stop("error: object must be of class 'intMsPGOcc'")
  }
  if (!tolower(fit.stat) %in% c("chi-squared", "freeman-tukey", "chi-square")) {
    stop("error: fit.stat must be either 'chi-squared' or 'freeman-tukey'")
  }
  fit.stat <- tolower(fit.stat)

  # Helper functions
  logit_inv <- function(z, a = 0, b = 1) {
    b - (b - a) / (1 + exp(z))
  }
  
  cat("Running posterior predictive check for intMsPGOcc output...\n")
  
  # Extract model components
  y <- object$y                    
  sites <- object$sites           
  species_list <- object$species  
  X_p <- object$X.p               
  
  n_data <- length(y)
  n_samples_total <- object$n.post * object$n.chains
  
  # Apply thinning
  sample_indices <- seq(1, n_samples_total, by = thin.by)
  n_samples <- length(sample_indices)
  
  cat("Using", n_samples, "samples (every", thin.by, "of", n_samples_total, "total samples)\n")
  cat("Calculating p-values for", n_data, "data sources\n")
  
  # Extract thinned samples
  z_samples <- object$z.samples[sample_indices, , , drop = FALSE]
  alpha_samples <- object$alpha.samples[sample_indices, , drop = FALSE]
  
  # Initialize storage - separate for each data source
  chi_sq_obs_by_source <- matrix(0, nrow = n_samples, ncol = n_data)
  chi_sq_rep_by_source <- matrix(0, nrow = n_samples, ncol = n_data)
  chi_sq_obs_total <- numeric(n_samples)
  chi_sq_rep_total <- numeric(n_samples)
  epsilon <- 1e-04
  
  # Process each MCMC sample
  for(s in 1:n_samples) {
    
    if(s %% 50 == 0) cat("Processing sample", s, "of", n_samples, "\n")
    
    sample_chi_obs_total <- 0
    sample_chi_rep_total <- 0
    
    # Process each data source
    for(d in 1:n_data) {
      
      sample_chi_obs_d <- 0  # Chi-square for this data source
      sample_chi_rep_d <- 0  # Chi-square replicated for this data source
      
      y_curr <- y[[d]]              
      sites_curr <- sites[[d]]     
      species_curr <- species_list[[d]]  
      X_p_curr <- X_p[[d]]         
      
      # Map species to full community indices
      if(is.character(species_curr)) {
        sp_indices <- match(species_curr, object$sp.names)
      } else {
        sp_indices <- species_curr
      }
      
      n_species_d <- length(sp_indices)
      
      # Get alpha coefficients for this data source
      alpha_start_idx <- 1
      if(d > 1) {
        for(prev_d in 1:(d-1)) {
          alpha_start_idx <- alpha_start_idx + 
            length(species_list[[prev_d]]) * ncol(X_p[[prev_d]])
        }
      }
      n_alpha_curr <- n_species_d * ncol(X_p_curr)
      alpha_curr <- alpha_samples[s, alpha_start_idx:(alpha_start_idx + n_alpha_curr - 1)]
      alpha_matrix <- matrix(alpha_curr, nrow = n_species_d, ncol = ncol(X_p_curr), byrow = TRUE)
      
      # Process each species in this data source
      for(i in 1:n_species_d) {
        
        sp_idx <- sp_indices[i]
        alpha_sp <- alpha_matrix[i, ]
        
        # Collect all valid observations and expected values for this species
        obs_vals <- c()
        exp_vals <- c()
        rep_vals <- c()
        site_rep_ids <- c()  # Track which site/rep combination each value comes from
        
        # Loop through sites and replicates
        for(j in 1:length(sites_curr)) {
          site_idx <- sites_curr[j]
          z_val <- z_samples[s, sp_idx, site_idx]
          
          for(k in 1:dim(y_curr)[3]) {
            if(!is.na(y_curr[i, j, k])) {
              
              # Observed value
              obs_vals <- c(obs_vals, y_curr[i, j, k])
              
              # Expected value
              obs_row_idx <- (j - 1) * dim(y_curr)[3] + k
              if(obs_row_idx <= nrow(X_p_curr)) {
                logit_p <- sum(X_p_curr[obs_row_idx, ] * alpha_sp)
                p_detect <- logit_inv(logit_p) * z_val
              } else {
                next
              }
              exp_vals <- c(exp_vals, p_detect)
              
              # Replicated value
              rep_vals <- c(rep_vals, rbinom(1, 1, p_detect))
                site_rep_ids <- c(site_rep_ids, k)  # Group by replicate
              
            }
          }
        }
        
        # Group the values and calculate chi-square contributions
        if(length(obs_vals) > 0) {
          
          unique_groups <- unique(site_rep_ids)
          
          for(group_id in unique_groups) {
            group_mask <- site_rep_ids == group_id
            
            obs_sum <- sum(obs_vals[group_mask])
            exp_sum <- sum(exp_vals[group_mask])
            rep_sum <- sum(rep_vals[group_mask])
            
            if(exp_sum > 0) {
              if(fit.stat %in% c("chi-squared", "chi-square")) {
                contrib_obs <- (obs_sum - exp_sum)^2 / (exp_sum + epsilon)
                contrib_rep <- (rep_sum - exp_sum)^2 / (exp_sum + epsilon)
              } else if(fit.stat == "freeman-tukey") {
                contrib_obs <- (sqrt(obs_sum) - sqrt(exp_sum))^2
                contrib_rep <- (sqrt(rep_sum) - sqrt(exp_sum))^2
              }
              
              # Add to data source specific totals
              sample_chi_obs_d <- sample_chi_obs_d + contrib_obs
              sample_chi_rep_d <- sample_chi_rep_d + contrib_rep
            }
          }
        }
      }
      
      # Store data source specific chi-square values
      chi_sq_obs_by_source[s, d] <- sample_chi_obs_d
      chi_sq_rep_by_source[s, d] <- sample_chi_rep_d
      
      # Add to overall totals
      sample_chi_obs_total <- sample_chi_obs_total + sample_chi_obs_d
      sample_chi_rep_total <- sample_chi_rep_total + sample_chi_rep_d
    }
    
    # Store overall totals
    chi_sq_obs_total[s] <- sample_chi_obs_total
    chi_sq_rep_total[s] <- sample_chi_rep_total
  }
  
  # Calculate overall Bayesian p-value
  bayesian_pvalue_overall <- mean(chi_sq_rep_total >= chi_sq_obs_total)
  
  # Calculate data source specific p-values
  bayesian_pvalues_by_source <- numeric(n_data)
  for(d in 1:n_data) {
    bayesian_pvalues_by_source[d] <- mean(chi_sq_rep_by_source[, d] >= chi_sq_obs_by_source[, d])
  }
  
  # Create data source names if not available
  source_names <- c("eDNA", "seine")
  names(bayesian_pvalues_by_source) <- source_names
  
  cat("Posterior predictive check complete!\n")
  cat("Overall Bayesian p-value:", round(bayesian_pvalue_overall, 4), "\n")
  cat("Data source specific p-values:\n")
  for(d in 1:n_data) {
    cat("  ", source_names[d], ":", round(bayesian_pvalues_by_source[d], 4), "\n")
  }
  
  # Return comprehensive results
  results <- list(
    # Overall results
    bayesian.pvalue = bayesian_pvalue_overall,
    chi_sq_obs_total = chi_sq_obs_total,
    chi_sq_rep_total = chi_sq_rep_total,
    mean_chi_sq_obs = mean(chi_sq_obs_total),
    mean_chi_sq_rep = mean(chi_sq_rep_total),
    
    # Data source specific results
    bayesian.pvalues.by.source = bayesian_pvalues_by_source,
    chi_sq_obs_by_source = chi_sq_obs_by_source,
    chi_sq_rep_by_source = chi_sq_rep_by_source,
    mean_chi_sq_obs_by_source = apply(chi_sq_obs_by_source, 2, mean),
    mean_chi_sq_rep_by_source = apply(chi_sq_rep_by_source, 2, mean),
    
    # Model info
    n_samples = n_samples,
    n_data_sources = n_data,
    group = group,
    fit.stat = fit.stat,
    source_names = source_names
  )
  
  return(results)
}

ppcOut <- ppcOcc_intMsPGOcc(out.1, "chi-squared", thin.by=100)
