data {
  int<lower=0> N;  // Total number of observations
  int<lower=0, upper=1> injury[N];  // Injury status (0 = no, 1 = yes)
  int<lower=1> K;  // Number of blocks
  int<lower=1, upper=K> block_id[N]; // Block ID 
  int<lower=1> G;  // Total number of groups 
  int<lower=1, upper=G> group[N]; // Group 
  int<lower=1> X;  // Total number of sexes
  int<lower=1, upper=X> sex[N]; // Sex 
}

parameters {
  vector[K] alpha_block_male_raw;   
  vector[K] alpha_block_female_raw;
  real mu_alpha_male;
  real mu_alpha_female;
  real<lower=0> sigma_alpha_male;
  real<lower=0> sigma_alpha_female;
  matrix[G, X] beta_group_sex_raw;
  real mu_beta;
  real<lower=0> sigma_beta;
  vector[G] beta_group;       
  vector[X] beta_sex;      
}

transformed parameters {
  vector[K] alpha_block_male;
  vector[K] alpha_block_female;
  matrix[G, X] beta_group_sex;
  alpha_block_male = mu_alpha_male + sigma_alpha_male * alpha_block_male_raw;
  alpha_block_female = mu_alpha_female + sigma_alpha_female * alpha_block_female_raw;
  beta_group_sex = mu_beta + sigma_beta * beta_group_sex_raw;
}

model {
  // Priors 
  alpha_block_male_raw ~ normal(0, 0.4);
  alpha_block_female_raw ~ normal(0, 0.4);
  mu_alpha_male ~ normal(0, 1);
  mu_alpha_female ~ normal(0, 1);
  sigma_alpha_male ~ exponential(1);
  sigma_alpha_female ~ exponential(1);
  
  for (g in 1:G) {
    for (x in 1:X) {
      beta_group_sex_raw[g, x] ~ normal(0, 0.4);
    }
  }
  
  mu_beta ~ normal(0, 1);
  sigma_beta ~ exponential(1);
  beta_group ~ normal(0, 0.4);
  beta_sex ~ normal(0, 0.4);

  // Likelihood
  //for (i in 1:N) {
    //real alpha_block_effect;
    //if (sex[i] == 2) {  
      //alpha_block_effect = alpha_block_male[block_id[i]];
    //} else {  
      //alpha_block_effect = alpha_block_female[block_id[i]];
    //}
    //real eta = alpha_block_effect +
              // beta_group[group[i]] +
              // beta_sex[sex[i]] +
              // beta_group_sex[group[i], sex[i]]; 
    //injury[i] ~ bernoulli_logit(eta);
  //}
}

generated quantities {
  vector[N] log_lik;  // Log likelihood for each observation
  int<lower=0, upper=1> injury_sim[N];  // Simulated injury outcomes
  matrix[G, 2] group_level_probs;  // Group-level probabilities for each sex

  // Calculate log likelihood and simulate injury outcomes
  for (i in 1:N) {
    real alpha_block_effect;
    if (sex[i] == 2) {  
      alpha_block_effect = alpha_block_male[block_id[i]];
    } else {  
      alpha_block_effect = alpha_block_female[block_id[i]];
    }
    real eta = alpha_block_effect +
               beta_group[group[i]] +
               beta_sex[sex[i]] +
               beta_group_sex[group[i], sex[i]]; 
    log_lik[i] = bernoulli_logit_lpmf(injury[i] | eta);
    injury_sim[i] = bernoulli_logit_rng(eta);  // Simulate new outcome based on eta
  }

  // Calculate group-level probabilities separately for males and females
  for (g in 1:G) {
    real avg_alpha_block_male = mean(alpha_block_male);  // Average alpha_block across all male blocks
    real avg_alpha_block_female = mean(alpha_block_female);  // Average alpha_block across all female blocks
    
    real eta_male = avg_alpha_block_male +
                    beta_group[g] +
                    beta_sex[2] +  
                    beta_group_sex[g, 2];
    real eta_female = avg_alpha_block_female +
                      beta_group[g] +
                      beta_sex[1] +  
                      beta_group_sex[g, 1]; 
                      
    group_level_probs[g, 2] = inv_logit(eta_male);
    group_level_probs[g, 1] = inv_logit(eta_female);
  }
}
