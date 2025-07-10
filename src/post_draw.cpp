#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector post_draw_gmix_cpp(const NumericVector &prior_shape,
                                 const NumericVector &prior_rate,
                                 const NumericVector &Qn,
                                 const IntegerMatrix &N, const NumericMatrix &E,
                                 int nsim) {
  int I = N.nrow();
  int J = N.ncol();
  int m = prior_shape.size();

  NumericVector lambda_sim(Dimension(I, J, nsim)); // flat array
  IntegerVector m_range = Rcpp::seq(0, m - 1);

  for (int i = 0; i < I; ++i) {
    for (int j = 0; j < J; ++j) {
      // 1. Get qij vector
      NumericVector qij(m);
      for (int k = 0; k < m; ++k) {
        qij[k] = Qn[i + I * j + I * J * k];
      }

      // generate mixture component index
      IntegerVector mix_idx_draws =
          sample(m_range, nsim, /*replace=*/true, /*probs=*/qij);

      // 3. Generate gamma samples
      for (int s = 0; s < nsim; ++s) {
        int comp_idx = mix_idx_draws[s];
        // now generate gamma
        double shape = N(i, j) + prior_shape[comp_idx];
        double rate = E(i, j) + prior_rate[comp_idx];
        lambda_sim[i + I * j + I * J * s] =
            R::rgamma(shape, 1 / rate); // Note: R::rgamma uses scale = 1/rate
      }
    }
  }

  return lambda_sim; // interpreted as array [nsim, I, J] in R
}

// [[Rcpp::export]]
NumericVector post_draw_discrete_cpp(const NumericVector &grid,
                                     const NumericVector &esti_prior,
                                     const IntegerMatrix &N,
                                     const NumericMatrix &E, int nsim) {
  int I = N.nrow();
  int J = N.ncol();
  int K = grid.size();

  NumericVector lambda_sim(Dimension(I, J, nsim)); // flat 3D array
  NumericVector log_esti_prior = log(esti_prior);

  for (int i = 0; i < I; ++i) {
    for (int j = 0; j < J; ++j) {
      int nij = N(i, j);
      double eij = E(i, j);

      // Step 1: compute unnormalized log posterior probs
      NumericVector log_post_prob(K);
      for (int k = 0; k < K; ++k) {
        double mean = grid[k] * eij;
        log_post_prob[k] = nij * std::log(mean) - mean + log_esti_prior[k];
      }
      double max_log_post_prob = Rcpp::max(log_post_prob);
      NumericVector post_prob = exp(log_post_prob - max_log_post_prob);

      NumericVector draws = Rcpp::sample(grid, nsim, /*replace=*/true,
                                         /*probs=*/post_prob);

      // Step 3: sample from grid with cumulative weights cum_prob
      for (int s = 0; s < nsim; ++s) {
        lambda_sim[i + I * j + I * J * s] = draws[s];
      }
    }
  }

  return lambda_sim; // shape: [I * J * nsim]
}
