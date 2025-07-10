#include <Rcpp.h>
#include <RcppEigen.h>
#include <Eigen/Core>
#include <vector>
#include <algorithm>
#include <cmath>
#include <limits>
// [[Rcpp::depends(RcppEigen)]]
#include <Rmath.h> // For R::digamma
//#include <boost/math/special_functions/digamma.hpp>  //boost::math::digamma(x)
#include <iostream>



struct NBmixResult {
  Eigen::MatrixXd r;
  Eigen::MatrixXd P;
  Eigen::MatrixXd h;
  int iter;
  double loglik;
};


Eigen::MatrixXd tcrossprod(const Eigen::VectorXd& a, const Eigen::VectorXd& b) {
  return a * b.transpose();  // Compute the outer product
}



Eigen::ArrayXXd mx_dnbinom_log(const Eigen::ArrayXXd& TN1, const Eigen::ArrayXXd& TR1,
                        const Eigen::ArrayXXd& TEh, const Eigen::VectorXd& P,
                        const int& len_N) {
  Eigen::ArrayXXd probs = (TEh + 1.0).inverse();
  Eigen::ArrayXXd log_probs = probs.log();
  Eigen::ArrayXXd log_1_minus_probs = (1.0 - probs).log();
  Eigen::ArrayXXd TN1_plus_TR1 = TN1 + TR1;
  Eigen::ArrayXXd lgamma_TN1_plus_TR1 = TN1_plus_TR1.unaryExpr([](double x) { return std::lgamma(x); });
  Eigen::ArrayXXd lgamma_TR1 = TR1.unaryExpr([](double x) { return std::lgamma(x); });
  Eigen::ArrayXXd lgamma_TN1_plus_1 = (TN1 + 1.0).unaryExpr([](double x) { return std::lgamma(x); });
  Eigen::ArrayXXd log_nb_lik = lgamma_TN1_plus_TR1
    - lgamma_TR1
    - lgamma_TN1_plus_1
    + (TR1 * log_probs)
    + (TN1 * log_1_minus_probs);
  Eigen::ArrayXXd r2 = (Eigen::VectorXd::Ones(len_N) * P.transpose()).array().log(); // Outer product
  return (log_nb_lik + r2);
}

Eigen::VectorXd vectorIndexing(const Eigen::VectorXd& v, const std::vector<int>& indices) {
  Eigen::VectorXd result(indices.size());
  for (std::vector<int>::size_type i = 0; i < indices.size(); ++i) {
    result(i) = v(indices[i]);
  }
  return result;
}

Eigen::ArrayXXd array_colIndexing(const Eigen::ArrayXXd& A, const std::vector<int>& indices) {
  Eigen::ArrayXXd reduced_A(A.rows(), indices.size());
  for (std::vector<int>::size_type i = 0; i < indices.size(); ++i) {
    reduced_A.col(i) = A.col(indices[i]);
  }
  return reduced_A;
}


// Core function
NBmixResult NBmix_s1_EM_g(const Eigen::VectorXd& N,
                          const Eigen::VectorXd& E,
                          Eigen::VectorXd& h,
                          Eigen::VectorXd& grid,
                          double alpha,
                          int maxi,
                          double eps,
                          bool dirichlet) {

  int k = grid.size();
  int len_N = N.size();
  int iter = 0;
  if(dirichlet == false){
    alpha = 1;
  }

  if (h.size() == 1) {
    h = Eigen::VectorXd::Constant(k, h[0]);
  }


  Eigen::VectorXd r = grid.array() / h.array() + 1.0;
  Eigen::VectorXd P = Eigen::VectorXd::Ones(k) / k;

  double loglik = -std::numeric_limits<double>::infinity();
  double loglik_new;


  Eigen::ArrayXXd TN1 = tcrossprod(N, Eigen::VectorXd::Ones(k)).array();
  Eigen::ArrayXXd TR1 = tcrossprod(Eigen::VectorXd::Ones(len_N), r).array();
  Eigen::ArrayXXd Delta = Eigen::MatrixXd::Ones(len_N, k).array();
  Eigen::ArrayXXd Tau = Eigen::MatrixXd::Ones(len_N, k).array();
  Eigen::ArrayXXd tmp = Eigen::MatrixXd::Ones(len_N, k).array();
  Eigen::ArrayXd tmp_max;
  Eigen::ArrayXd r_nume;
  Eigen::ArrayXd r_deno;
  Eigen::ArrayXd h0_nume;
  Eigen::ArrayXd h0_deno;
  Eigen::ArrayXXd TEh = tcrossprod(E, h).array();
  Eigen::ArrayXXd T1h = tcrossprod(Eigen::VectorXd::Ones(len_N), h).array();
  Eigen::ArrayXXd TE1 = tcrossprod(E, Eigen::VectorXd::Ones(k)).array();
  std::vector<int> filtered_seq;
  Eigen::VectorXd h0;
  Eigen::Array<bool, Eigen::Dynamic, 1> indip;
  double h_dif = 1.0;
  double h_tol = 1e-10;

  Eigen::VectorXd v;

  bool stop = false;

  while(!stop){
    int k = P.size();
    //E-step
    // updating Delta
    Delta = TR1 * ( (TR1 + TN1).unaryExpr([](double x) { return R::digamma(x); })
      - TR1.unaryExpr([](double x) { return R::digamma(x); }));
    Delta = (TR1 < 1e-200).select(1, Delta);
    //updating Tau
    if(iter==0){
      tmp = mx_dnbinom_log(TN1, TR1, TEh, P, len_N);
      Tau = (tmp.colwise() - tmp.rowwise().maxCoeff()).exp();
      Tau = Tau.colwise() / Tau.rowwise().sum();
    }else{
      Tau = (tmp.colwise() - tmp_max).exp();
      Tau = Tau.colwise() / Tau.rowwise().sum();
    }



    //CM-steps
    P = ((Tau.colwise().sum() + alpha - 1.0) / (len_N + k*alpha - k)).matrix();




    indip = (P.array() < 0);
    P = indip.select(0.0, P).matrix();
    P = P / P.sum();

    if(dirichlet == true){
      filtered_seq = std::vector<int>();
      for (auto i = 0; i < indip.size(); ++i) {
        if (!indip(i)) {
          filtered_seq.push_back(i);
        }
      }
      P = vectorIndexing(P, filtered_seq);
      r = vectorIndexing(r, filtered_seq);
      h = vectorIndexing(h, filtered_seq);

      Tau = array_colIndexing(Tau, filtered_seq);
      Tau = Tau.colwise() / Tau.rowwise().sum();

      Delta = array_colIndexing(Delta, filtered_seq);
      TEh = array_colIndexing(TEh, filtered_seq);
      TE1 = array_colIndexing(TE1, filtered_seq);
      TN1 = array_colIndexing(TN1, filtered_seq);
      TR1 = array_colIndexing(TR1, filtered_seq);
      T1h = array_colIndexing(T1h, filtered_seq);

    }

    //updating r_k
    r_nume = -(Tau*Delta).colwise().sum();
    r_deno = (-Tau*(1.0 + TEh).log()).colwise().sum();
    r = r_nume / r_deno;

    r = r.array().max(1.0);



    TR1 = tcrossprod(Eigen::VectorXd::Ones(len_N), r).array();

    //updating h_k

    while(h_dif>h_tol){
      h0_nume = (TN1*Tau).colwise().sum();
      h0_deno = ((TN1+TR1)*Tau/(1/TE1+T1h)).colwise().sum();
      h0 = h0_nume / h0_deno;

      h_dif = (h - h0).array().abs().maxCoeff();
      h = h0;
      T1h = tcrossprod(Eigen::VectorXd::Ones(len_N), h).array();
    }


    //check parameters
    v = r.array() * h.array().square();
    for (int i = 0; i < v.size(); ++i) {
      if (v(i) > 1.0) {
        h(i) = std::sqrt(1.0 / r(i));
      }
    }

    h = h.array().max(1e-10);
    r = r.array().min(1e+8);
    h = h.array().min(1);
    TR1 = tcrossprod(Eigen::VectorXd::Ones(len_N), r).array();
    T1h = tcrossprod(Eigen::VectorXd::Ones(len_N), h).array();
    TEh = tcrossprod(E, h).array();

    //compute log-likelihood
    tmp = mx_dnbinom_log(TN1, TR1, TEh, P, len_N);
    tmp_max = tmp.rowwise().maxCoeff();
    loglik_new = ((tmp.colwise() - tmp_max).exp().rowwise().sum().log() + tmp_max).sum();


    if ((std::abs(loglik - loglik_new) <= eps) && (iter > 100)) {
      stop = true;
    }
    loglik = loglik_new;
    iter++;
    if(iter>= maxi){
      stop = true;
    }


  }

  NBmixResult result;
  result.r = r;
  result.P = P;
  result.h = h;
  result.iter = iter;
  result.loglik = loglik;
  return result;
}



// [[Rcpp::export]]
Rcpp::List NBmix_s1_EM_g_Rcpp(const Eigen::VectorXd& N,
                              const Eigen::VectorXd& E,
                              Eigen::VectorXd& h,
                              Eigen::VectorXd& grid,
                              double alpha,
                              int maxi,
                              double eps,
                              bool dirichlet) {

  // Call the core function using Rcpp::as to convert R types to std::vector<int>
  NBmixResult result = NBmix_s1_EM_g(N, E, h, grid, alpha, maxi, eps, dirichlet);

  // Return results as Rcpp List
  return Rcpp::List::create(
    Rcpp::Named("r") = result.r,
    Rcpp::Named("omega") = result.P,
    Rcpp::Named("h") = result.h,
    Rcpp::Named("iter") = result.iter,
    Rcpp::Named("loglik") = result.loglik
  );
}





