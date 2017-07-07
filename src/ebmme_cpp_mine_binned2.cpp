#include "utils.h"
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;


// [[Rcpp::export]]
RcppExport SEXP newtissue(SEXP cat_, SEXP block_id_,
                        SEXP nclust_,
                        SEXP bins_,
                        SEXP fs_binned_, 
                        SEXP alpha_,
                        SEXP inner_iters_)
{
    mat bins = as<mat>(bins_);
    int m = bins.n_rows;
    int inner_iters = as<int>(inner_iters_); 
 
    vec cat = as<vec>(cat_);
    vec cat_id = unique(cat);
    int ncat = cat_id.n_elem; 
    int nclust = as<int>(nclust_);
    
    mat p = zeros<mat>(nclust,m);
    vec block_id  = as<vec>(block_id_);
    vec ublock_id = unique(block_id);
    int B = ublock_id.n_elem;
    
    cube fs = zeros<cube>(m, B, nclust);
    
    mat f = zeros<mat>(m, nclust);
    cube fs_binned = as<cube>(fs_binned_);
    
    
    vec alpha = as<vec>(alpha_);
    mat a = zeros<mat>(nclust,ncat);
    
    mat w_sum = zeros<mat>(nclust,ncat);
    
    for (int i = 0; i < m; i++){
        for (int k = 0; k < nclust; k++){
            for (int b=0; b < B; b++){
                fs(i, b, k) = fs_binned(bins(i, b), b, k);
            }
        }
    }
    for (int i = 0; i < m; i++){
        for (int k = 0; k < nclust; k++){
            f(i, k) = accu(log(fs.slice(k).row(i)));
        }
    } 
    for(int i = 0; i < inner_iters; i++){
      vec tsum=zeros<vec>(m);
      for (int ct=0; ct < ncat; ct++)
        for (int k=0; k < nclust; k++)
        {
          w_sum(k,ct)=0;
          for(int j = 0; j < m; j++)
            if (cat(j)==ct) w_sum(k,ct)=w_sum(k,ct)+p(k,j);
        }
            
        for (int k=0; k < nclust; k++){
          for (int ct=0; ct < ncat; ct++){
            a(k,ct) = alpha(k) + w_sum(k,ct);
          }
        }
        for(int j = 0; j < m; j++){
          for (int k=0; k < nclust; k++){
            p(k,j) = digamma(a(k,cat(j))) + f(j,k);
          }
          double maxv = max(p.col(j));
          for (int k=0; k < nclust; k++){
            p(k,j) = p(k,j) - maxv;
          }
        }
        for(int j = 0; j < m; j++){
          double s = accu(exp(p.col(j)));
          for (int k=0; k<nclust; k++){
            p(k,j) = exp(p(k,j)) / s;
          }
        }
    }  
  return List::create(Named("p") = p, Named("f") = f, Named("fs")=fs, Named("a")=a);
}



// [[Rcpp::export]]
RcppExport SEXP ebmme_cpp_binned(SEXP dat_, SEXP cat_, SEXP block_id_,
                                 SEXP H1_inv_, SEXP p_, SEXP alpha_,
                                 SEXP iters_, SEXP inner_iters_,SEXP nclust_,
                                 SEXP bins_, SEXP bin_data_, SEXP kde_nbins_)//, SEXP playvec_)
{
  
  mat dat = as<mat>(dat_);
  vec cat = as<vec>(cat_);
  mat bins = as<mat>(bins_);
  mat bin_data = as<mat>(bin_data_);
  int m = dat.n_rows;
  
  vec cat_id = unique(cat);
  int ncat = cat_id.n_elem;
  
  int iters = as<int>(iters_);
  int inner_iters = as<int>(inner_iters_);
  int nclust = as<int>(nclust_);
  int nbins = as<int>(kde_nbins_);
  
  mat H1_inv = as<mat>(H1_inv_);
  mat p = as<mat>(p_);
  
  vec alpha = as<vec>(alpha_);
  mat a = zeros<mat>(nclust,ncat);
  cube all_as = zeros<cube>(nclust,ncat,iters);

  vec alpha0=zeros<vec>(nclust);
  alpha0 = alpha;
  
  mat w_sum = zeros<mat>(nclust,ncat);
  
  vec block_id  = as<vec>(block_id_);
  vec ublock_id = unique(block_id);
  int B = ublock_id.n_elem;
  
  mat f = zeros<mat>(m, nclust);
  
  vec df = zeros<vec>(nclust);
  mat d2f = zeros<mat>(nclust, nclust);
  mat d2f_inv = zeros<mat>(nclust, nclust);
  
  vec lb = zeros<vec>(iters);
  
  cube fs_binned_save;
    
    for (int iter = 0; iter < iters; iter++){
        Rcpp::checkUserInterrupt();
        Rcout << '.';
        cube fs = zeros<cube>(m, B, nclust);
        cube fs_binned = zeros<cube>(nbins, B, nclust);
        vec p_sum = sum(p, 1);
        for (int b = 0; b < B; b++){
            mat p_binned = zeros<mat>(nclust, nbins);
            for (int i = 0; i < m; i++){
                for (int k = 0; k < nclust; k++){
                    p_binned(k, bins(i, b)) += p(k, i);
                }
            }
            uvec index  = find(block_id == ublock_id(b));
            mat H1_inv_b = H1_inv.submat(index, index);
            for(int i = 0; i < nbins; i++){
                rowvec zi = bin_data.row(i);
                for(int j = 0; j < nbins; j++){
                    rowvec zj = bin_data.row(j);
                    vec x = zi(index) - zj(index);
                    if (!x.has_nan()){
                        for (int k = 0; k < nclust; k++){
                            fs_binned(i, b, k) += p_binned(k, j) * exp(-.5 * x(0) * H1_inv_b(0,0) * x(0));
                        }
                    }
                }
            }
        }
        
        for (int k = 0; k < nclust; k++){
            for (int b=0; b < B; b++){
                uvec index  = find(block_id == ublock_id(b));
                mat H1_inv_b = H1_inv.submat(index, index);
                double det_H1_inv_b = arma::det(H1_inv_b);
                for (int i = 0; i < nbins; i++){
                    fs_binned(i, b, k) *= sqrt(det_H1_inv_b) / p_sum(k);
                }
            }
        }
        fs_binned_save = fs_binned;
        for (int i = 0; i < m; i++){
            for (int k = 0; k < nclust; k++){
                for (int b=0; b < B; b++){
                    fs(i, b, k) = fs_binned(bins(i, b), b, k);
                }
            }
        }
        
        for (int i = 0; i < m; i++){
            for (int k = 0; k < nclust; k++){
                f(i, k) = prod(fs.slice(k).row(i));
            }
        }
        
        // the inner loop to estimate the variational parameters
        for(int i = 0; i < inner_iters; i++){
            Rcpp::checkUserInterrupt();
            vec tsum=zeros<vec>(m);
            for (int ct=0; ct < ncat; ct++)
                for (int k=0; k < nclust; k++)
                {
                    w_sum(k,ct)=0;
                    for(int j = 0; j < m; j++)
                        if (cat(j)==ct) w_sum(k,ct)=w_sum(k,ct)+p(k,j);
                }
            
            for (int k=0; k < nclust; k++){
                for (int ct=0; ct < ncat; ct++){
                    a(k,ct) = alpha(k) + w_sum(k,ct);
                    all_as(k,ct,iter) = a(k,ct);
                }
            }
            for (int k=0; k < nclust; k++)
                for(int j = 0; j < m; j++)
                {
                    p(k,j) = exp(digamma(a(k,cat(j)))) * f(j,k);
                    tsum(j)=tsum(j)+p(k,j);
                }
            
            
            for (int k=0; k<nclust; k++)
                for(int j = 0; j < m; j++)
                    p(k,j) = p(k,j)/tsum(j);
        }
        
        // Newton-Raphson optimization step
        /* alpha=alpha0;
         for(int j = 0; j < inner_iters; j++)
         {
         df = d_object_f(alpha, a);
         d2f = d2_object_f(alpha);
         d2f_inv = inv(d2f);
         alpha = alpha - d2f_inv * df;
         }*/
        
        
        // Calculate lower bound:
        //lb(iter) = lgamma(accu(alpha)) - accu( lgamma_vec(alpha) ) +
        //accu( (alpha - 1) % (digamma_vec(a) - digamma(accu(a))) ) +
        // p_sum * (digamma(a(1)) - digamma(accu(a))) +
        //q_sum * (digamma(a(0)) - digamma(accu(a))) +
        //accu( (1 - p) % log(f.col(1)) ) +
        //accu( p % log(f.col(0)) ) +
        //( lgamma(a(0)) + lgamma(a(1)) - lgamma(accu(a)) ) -
        //( a(1) - 1 ) * ( digamma(a(1)) - digamma(accu(a)) ) -
        //( a(0) - 1 ) * ( digamma(a(0)) - digamma(accu(a)) ) -
        //accu( p % log(p) ) - accu( (1 - p) % log(1 - p) );
        
    }
    //fs_binned_save.reshape(5*4, 3, 1);
    //mat fs_binned = fs_binned_save.slice(0);
    
    Rcout << '\n';
    return List::create(Named("p") = p, Named("f") = f, Named("alpha") = alpha, Named("lb")=lb, Named("df")=df, Named("d2f")=d2f, Named("d2f_inv")=d2f_inv, Named("a")=a, Named("fs_binned")=fs_binned_save, Named("all_as")=all_as);
    
    
}

// [[Rcpp::export]]
RcppExport SEXP passcube(SEXP cube_)
{
    cube d = as<cube>(cube_);
    return List::create(Named("cube") = d);
}


// [[Rcpp::export]]
RcppExport SEXP predictlogsum(SEXP cat_, SEXP block_id_,
                        SEXP nclust_,
                        SEXP bins_,
                        SEXP fs_binned_, SEXP a_)
{
    mat bins = as<mat>(bins_);
    int m = bins.n_rows;
    
    vec cat = as<vec>(cat_);
    
    int nclust = as<int>(nclust_);
    
    mat p = zeros<mat>(nclust,m);
    vec block_id  = as<vec>(block_id_);
    vec ublock_id = unique(block_id);
    int B = ublock_id.n_elem;
    
    cube fs = zeros<cube>(m, B, nclust);
    
    mat f = zeros<mat>(m, nclust);
    cube fs_binned = as<cube>(fs_binned_);
    mat a = as<mat>(a_);
    for (int i = 0; i < m; i++){
        for (int k = 0; k < nclust; k++){
            for (int b=0; b < B; b++){
                fs(i, b, k) = fs_binned(bins(i, b), b, k);
            }
        }
    }
    for (int i = 0; i < m; i++){
        for (int k = 0; k < nclust; k++){
            f(i, k) = accu(log(fs.slice(k).row(i)));
        }
    }
    for(int j = 0; j < m; j++){
        for (int k=0; k < nclust; k++){
            p(k,j) = digamma(a(k,cat(j))) + f(j,k);
        }
        double maxv = max(p.col(j));
        for (int k=0; k < nclust; k++){
            p(k,j) = p(k,j) - maxv;
        }
    }
    for(int j = 0; j < m; j++){ 
      double s = accu(exp(p.col(j)));
      for (int k=0; k<nclust; k++){
        p(k,j) = exp(p(k,j)) / s;            
      }
    }
  return List::create(Named("p") = p, Named("f") = f, Named("fs")=fs, Named("a")=a);
}

