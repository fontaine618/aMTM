// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(BH)]]

#include <RcppArmadillo.h>
#include <boost/math/special_functions/erf.hpp>
#include <boost/math/distributions/normal.hpp>

using namespace Rcpp;

inline double evalTarget(Function target, arma::rowvec x, List parms) {
   return as<double>( wrap( target(x, parms) ) );
}

// Cholesky update
// Given the lower triangular matrix L obtained from the Cholesky decomposition of A=LL',
// updates L such that it corresponds to the decomposition of A + u*u'.
//

inline arma::mat chol_update(arma::mat& L, arma::vec& u) {
   unsigned int n = u.n_elem - 1;
   for (arma::uword i = 0; i < n; i++) {
      double r = sqrt(L(i,i) * L(i,i) + u(i) * u(i));
      double c = r / L(i, i);
      double s = u(i) / L(i, i);
      L(i, i) = r;
      L(arma::span(i + 1, n), i) =
         (L(arma::span(i + 1, n), i) + s * u.rows(i + 1, n)) / c;
      u.rows(i + 1, n) = c * u.rows(i + 1, n) -
         s * L(arma::span(i + 1, n), i);
   }
   L(n, n) = sqrt(L(n, n) * L(n, n) + u(n) * u(n));
   return L;
}


// Cholesky downdate
// Given the lower triangular matrix L obtained from the Cholesky decomposition of A=LL',
// updates L such that it corresponds to the decomposition of A - u*u'.
//
// NOTE: The function does not check that the downdating produces a positive definite matrix!

inline arma::mat chol_downdate(arma::mat& L, arma::vec& u) {
   unsigned int n = u.n_elem - 1;
   for (arma::uword i = 0; i < n; i++) {
      double r = sqrt(L(i,i) * L(i,i) - u(i) * u(i));
      double c = r / L(i, i);
      double s = u(i) / L(i, i);
      L(i, i) = r;
      L(arma::span(i + 1, n), i) =
         (L(arma::span(i + 1, n), i) - s * u.rows(i + 1, n)) / c;
      u.rows(i + 1, n) = c * u.rows(i + 1, n) -
         s * L(arma::span(i + 1, n), i);
   }
   L(n, n) = sqrt(L(n, n) * L(n, n) - u(n) * u(n));return L;
}
// [[Rcpp::export]]
List aMTMsample(Function target,
                int N,
                int d,
                int K,
                arma::vec x0,
                arma::cube sig0,
                arma::mat mu0,
                arma::vec lam0,
                int adapt,
                int proposal,
                double accrate,
                double gamma,
                List parms,
                double beta) {
   //--------------------------------------------------------------------
   // DECLARE VARIABLES AND INITIALIZE
   arma::vec w(K);            //vector of weights
      w.ones();
   double sw=0;               //sum of weights
   arma::vec wt(K);           //vector of reverse weights
      wt.ones();
   double swt=0;              //sum of reverse weights
   arma::vec wb(K);           //vector of standardized weights, i.e. selection probabilities
      wb.ones();
   arma::mat Sy(N,K);         //contains the selected candidate (1 if selected, 0 otherwise)
      Sy.zeros();
   arma::vec sel(N);          //vector containing which candidate was selected (k_n)
      sel.zeros();
   arma::vec acc(N);          //vector of acceptatnce probability (alpha_n^(k_n))
      acc.zeros();
   arma::vec alpha(K+1);      //running acceptance probability for each candidate and globally (K+1)
      alpha.zeros();
   arma::mat X(N,d);          //matrix of the sampled chain
      X.zeros();
      X.row(0) = x0.t();      //initialize to initial value
   arma::mat mu=mu0;          //mean parameters, initialize to initial value
   arma::vec lam=lam0;        //scale parameter, initialize to initial value
   arma::cube S(d,d,K);       //square root of variances, itilialize to initial value
   for(int k=0;k<K;k++){
      S.slice(k) = arma::chol(sig0.slice(k)).t();
   }
   arma::mat U(d,K);          //standard normal vectors for candidates
      U.zeros();
   arma::mat Ut(d,K);         //standard normal vectors for reference points
      Ut.zeros();
   arma::mat Y(d,K);          //candidates
      Y.zeros();
   arma::mat Yt(d,K);         //reference points
      Yt.zeros();
   arma::vec u(d);            //uniform vector for PIT
      u.zeros();
   arma::rowvec y;            //to pass to R when evaluating target
      y.zeros();
   double gam=0;              //adaptation step
   double un=0.0;             //uniform RV for selection
   double a=0.0;              //temp acceptance probability
   int s=0;                   //selection index
   //--------------------------------------------------------------------
   // SPECIFIC INITIALIZATIONS AND PRECOMPUTATIONS FOR ADAPTATION
   //if we do not use a scale parameter
   if(adapt==2){
      lam.ones();
   }
   //for QMC
   boost::math::normal norm;  //normal distribution object
   int ai = floor(K/3);       //Koborov integer parameter
   arma::vec Ua(d);           //base vector for Koborov rule (others are multiple of this one and mod1)
   if(proposal == 2){
      for(int i=0;i<d;i++)Ua(i) = pow(ai,i) /K;
   }
   //for extremely antithetic

   //--------------------------------------------------------------------
   // MCMC ITERATION
   for(int n=1;n<N;n++){
      //--------------------------------------------------------------------
      // MTM SAMPLING STEP
      // Candidates sampling normal standard by correlation structure
         switch (proposal){
         case 0:
            //independent
            for(int k=0;k<K;k++){u = arma::randn(d);U.col(k) = u;}
            break;
         case 1:
            //common RV
            u = arma::randn(d);U = arma::repmat(u,1,K);
            break;
         case 2:
            //QMC
            u = arma::randu(d);
            for(int k=0;k<K;k++){
               U.col(k) = (u + k*Ua) - floor(u + k*Ua);
               for(arma::uword i=0; i < u.n_elem;i++) U(i,k) = quantile(norm,U(i,k));
            }
            break;
         }
      // compute the candidates and weights
         for(int k=0;k<K;k++){
            Y.col(k) = X.row(n-1).t() + sqrt(lam(k)) * S.slice(k) * U.col(k);
            y = Y.col(k).t();
            w(k) = evalTarget(target,y,parms);
            w(k) = w(k) + beta * sum(-0.5*(log(2 * M_PI) + U.col(k)%U.col(k)));
         }
         w = exp(w);
      // compute weights
         sw = sum(w);
         if(sw>0.0)wb = w/sw;
         else wb=arma::ones(K)/K;
      // proposal selection
         un = arma::randu<double>();
         sw=0;s=0;
         for(int k=0;k<K;k++){
            sw=sw+wb(k);
            if(sw>un && un>sw-wb(k))s=k;
         }
         Sy(n,s)=1.0;sel(n)=s;
      // reference points sampling normal standard by correlation structure
         switch (proposal){
         case 0:
            //independent
            for(int k=0;k<K;k++){u = arma::randn(d);Ut.col(k) = u;}
            break;
         case 1:
            //common RV
            //u = arma::randn(d);Ut = arma::repmat(u,1,K); (was wrong)
            Ut=-U;
            break;
         case 2:
            //QMC
            for(int k=0;k<K;k++){
               Ut.col(k) = (-U.col(s) + k*Ua) - floor(-U.col(s) + k*Ua);
               for(arma::uword i=0; i < u.n_elem;i++) Ut(i,k) = quantile(norm,Ut(i,k));
            }
            break;
         }
      // compute the reference points and weights
         for(int k=0;k<K;k++){
            if(k==s) Yt.col(k) = X.row(n-1).t();
            else Yt.col(k) = Y.col(s) + sqrt(lam(k)) * S.slice(k) * Ut.col(k);
            y = Yt.col(k).t();
            wt(k) = evalTarget(target,y,parms);
            wt(k) = wt(k) + beta * sum(-0.5*(log(2 * M_PI) + Ut.col(k)%Ut.col(k)));
         }
         wt = exp(wt);
         swt = sum(wt);
      // MTM acceptatance probability
         if(swt==0.0)a=1.0;
         else a=sum(w)/swt;
         if(a>1.0)a=1.0;
      // acceptation/rejection of the proposal
         un = arma::randu<double>();
         if(un<a){
            X.row(n)=Y.col(s).t();
            acc(n)=1;
         }else X.row(n)=X.row(n-1);

      //--------------------------------------------------------------------
      // ADAPTATION STEP
      // adaptation rate
         gam = std::pow((double)n,-gamma);
         if(gam>0.99)gam=0.99;
      // adaptation of kernel s, adapt = 0 => no adaptation

         switch (adapt){
            case 1:
               // alg 1
               u = (X.row(n).t() - mu.col(s));
               mu.col(s) = mu.col(s) + gam * u;
               u = u * sqrt(gam / (1-gam));
               chol_update(S.slice(s),u);
               S.slice(s)=S.slice(s)*sqrt(1-gam);
               lam(s) = exp(log(lam(s)) + gam * (a-accrate));
            break;
            case 2:
               // alg 2
               u = S.slice(s) * U.col(s);
               u = u* sqrt(gam * std::pow(fabs(a-accrate), 1.0/1.0)) /
                  arma::norm(U.col(s));
               if(a-accrate > 0.0) {
                  chol_update(S.slice(s),u);
               }
               else{
                  chol_downdate(S.slice(s),u);
               }
               if(any(S.slice(s).diag())<1e-7) S.slice(s).diag() += 1e-7;
            break;
            case 3:
               // alg 3
               u = (X.row(n).t() - X.row(n-1).t());
               u = u * sqrt(gam / (1-gam));
               chol_update(S.slice(s),u);
               S.slice(s)=S.slice(s)*sqrt(1-gam);
               lam(s) = exp(log(lam(s)) + gam * (a-accrate));
            break;
         }
         alpha(s) = alpha(s) + gam * (a-alpha(s));
         alpha(K) = alpha(K) + gam * (a-alpha(K));


   }
   // return the covariances
      arma::cube Sig(d,d,K);
      for(int k=0;k<K;k++){
         Sig.slice(k) = S.slice(k) * S.slice(k).t();
      }
   // OUTPUT LIST
      return List::create(X,sel,acc,alpha,mu,Sig,
                        lam,accrate,N,K,d,Sy);
}

