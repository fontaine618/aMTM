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
List aMTMsample(Function target,             // target density
                int N,                       // sample size
                int d,                       // dimension
                int K,                       // number of proposals
                arma::vec x0,                // initial state
                arma::cube sig0,             // initial covariances
                arma::mat mu0,               // initial means
                arma::vec lam0,              // initial scale parameter
                int adapt,                   // update function (0= nodapt, 1=AM, 2=ASWAM, 3=RAM)
                int global,                 // adapt global component 
                int scale,                  // adapt non selected
                int local,                  // use local steps (AM and ASWAM only)
                int proposal,                // type of proposal (0=independant,1=common RV,2=QMC,3=EA)
                double accrate,              // target acceptance rate (ASWAM and RAM)
                double gamma,                // adaptation step decrease power
                List parms,                  // list of parameters to pass to target
                double beta) {               // for the weight function
   //--------------------------------------------------------------------
   // DECLARE VARIABLES AND INITIALIZE
   arma::vec w(K);            //vector of weights
      w.ones();
   double sw=0.0;               //sum of weights
   arma::vec wt(K);           //vector of reverse weights
      wt.ones();
   double swt=0.0;              //sum of reverse weights
   arma::vec wb(K);           //vector of standardized weights, i.e. selection probabilities
      wb.ones();
   arma::vec Sy(K);           //vector of running selection proportions
      Sy.fill(1.0/K);
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
   arma::vec lam;             //scale parameter
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
   double gam=0.0;            //adaptation step
   double un=0.0;             //uniform RV for selection
   double a=0.0;              //temp acceptance probability
   int s=0;                   //selection index
   //--------------------------------------------------------------------
   // SPECIFIC INITIALIZATIONS AND PRECOMPUTATIONS FOR ADAPTATION
   // initialize scale parameter
   switch(adapt){
   case 1: //AM
      lam.fill((2.38)*(2.38)/d);
   case 2: //ASWAM
      lam.ones();
   case 3: //RAM
      lam = lam0;
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
         sw=0.0;s=0;
         for(int k=0;k<K;k++){
            sw=sw+wb(k);
            if(sw>un && un>sw-wb(k))s=k;
         }
         //Sy(n,s)=1.0;
         sel(n)=s;
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
      // adaptation cycling through the proposals
         for(int k=0;k<K;k++){
         // adapt covariance if selected of if first and global adaptation is enables
            if(s == k || (k==1 && global==1)){
            // AM and ASWAM are similar so we group them
               if(adapt == 1 || adapt == 2){
               // the update steps depends on the local trigger
                  switch(local){
                     case 1:
                        u = (X.row(n).t() - X.row(n-1).t());
                        break;
                     case 0:
                        u = (X.row(n).t() - mu.col(k));
                        mu.col(k) = mu.col(k) + gam * u;
                        break;
                  }
               // update the covariance
                  u = u * sqrt(gam / (1-gam));
                  chol_update(S.slice(k),u);
                  S.slice(k)=S.slice(k)*sqrt(1-gam);
               // update the scaling parameter if ASWAM
                  if(adapt == 2){
                     lam(k) = exp(log(lam(k)) + gam * (a-accrate));
                  }
               }
            // RAM update
               if(adapt == 3){
                  u = S.slice(k) * U.col(k);
                  u = u* sqrt(gam * std::pow(fabs(a-accrate), 1.0/1.0)) /
                     arma::norm(U.col(k));
                  if(a-accrate > 0.0) {
                     chol_update(S.slice(k),u);
                  }
                  else{
                     chol_downdate(S.slice(k),u);
                  }
                  if(any(S.slice(k).diag())<1e-7) S.slice(k).diag() += 1e-7;
               }
            // adapt = 0 we do nothing
            }
         // adapt scaling of not selected if adaptation is enabled
            // Sy[k] = Sy[k] + gam * (double(s==k) - 1.0/K);
            if(n>100){
               Sy(k) = 0.0;
               for(int i = n-99;i<=n;i++){
                  Sy(k) += (sel(i) == k)/100.0;
               }
            }
            //if(Sy[k]<0.0)Sy[k]=0.0;
            //if(Sy[k]>1.0)Sy[k]=1.0;
            if(s!=k && scale==1 && Sy(k) < 0.1/K){
               lam(k) = exp(log(lam(k)) + gam );
            }
         }
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

