# include <RcppEigen.h>
# include <Eigen/Geometry>
# include <Eigen/Eigenvalues>
# include <cmath>
// [[Rcpp::depends("RcppEigen")]]
using namespace std;
using namespace Rcpp;
using namespace Eigen;

double Eigennorm (MatrixXd A) {
  MatrixXd tAA = A.adjoint() * A;
  double sv = ::sqrt(tAA.eigenvalues().real()(0));
  return sv;
}

Eigen::MatrixXd MatFind (Eigen::MatrixXd A, double ZeroThres) {
  Eigen::MatrixXd B;
  B.setZero(A.rows(), A.cols());
  Eigen::MatrixXd Atrunc = (A.array() < ZeroThres).select(B, A);
  return(Atrunc);
}

Eigen::MatrixXd MatFindlb (Eigen::MatrixXd A, double lb) {
  Eigen::MatrixXd B;
  B.setOnes(A.rows(), A.cols());
  B *= lb;
  Eigen::MatrixXd Alb = (A.array() < lb).select(B, A);
  return(Alb);
}

Rcpp::List SMatC (Eigen::MatrixXd X, Eigen::VectorXi clu) {
  int curidx = 0;
  int K = clu.size();
  Eigen::MatrixXd subX, Sw, Sb;
  Sw.setZero(X.rows(), X.rows());
  Sb.setZero(X.rows(), X.rows());
  Eigen::VectorXd AllRowMean, subRowMean;
  
  AllRowMean = X.rowwise().sum()/X.cols();
  for (int k=0; k<K; k++) {
    subX = X.middleCols(curidx, clu[k]);
    subRowMean = subX.rowwise().sum()/clu[k];
    subX = subX.array().colwise() - subRowMean.array();
    Sw += subX * subX.adjoint();
    subRowMean -= AllRowMean;
    Sb += clu[k] * (subRowMean * subRowMean.adjoint());
    curidx += clu[k];
  }
  
  Rcpp::List Smats;
  Smats["Sw"] = Sw;
  Smats["Sb"] = Sb;
  return (Smats);
}



// [[Rcpp::export()]]
Eigen::MatrixXd PNMF_EucDistC(Eigen::MatrixXd X, Eigen::MatrixXd W_init, double tol, int maxIter, bool verboseN, double zerotol) {
  // initialization
  Eigen::MatrixXd W = W_init;
  Eigen::MatrixXd W_old = W_init;
  
  Eigen::MatrixXd XX = X * X.transpose();
  Eigen::MatrixXd XXW, SclFactor;
  
  double diffW;
  int iter = 0;
  
  // iterations
  for (iter = 0; iter < maxIter; iter++) {
    W_old = W;
    
    XXW = XX * W;
    SclFactor = W * (W.transpose() * XXW) + (XXW * (W.transpose() * W));
    // QuotientLB
    SclFactor = MatFindlb(SclFactor, zerotol);
    // SclFactor = W * (W.transpose() * (XX * W)) + XX * W * CrossProdC(W, nG, nD);
    SclFactor = XXW.cwiseQuotient(SclFactor);
    W = W.cwiseProduct(SclFactor);
    
    W /= Eigennorm(W);
    W = MatFind(W, zerotol);
    
    diffW = (W_old - W).norm() / W_old.norm();
    if (diffW < tol) {
      break;
    }
  }
  
  if (verboseN) {
    Rprintf("%d iterations used.\n", (iter+1));
  }
  return (W);
  
}

// [[Rcpp::export()]]
Eigen::MatrixXd PNMF_KLC(Eigen::MatrixXd X, Eigen::MatrixXd W_init, double tol, int maxIter, bool verboseN, double zerotol) {
  // initialization
  Eigen::MatrixXd W = W_init;
  Eigen::VectorXd Xsum = X.rowwise().sum();
  Eigen::VectorXd XW, sW;
  Eigen::MatrixXd WWX, Z, SclFactor, W_old, bsxfuntimes, bsxfunplus;
  
  double diffW;
  int iter = 0;
  
  // iterations
  for (iter = 0; iter < maxIter; iter++) {
    W_old = W;
    // QuotientLB
    WWX = MatFindlb(W * (W.transpose() * X), zerotol);
    Z = X.cwiseQuotient(WWX);
    
    sW = W.colwise().sum();
    bsxfuntimes = Xsum * sW.adjoint();
    XW = W.adjoint() * Xsum;
    bsxfunplus = bsxfuntimes.array().rowwise() + XW.transpose().array();
    // QuotientLB
    bsxfunplus = MatFindlb(bsxfunplus, zerotol);
    SclFactor = Z * (X.adjoint() * W) + X * (Z.adjoint() * W);
    SclFactor = SclFactor.cwiseQuotient(bsxfunplus);
    W = W.cwiseProduct(SclFactor);
    
    W /= Eigennorm(W);
    W = MatFind(W, zerotol);
    
    diffW = (W_old - W).norm() / W_old.norm();
    if (diffW < tol) {
      break;
    }
  }
  
  if (verboseN) {
    Rprintf("%d iterations used.\n", (iter+1));
  }
  return (W);
  
}



// [[Rcpp::export()]]
Eigen::MatrixXd DPNMFC(Eigen::MatrixXd X, Eigen::MatrixXd W_init, double tol, int maxIter, bool verboseN, double zerotol, Eigen::MatrixXd Xord, Eigen::VectorXi clunum, double mu, double lambda) {
  // initialization
  Eigen::MatrixXd W = W_init;
  //VectorXd ;
  Eigen::MatrixXd XX = X * X.transpose();
  Eigen::MatrixXd W_old, XXW, SclFactor, Sb, Sw, Sb_pos, Sw_pos;
  
  double diffW;
  int iter = 0;
  
  // calculation of Sb, Sw
  Rcpp::List Smats = SMatC(Xord, clunum);
  Sw = Smats["Sw"];
  Sb = Smats["Sb"];
  
  Sb_pos = lambda * Sw - Sb;
  Sw_pos = MatFind(Sb_pos, 0.0);
  Sb_pos = Sw_pos - Sb_pos;
  
  // iterations
  for (iter = 0; iter < maxIter; iter++) {
    W_old = W;
    
    XXW = XX * W;
    SclFactor = W * (W.transpose() * XXW) + (XXW * (W.transpose() * W)) + mu*(Sw_pos * W);
    // QuotientLB
    SclFactor = MatFindlb(SclFactor, zerotol);
    SclFactor = (2.0 * XXW + mu*(Sb_pos * W)).cwiseQuotient(SclFactor);
    W = W.cwiseProduct(SclFactor);
    
    W /= Eigennorm(W);
    W = MatFind(W, zerotol);
    
    diffW = (W_old - W).norm() / W_old.norm();
    if (diffW < tol) {
      break;
    }
  }
  
  if (verboseN) {
    Rprintf("%d iterations used.\n", (iter+1));
  }
  return (W);
  
}










