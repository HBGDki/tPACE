#include <RcppEigen.h>
#include <map>          // to map kernels to integers for the switch
#include <string>       // to read in the kernel name
#include <vector>       // to use vectors
#include <algorithm>    // to get the intersect, sort, lower_bound, upper_bound
#include <ANN/ANN.h>

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]

Eigen::MatrixXd RmullwlskCCtree( const std::string kernel_type, const Eigen::Map<Eigen::MatrixXd> & tPairs, const Eigen::Map<Eigen::MatrixXd> & cxxn, const Eigen::Map<Eigen::VectorXd> & win,  const Eigen::Map<Eigen::VectorXd> & xgrid, const Eigen::Map<Eigen::VectorXd> & ygrid ) { 
  // Assumes the first row of tPairs is sorted in increasing order.
  // tPairs : xin (in MATLAB code)
  // cxxn : yin (in MATLAB code)
  // xgrid: out1 (in MATLAB code)
  // ygrid: out2 (in MATLAB code)

  const double eps = 0.0; // for approximate NN search
  const double invSqrt2pi=  1./(sqrt(2.*M_PI));

  // Map the kernel name so we can use switches  
  std::map<std::string,int> possibleKernels; 
  possibleKernels["epan"]    = 1;   possibleKernels["rect"]    = 2;
  possibleKernels["gauss"]   = 3;   possibleKernels["gausvar"] = 4; 
  possibleKernels["quar"]    = 5; 
   
  // The following test is here for completeness, we mightwant to move it up a 
  // level (in the wrapper) in the future. 

  // If the kernel_type key exists set KernelName appropriately
  int KernelName = 0;
  if ( possibleKernels.count( kernel_type ) != 0){ 
    KernelName = possibleKernels.find( kernel_type )->second; //Set kernel choice
  } else {
  // otherwise use "epan"as the kernel_type 
    //Rcpp::Rcout << "Kernel_type argument was not set correctly; Epanechnikov kernel used." << std::endl;
    Rcpp::warning("Kernel_type argument was not set correctly; Epanechnikov kernel used.");
    KernelName = possibleKernels.find( "epan" )->second;
  }

  // Start the actual smoother here  
  const unsigned int xgridN = xgrid.size();  
  const unsigned int ygridN = ygrid.size();  
  const unsigned int n = tPairs.cols();
  
  double* tDat = new double[n * 2];
  memcpy(tDat, tPairs.data(), n * 2 * sizeof(double));
  // Eigen::Map<Eigen::MatrixXd>(tDat, tPairs.rows(), tPairs.cols()) = tPairs;
  
  Eigen::MatrixXd mu(xgrid.size(), ygrid.size());
  mu.setZero();    

  // if (KernelName != 3) {
    // build tree
    ANNpointArray	dataTree = new ANNpoint[n];
    for (unsigned int i = 0; i < n; i++) {
      dataTree[i] = tDat + i * 2; // 2-D input xin
    }
    ANNkd_tree	*kd_tree = new ANNkd_tree(dataTree, n, 2);	

    // get a buffer to store the nearest neighbors
    ANNidxArray	nnIndex = new ANNidx[n];		// Near neighbor indices
    unsigned int NN = 0; // number of nearest neighbors.

  // }

  for (unsigned int j = 0; j != ygridN; ++j) { 
    for (unsigned int i = 0; i != xgridN; ++i) {  
      double outPt[2] = {xgrid[i], ygrid[j]};
      
      //locating local window (LOL) (bad joke)
      std::vector <unsigned int> indx; 
      //if the kernel is not Gaussian
      if ( KernelName != 3) { 
        // search tree
        NN = kd_tree->annkFRSearch(outPt, 1, n, nnIndex, NULL); //, eps);
        for (int k = 0; k < NN; k++) 
          indx.push_back(nnIndex[k]);
        // std::sort(indx.begin(), indx.end()); // may be faster

      } else{ // just get the whole deal
        for (unsigned int y = 0; y != tPairs.cols(); ++y){
          indx.push_back(y);
        }
      }  

      unsigned int indxSize = indx.size();
      Eigen::VectorXd lw(indxSize);  
      Eigen::VectorXd ly(indxSize);
      Eigen::MatrixXd lx(2,indxSize);

      for (unsigned int u = 0; u !=indxSize; ++u){ 
        lx.col(u) = tPairs.col(indx[u]); 
        lw(u) = win(indx[u]); 
        ly(u) = cxxn(indx[u]); 
      }

      // check enough points are in the local window 
      unsigned int meter=1;  
      for (unsigned int u =0; u< indxSize; ++u) { 
        for (unsigned int t = u + 1; t < indxSize; ++t) {
          if ( (lx(0,u) !=  lx(0,t) ) || (lx(1,u) != lx(1,t) ) ) {
            meter++;
          }
        }
        if (meter >= 3) { 
          break; 
        }
      }
   
      //computing weight matrix 
      if (meter >=  3) { 
        Eigen::VectorXd temp(indxSize);
        Eigen::MatrixXd llx(2, indxSize );  
        llx.row(0) = (lx.row(0).array() - xgrid(i));
        llx.row(1) = (lx.row(1).array() - ygrid(j));

        //define the kernel used 

        switch (KernelName){
          case 1: // Epan
            temp=  ((1-llx.row(0).array().pow(2))*(1- llx.row(1).array().pow(2))).array() * 
                   ((9./16)*lw).transpose().array(); 
            break;  
          case 2 : // Rect
            temp=(lw.array())*.25 ; 
            break;
          case 3 : // Gauss
            temp = ((-.5*(llx.row(1).array().pow(2))).exp()) * invSqrt2pi  *   
                   ((-.5*(llx.row(0).array().pow(2))).exp()) * invSqrt2pi  *
                   (lw.transpose().array()); 
            break;
          case 4 : // GausVar
            temp = (lw.transpose().array()) * 
		   ((-0.5 * llx.row(0).array().pow(2)).array().exp() * invSqrt2pi).array() *
                   ((-0.5 * llx.row(1).array().pow(2)).array().exp() * invSqrt2pi).array() * 
                   (1.25 - (0.25 * (llx.row(0).array().pow(2))).array())  * 
                   (1.50 - (0.50 * (llx.row(1).array().pow(2))).array()); 
            break;
          case 5 :  // Quar
              temp = (lw.transpose().array()) * 
                     ((1.-llx.row(0).array().pow(2)).array().pow(2)).array() *
                     ((1.-llx.row(1).array().pow(2)).array().pow(2)).array() * (225./256.);
            break;
        } 

        // make the design matrix
        Eigen::MatrixXd X(indxSize ,3);
        X.setOnes();    
        X.col(1) = llx.row(0);
        X.col(2) = llx.row(1);
        Eigen::LDLT<Eigen::MatrixXd> ldlt_XTWX(X.transpose() * temp.asDiagonal() *X);

        Eigen::VectorXd beta = ldlt_XTWX.solve(X.transpose() * temp.asDiagonal() * ly);  
        mu(i,j)=beta(0); 
      } else if(meter < 3){
          Rcpp::warning("No enough points in local window, please increase bandwidth.");
          mu(i, j) = NAN;
      }
    }
  }

  delete[] tDat;
  // if (KernelName != 3) {
    delete kd_tree;
    delete[] dataTree;
    delete[] nnIndex;
    annClose();
  // }
  return ( mu ); 
}

 
