//
// Created by andreas on 01.05.15.
//

#include "LanczosSolver.h"

#include <Eigen/Sparse>
#include <Eigen/Eigenvalues>

using namespace std;
using namespace Eigen;

//template<typename MatrixMult>
//void LanczosSolver::compute(const MatrixMult& mult, int dim, int nev, bool evec) {
//    // SETUP
//    const int MAX_ITER(100);
//    vector<VectorXd> KrylovBasis;
//    KrylovBasis.reserve(MAX_ITER+1); // allocate enough memory for recovery of eigenvectors
//    KrylovBasis.push_back(VectorXd::Zero(dim)); // zero vector in iteration;
//    KrylovBasis.push_back(VectorXd::Random(dim)); // initial vector
//    KrylovBasis.back().normalize(); // needs to be normalized
//    SparseMatrix<double> matrix(MAX_ITER, MAX_ITER); // the Lanczos matrix to be solved
//    matrix.reserve(3); // the matrix is tridiagonal
//    // CALCULATE MATRIX
//    double tempA(0), tempB(0);
//    for (int i(0); i!=MAX_ITER-1; ++i) {
//        tempW = mult(KrylovBasis.back().adjoint());
//        tempA = tempW*KrylovBasis.back();
//        tempW = tempW-KrylovBasis.back()*tempA-KrylovBasis[i]*tempB;
//        tempB = tempW.norm();
//        tempW.normalize();
//        KrylovBasis.push_back(tempW);
//        matrix.insert(i,i) = tempA;
//        matrix.insert(i,i+1) = tempB;
//        matrix.insert(i+1,i) = tempB;
//    }
//    tempW=mult(KrylovBasis.back()).adjoint;
//    matrix.insert(MAX_ITER, MAX_ITER)=tempW*KrylovBasis.back();
//    // SOLVE
//    EigenSolver<MatrixXd> solver(matrix, evec);
//    eigenValues=solver.eigenvalues();
//    eigenVectors=solver.eigenvectors();
//}
