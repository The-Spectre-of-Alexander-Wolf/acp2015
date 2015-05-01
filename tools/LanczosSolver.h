//
// Created by andreas on 01.05.15.
//

#ifndef ACP_LANCZOSSOLVER_H
#define ACP_LANCZOSSOLVER_H

#include <vector>
#include <complex>
#include <Eigen/Dense>

class LanczosSolver {
public:
    LanczosSolver() {}
    // Matrix must provide a matrix-vector multiplication
    template<typename Matrix>
    void compute(Matrix mult, int dim, int nev, bool evec=false) {
        // SETUP
        const int MAX_ITER(100);
        std::vector<Eigen::VectorXd> KrylovBasis;
        KrylovBasis.reserve(MAX_ITER+1); // allocate enough memory for recovery of eigenvectors
        KrylovBasis.push_back(Eigen::VectorXd::Zero(dim)); // zero vector in iteration;
        KrylovBasis.push_back(Eigen::VectorXd::Random(dim)); // initial vector
        KrylovBasis.back().normalize(); // needs to be normalized
        Eigen::MatrixXd matrix(MAX_ITER, MAX_ITER); // the Lanczos matrix to be solved
//        matrix.reserve(3); // the matrix is tridiagonal
        // CALCULATE MATRIX
        double tempA(0), tempB(0);
        for (int i(0); i!=MAX_ITER-1; ++i) {
            tempW = mult(KrylovBasis.back());
            tempA = abs(tempW.dot(KrylovBasis.back()));
            tempW = tempW-KrylovBasis.back()*tempA-KrylovBasis[i]*tempB;
            tempB = tempW.norm();
            tempW.normalize();
            KrylovBasis.push_back(tempW);
            matrix(i,i) = tempA;
            matrix(i,i+1) = tempB;
            matrix(i+1,i) = tempB;
        }
        tempW=mult(KrylovBasis.back());
        matrix(MAX_ITER, MAX_ITER)=abs(tempW.dot(KrylovBasis.back()));
        // SOLVE
        Eigen::EigenSolver<Eigen::MatrixXd> solver(matrix, evec);
        eigenValues=solver.eigenvalues();
//        eigenVectors=solver.eigenvectors();
    }
    std::vector<std::complex<double> > eigenvalues();
    Eigen::VectorXcd eigenValues;
    Eigen::MatrixXcd eigenVectors;
private:
    Eigen::VectorXd tempW;
};


#endif //ACP_LANCZOSSOLVER_H
