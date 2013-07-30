#include </u/liliesiu/MFC_code/Eigen/Dense>
#include </u/liliesiu/MFC_code/Eigen/Eigenvalues>
#include <fstream>

#ifndef LAMBDA_H
#define LAMBDA_H

using namespace Eigen;

class Lambda
{
public:
    int M;
    int kTrunc;
    double **lambda;
    double *fields; 

    Lambda();
    Lambda(const Lambda& l);
    Lambda(int giveM, int giveK, double **giveLambda, double *giveFields);
    Lambda(int giveM, int giveK, double **lambdaLimitsUp, double **lambdaLimitsDown);
    Lambda(int giveM, int giveK);
    ~Lambda();
    Lambda & operator= (const Lambda& l);
    
    double epsilon();
    MatrixXd lambdaMatrix();
    Lambda lSum(Lambda l);
    Lambda lProd(double scalar);
    void printLambda(std::ofstream &file);
    void printEpsilon(std::ofstream &file);
    void printFields(std::ofstream &file); 

private:
    void toolsIncrease(int *a, int position);
    int toolsIndex(const int* a, const int number);
    void generatePerm(const int* a, const int number, const  int index);
    
};

#endif
    
