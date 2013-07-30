#include </u/liliesiu/MFC_code/Eigen/Dense>
#include </u/liliesiu/MFC_code/Eigen/Eigenvalues>
#include <vector>
#include <fstream>

#include "tools.h"
#include "lambda.h"
#include "ObsParameters.h"

#ifndef FLOW_H
#define FLOW_HA

using namespace Eigen;
using std::vector;

class Flow
{
public:
    double dN;
    int M;
    int kTrunc;
    int nrIter;
    Lambda initial;
    vector<Lambda> memorized;
    vector<double> memorizedEps; 
    ObsParam observable;

    Flow();
    Flow(int M, int kTrunc,
         double **lambdaLimitsUp, double **lambdaLimitsDown);
    Flow(int M, int kTrunc, Lambda l);


    void generate();
    void printFlow(std::ofstream &file); 
private:
    double mean;
    double stddev;
    vector<double> means;
    vector<double> stddevs; 
    
    void toolsAdd(int *a, int position);
    void f(Lambda l, Lambda &dl);
    enum solutionType checkMatrix();
    enum solutionType checkFixedLambda();
    enum solutionType checkSol();
    void rk45(Lambda &RKresult);
};

#endif
