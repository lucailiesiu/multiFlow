#include <fstream>
#include "tools.h"
#include "lambda.h"

#ifndef OBSPARAM_H
#define OBSPARAM_H

struct ObsParam
{
    enum solutionType solution;
    VectorXd eigenV;
    MatrixXd Hessian, U;
    double H, tau, k; 
    double n;
    double nT;
    double r;
    double alpha;
    double fNL;

    ObsParam();
    ObsParam(Lambda &l, enum solutionType sol);
    void Hess(Lambda l);
    void spectralGenerate(Lambda l);
    void spectralIndex(Lambda l);
    void spectralTensorIndex(Lambda l);
    void spectralTilt(Lambda l);
    void spectralAlpha(Lambda l);
    void spectralfNL(Lambda l);
    void printObs(std::ofstream &file);
};
#endif

