#include </u/liliesiu/MFC_code/Eigen/Dense>
#include </u/liliesiu/MFC_code/Eigen/Eigenvalues>
#include <stdlib.h>
#include <cmath>
#include <fstream>

#include "ObsParameters.h"
#include "lambda.h"
#include "tools.h"


void ObsParam::spectralIndex(Lambda l)
{
    int M = l.M;
    double epsilon = l.epsilon();
    double sumDown = 0, sumUp = 0;
    for (int i = 0; i < M; i++)
        for (int j = 0; j < M; j++)
            for (int a = 0; a < M; a++)
            {
                double sumC = (1 - (epsilon + eigenV(a)) * CE + eigenV(a));
                double sum = (U(i,a)*U(j,a))/(l.lambda[0][i]*l.lambda[0][j]) * pow(sumC, 2);
                sumDown += sum;
                double sumK = 0;
                for (int k = 0; k < M; k++)
                {
                    sumK += l.lambda[0][k] * ((l.lambda[1][l.M * i + k]/l.lambda[0][i]) +
                                              (l.lambda[1][l.M * j + k]/l.lambda[0][j]));
                }
                sumUp += (U(i,a)*U(j,a))/(l.lambda[0][i]*l.lambda[0][j]) *
                    sumC * ((sumK - 2 * epsilon) * sumC + 2 * epsilon);
            }
    n = 1 - ((1/(1 - epsilon)) * (sumUp/sumDown) + 2 * epsilon);
}

void ObsParam::spectralTensorIndex(Lambda l)
{
    double epsilon = l.epsilon();
    double sumNT = 0;
    for (int i = 0; i < l.M; i++)
        for (int j = 0; j < l.M; j++)
        {
            sumNT += l.lambda[0][i] * l.lambda[0][j] * l.lambda[1][l.M * i + j];
        }
    nT = 1 - 2 * epsilon - (3 + C) * epsilon * epsilon - (1 + C) * sumNT;
}

void ObsParam::spectralTilt(Lambda l)
{
    int M = l.M;
    double epsilon = l.epsilon();
    double sumDown = 0;
    for (int i = 0; i < M; i++)
        for (int j = 0; j < M; j++)
            for (int a = 0; a < M; a++)
            {
                double sumC = (1 - (epsilon + eigenV(a)) * CE + eigenV(a));
                double sum = (U(i,a)*U(j,a))/(l.lambda[0][i]*l.lambda[0][j]) * pow(sumC, 2);
                sumDown += sum;
            }
    r = 4 * (1 - epsilon * CE)/sumDown ;
}

void ObsParam::spectralAlpha(Lambda l)
{
    double epsilon = l.epsilon();
    double dedN = - epsilon * epsilon;
    for (int i = 0; i < l.M; i++)
        for (int j = 0; j < l.M; j++)
        {
            dedN += l.lambda[0][i] * l.lambda[0][j] * l.lambda[1][l.M * i + j];
        }
    alpha = 2 * dedN * (2 + n/(1 - epsilon))/((1-epsilon) * (1 - epsilon)) ;
}

void ObsParam::spectralfNL(Lambda l)
{
    fNL = 0;
}

void ObsParam::Hess(Lambda l)
{
    int M = l.M;
    double epsilon = l.epsilon();
    Hessian = MatrixXd(l.M, l.M);
    for (int i = 0; i < M; i++)
        for (int j = 0; j < M; j++)
        {
            Hessian(i,j) = 3 * ((l.lambda[0][i] * l.lambda[0][j]) + l.lambda[1][i * l.M + j]) * (1 - (epsilon/3));
            for (int k = 0; k < M; k++)
                for (int a = 0; a < M; a++)
                    Hessian(i,j) -= 2 * ((l.lambda[0][i]/l.lambda[0][j]) + (l.lambda[0][j]/l.lambda[0][i]))
                        * (l.lambda[0][k] * l.lambda[0][a] * l.lambda[1][k * l.M + a] - epsilon * epsilon); 
        }
    SelfAdjointEigenSolver<MatrixXd> es(Hessian);

    eigenV = es.eigenvalues();
    U = es.eigenvectors();
}

void ObsParam::printObs(std::ofstream &file)
{
    file << '(';
    if (solution == FIXED)
        file << 1 << ',';
    if (solution == ZEROFIXED)
        file << 0 << ',';
    if (solution == INSUF)
        file << 2 << ',';
    file << n << ',' << r << ',' << nT << ',' << alpha << ','
         << fNL << ')' << std::endl;
}

ObsParam::ObsParam() { }

ObsParam::ObsParam(Lambda &l, enum solutionType sol)
{

    solution = sol;


    H = l.fields[l.M];
    tau = l.fields[l.M + 1]; k = 1/tau; 
    Hess(l);

    spectralIndex(l);
    spectralTensorIndex(l);
    spectralTilt(l);
    spectralAlpha(l);
    spectralfNL(l);
}

