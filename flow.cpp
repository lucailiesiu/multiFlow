#include </u/liliesiu/MFC_code/Eigen/Core>
#include </u/liliesiu/MFC_code/Eigen/Dense>
#include </u/liliesiu/MFC_code/Eigen/Eigenvalues>
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>

#include "flow.h"
#include "tools.h"
#include "lambda.h"
#include "ObsParameters.h"


using namespace Eigen;
using std::vector;

Flow::Flow() { dN = DELTAN; nrIter = 0;}

Flow::Flow(int giveM, int giveK,
           double **lambdaLimitsUp, double **lambdaLimitsDown)
{
    dN = DELTAN;
    nrIter = 0;
    M = giveM;
    kTrunc = giveK;
    Lambda *lInit = new Lambda (M, kTrunc, lambdaLimitsUp, lambdaLimitsDown);
    initial = *lInit;
}

Flow::Flow(int M, int kTrunc, Lambda l)
{
    dN= DELTAN;
    nrIter = 0;
    this->M = M; 
    this->kTrunc = kTrunc;
    initial = Lambda(l);
}


void Flow::toolsAdd(int *a, int position)
{
    if (a[position] < M - 1)
    {
        a[position] += 1;
    }
    else
    {
        a[position] = 0;
        toolsAdd(a, position - 1);
    }
   
}

void Flow::f(Lambda l, Lambda &dl)
{
   
    Lambda *zeroLam = new Lambda(M, kTrunc); 
    dl = *zeroLam;
    double epsilon = l.epsilon();

    for (int i = 0; i < M; i++)
    {
        for(int j = 0; j < M; j++)
        {
            dl.lambda[0][i] += l.lambda[0][j] * l.lambda[1][M * i + j];
        }
        dl.lambda[0][i] -= l.lambda[0][i] * epsilon;
    }
    int *a;
    a = new int[1];
    for (int i = 1; i < kTrunc; i++)
    {
        // Initialize a
        delete[] a;
        a = new int[i + 1];
        for (int j = 0; j < i + 1; j++) { a[j] = 0; }
        for (int j = 0; j < pow(M, i+1); j++)
        {
            // First term in multi-field ODE
            if (i != kTrunc - 1)
            {
                for (int k = 0; k < M; k++)
                {
                    dl.lambda[i][j] += l.lambda[i + 1][M * j + k];
                }
            }
            //Second term in multi-field ODE
                      
            for (int ind = 0; ind < M; ind++)
                for(int k = 2; k < i + 1; k++)
                    dl.lambda[i][j] += l.lambda[i][j] * l.lambda[1][ind * M + a[k]] * l.lambda[0][ind]/l.lambda[0][a[k]];
            dl.lambda[i][j] -= i * l.lambda[i][j] * epsilon;
            if (j!= pow(M, i + 1) - 1) toolsAdd(a,i);
        }
    }
    //ODE for the fields 
    for (int i = 0; i < M; i++)
    {
        dl.fields[i] = - FIELD_SCALE * l.lambda[0][i];
    }
    //ODE for H
    dl.fields[M] = epsilon*l.fields[M];
    //ODE for Tau
    dl.fields[M + 1] = ((double)1.00)/(INITIAL_SCALE * exp(dN * nrIter) * l.fields[M]);
}


enum solutionType Flow::checkMatrix()
{
    Lambda last = memorized.at(nrIter);
    Lambda blast  = memorized.at(nrIter - 1);
    ObsParam pastObs(blast, NOTFOUND); 
    ObsParam currentObs(last, NOTFOUND);
    MatrixXd dUdt = (currentObs.U - pastObs.U) * (1/dN);
    MatrixXd timeVar = (3.0/2.0) * (currentObs.U).transpose() * dUdt;
    enum solutionType aux = SLOWMOVING; 
    for (int i = 0; i < M; i++)
        for (int j = 0; j < M; j++)
            if ((i != j && timeVar(i, j) > 0.05) || (i == j && timeVar(i,j) > 0.05))
                aux = NOTFOUND;
    return aux; 
}

enum solutionType Flow::checkFixedLambda()
{
    MatrixXd meanMatrix = MatrixXd::Zero(M, M);
    MatrixXd stddevMatrix  = MatrixXd::Zero(M, M);
    double min =  MINSAVE/dN;
    for (int i = nrIter - (int)min  - 1; i < nrIter; i++)
    {
        meanMatrix += memorized.at(i).lambdaMatrix();
    }
    meanMatrix = meanMatrix * (1/( min));
    for (int i = nrIter - (int)min  - 1; i < nrIter; i++)
    {
        stddevMatrix += ((memorized.at(i).lambdaMatrix()
                          - meanMatrix).array().abs()).matrix();
    }
    stddevMatrix = stddevMatrix * (1/((double) min));
    std::cout << std::endl <<  stddevMatrix << std::endl;
    for (int i = 0; i < M; i++)
        for (int j = 0; j < M; j++)
            if (stddevMatrix(i,j) > 0.001)
                return NOTFOUND;
    return SLOWMOVING; 
}

enum solutionType Flow::checkSol()
{
    double epsilon = memorizedEps.back();
    double min = MINSAVE/dN;
    // We allow a minimum number of E-folds 
    if (nrIter < (int) min)
    { 
        return NOTFOUND;
    }

    // When we reach the minimum number of E-folds
    // we compute the mean and stddev.
    if (nrIter == (int) min)
    {
        mean = 0; stddev = 0;
        for (std::vector<double>::iterator it = memorizedEps.begin(); it != memorizedEps.end(); it++)
            mean += (*it);
        mean = mean / (min);
        for (std::vector<double>::iterator it = memorizedEps.begin(); it != memorizedEps.end(); it++)
            stddev += std::abs((*it) - mean);
        stddev = stddev / (min);

        // Insert the means and stddevs in a vector for tracking convenience
        for (int i =0; i <= (int) min; i++)
        {
            means.push_back(mean);
            stddevs.push_back(stddev);
        }
        return NOTFOUND; 
    }

    // If inflation stops before the minimum number of E-folds
    if ((epsilon > 1) && (nrIter < (CLASSICEFOLDS/dN)))
    {
        return INSUF;
    }

    // If we have exceded the max number of e-folds and found
    // no convergent behaviour.
    if (nrIter > (MAXSAVE/(dN)))
        {
            return NONTRIVIAL;
        }

    // Recompute mean and stddev
    mean += (epsilon - means.at(nrIter - (int)min - 1))/(min);
    stddev += (std::abs(epsilon - mean) - stddevs.at(nrIter - (int)min - 1))/(min);
    means.push_back(mean);
    stddevs.push_back(stddev);
    // Informative print statement for means and stddev
    if (nrIter % PRINTPARSE == 0) 
        std::cout << mean << " " << stddev << std::endl;
    // Condition for fixed point
    if (nrIter % PARSEFIXED == 0)
    {
        if (stddev < 0.0005)
        {
            enum solutionType matrixCheck= checkMatrix();
            enum solutionType lambdaCheck = checkFixedLambda();
            if (matrixCheck == SLOWMOVING && lambdaCheck)
            {
                // Condition for fixed point at 0
                if (mean < 0.005)
                    return ZEROFIXED;
                return FIXED;
            }
        }
    }
    // If no conditions have been satisfied,
    // behaviour still not found
    return NOTFOUND;
}



void Flow::rk45(Lambda &RKresult)
{
    Lambda l = memorized.back();
    Lambda a1, a2, a3, a4, A2, A3;
    f(l, a1);
    f(l.lSum(a1.lProd(dN/4)), a2);
    f(l.lSum(a2.lProd(dN/2)), a3);
    f(l.lSum(a3.lProd(dN)), a4);
        
    A2 = a2.lProd(2.00);
    A3 = a3.lProd(2.00);

    Lambda sum = a1.lSum(A2.lSum(A3.lSum(a4)));
    RKresult = l.lSum(sum.lProd(dN/6));
}


void Flow::generate()
{
    Lambda save;
    enum solutionType solType = NOTFOUND;
    std::cout << "Solving ODE system for:" << std::endl;
    memorized.push_back(initial);
    memorizedEps.push_back(initial.epsilon());
    while(solType == NOTFOUND)
    {
        rk45(save);
        memorized.push_back(save);
        memorizedEps.push_back(save.epsilon());
        solType = checkSol(); 
        nrIter++;
        if (nrIter % PRINTPARSE == 0)
        {
            std::cout << nrIter << std::endl;
        }
    }
    (*this).observable = ObsParam(save, solType);
}


void Flow::printFlow(std::ofstream &file)
{
    int count = 0; 
    observable.printObs(file);
    
    for (std::vector<Lambda>::iterator it = memorized.begin();
         it != memorized.end(); it++)
    {
        if (count % PARSE == 0)
        {
            (*it).printLambda(file);
            (*it).printFields(file);
        }
        count ++;
    }
}

