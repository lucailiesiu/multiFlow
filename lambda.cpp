#include </u/liliesiu/MFC_code/Eigen/Dense>
#include </u/liliesiu/MFC_code/Eigen/Eigenvalues>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <time.h>
#include <fstream>

#include "tools.h"
#include "lambda.h"

using namespace Eigen;

Lambda::Lambda() {M = 0; kTrunc = 0;}

Lambda::Lambda(const Lambda& l)
{
    M = l.M; kTrunc = l.kTrunc;
    lambda = new double*[kTrunc];
    for (int i = 0; i < kTrunc; i++)
    {
        lambda[i] = new double[(int)pow(M,i+1)];
        for (int j = 0; j < (pow(M, i+1)); j++)
        {
                lambda[i][j] = l.lambda[i][j];
        }
    
    }
    fields = new double[M + 5];
    for (int i = 0; i < M + 5; i++)
        fields[i] = l.fields[i];
}


Lambda::Lambda(int giveM, int giveK, double **giveLambda, double *giveFields)
{
    (*this).M = giveM;
    (*this).kTrunc = giveK;
    lambda = new double*[giveK];
    for (int i = 0; i < giveK; i++)
    {
        lambda[i] = new double[(int)pow(giveM, i+1)];
        for (int j = 0; j < pow(giveM, i+1); j++)
        {

            lambda[i][j] = giveLambda[i][j];
        }
    }
    fields = new double[giveM + 5];
    for (int i = 0; i < giveM + 5; i++)
        fields[i] = giveFields[i];
    
}
 
 
Lambda::Lambda(int giveM, int giveK, double **lambdaLimitsUp,
                       double **lambdaLimitsDown)
{

    srand (time(NULL));
    M = giveM;
    kTrunc = giveK;
    lambda = new double*[kTrunc];
    for (int i = 0; i < giveK; i++)
    {
        lambda[i] = new double[(int)pow(M, i+1)];
        int* a;
        a = new int[i + 1];
        for (int j = 0; j < i + 1; j++) a[j] = 0;
        while(a[0] != M - 1)
        {
            int j = toolsIndex(a, i + 1);
            lambda[i][j] = (lambdaLimitsUp[i][j] - lambdaLimitsDown[i][j])
                * ((float)rand()/(float)RAND_MAX) + lambdaLimitsDown[i][j];
            generatePerm(a, i + 1, j); 
            toolsIncrease(a, i);
        }
        int j = toolsIndex(a, i + 1);
        lambda[i][j] =  (lambdaLimitsUp[i][j] - lambdaLimitsDown[i][j])
            * ((float)rand()/(float)RAND_MAX) + lambdaLimitsDown[i][j];
        generatePerm(a, i + 1, j);
    }
    double *giveFields;
    giveFields = new double[giveM + 5];
    for (int i = 0; i < M + 5; i++) giveFields[i] = 0;
    giveFields[giveM] = 1.00;
    fields = giveFields;
}


Lambda::Lambda(int giveM, int giveK)
{
    M = giveM;
    kTrunc = giveK;
    lambda = new double*[giveK];
    for (int i = 0; i < giveK; i++)
    {
        lambda[i] = new double[(int)pow(giveM, i+1)];
        for (int j = 0; j < pow(giveM, i+1); j++)
        {
            lambda[i][j] = 0;
        }
    }
    double *giveFields; 
    giveFields = new double[giveM + 5];
    for (int i = 0; i < M + 5; i++) giveFields[i] = 0;
    giveFields[giveM] = 1.00;
    fields = giveFields;
}


Lambda::~Lambda()
{
    if (M == 0 || kTrunc == 0)
        return;
    for (int i=0; i < kTrunc; i++)
    {
        delete [] lambda[i];
    }
    delete[] lambda;
    delete[] fields;
    lambda = NULL;
    fields = NULL;
}


Lambda &Lambda::operator= (const Lambda& l)
{

    if (this != &l)
    {

        this->~Lambda();
        kTrunc = l.kTrunc;
        M = l.M;
        lambda = new double*[kTrunc];
        for (int i = 0; i < kTrunc; i++)
        {
            lambda[i] = new double[(int)pow(M, i+1)];
            for (int j = 0; j < pow(M, i+1); j++)
            {

                lambda[i][j] = l.lambda[i][j];
            }

        }
        
        fields = new double[M + 5];
        for (int i = 0; i < M + 5; i++)
            fields[i] = l.fields[i];
    }
    return *this;
}

void Lambda::toolsIncrease(int* a, int position)
{
     if (a[position]  < M - 1)
    {
        a[position] += 1;
    } else
    {
        toolsIncrease(a, position - 1);
        a[position] = a[position - 1];
    }
}

int Lambda::toolsIndex(const int* a, const int number)
{
    int ind = 0; 
    for(int i = 0; i < number; i++)
    {
        ind += int(a[i] * pow(M, number - i - 1));
    }
    return ind; 
}

double Lambda::epsilon()
{
    double eps = 0;
    for (int i = 0; i < M; i++)
        eps += pow(lambda[0][i], 2);
    return eps;
}

void Lambda::generatePerm(const int *a, const int number, const int index)
{
    int *b;
    b = new int[number];
    for (int i = 0; i < number; i++) b[i] = a[i];
    do {
        lambda[number - 1][toolsIndex(b, number)] = lambda[number - 1][index];
    } while (std::next_permutation(b, b + number));
    
}

MatrixXd Lambda::lambdaMatrix()
{
    int M = (*this).M;
    MatrixXd H(M, M);
    for (int i = 0; i < M; i++)
        for (int j = 0; j < M; j++)
            H(i,j) = (*this).lambda[1][M * i + j];
    return H; 
}

Lambda Lambda::lSum(Lambda l)
{
    if ((M = l.M) && (kTrunc == l.kTrunc))
    {
        Lambda sum = Lambda(l.M, l.kTrunc);
        for (int i = 0; i < kTrunc; i++)
        {
            for (int j = 0; j < pow(M, i+1); j++)
            {
                sum.lambda[i][j] = l.lambda[i][j] + lambda[i][j];
            }
        }
        for (int i = 0; i < M + 5; i++)
        {
            sum.fields[i] = l.fields[i] + fields[i];
        }
        return sum; 
    }
    return Lambda(0,0);
}

Lambda Lambda::lProd(double scalar)
{
    Lambda prod = Lambda(M, kTrunc);
    for (int i = 0; i < kTrunc; i++)
    {
        for (int j = 0; j < pow(M, i+1); j++)
        {
            prod.lambda[i][j] = scalar * (*this).lambda[i][j];
        }
    }
    for (int i = 0; i < M + 4; i++)
    {
        prod.fields[i] = scalar * (*this).fields[i];
    }
    return prod; 
}


void Lambda::printLambda(std::ofstream &file)
{
    file << "[";
    for (int i = 0; i < kTrunc; i++)
    {
        file << "[";
        for (int j = 0; j < pow(M, i+1); j++)
        {
            if (j < pow(M, i + 1) - 1)
            {
                file << lambda[i][j];
                file << ", ";
             }
            else
                file << lambda[i][j];
        }
        if (i < kTrunc - 1)
            file << "], ";
        else
            file << "]";
    }
    file << "]" << std::endl;
}

void Lambda::printEpsilon(std::ofstream &file)
{
    file << "[";
    for (int i = 0; i < M; i++)
    {
        if (i < M)
        {
            file << lambda[0][i];
            file << ", ";
        }
        else
            file << lambda[0][i];
    }
    file << "]" << std::endl;
}

void Lambda::printFields(std::ofstream &file)
{
    file << "[";
    for (int i = 0; i < M + 4; i++)
    {
        if (i < M + 3)
        {
            file << fields[i];
            file << ", ";
        }
        else
            file << fields[i];
    }
    file << "]" << std::endl; 
}
