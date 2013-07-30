
#include <sys/stat.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include <cmath>

#include "flow.h"
#include "lambda.h"
#include "ObsParameters.h"

bool fileExists(const std::string& filename)
{
    struct stat buf;
    if (stat(filename.c_str(), &buf) != -1)
    {
        return true;
    }
    return false;
}

double** up(int M, int Ntrunc)
{
    double **a = new double*[Ntrunc];
    for (int i = 0; i < Ntrunc - 1; i++)
    {
        a[i] = new double[(int) pow(M, i+1)]; 
    }
    for(int i = 0; i < M; i++)
    {
        a[0][i] = 0.8;
    }
    for (int i = 1; i < Ntrunc - 1; i++)
        for(int j = 0; j < pow(M, i+1); j++)
        {
            a[i][j] = 0.1 * pow(0.1, i-1);
        }
    return a; 
}

double** down(int M, int Ntrunc)
{
    double **a = new double*[Ntrunc];
    for (int i = 0; i < Ntrunc - 1; i++)
    {
        a[i] = new double[(int) pow(M, i+1)];
            }
    for(int i = 0; i < M; i++)
    {
        a[0][i] = -0.8;
    }
    for (int i = 1; i < Ntrunc - 1; i++)
        for(int j = 0; j < pow(M, i+1); j++)
        {
            a[i][j] = -0.1 * pow(0.1, i-1);
        }
    return a;
}


int main()
{
   
    int M = 4;
    int Ntrunc = 6;
    double** lambdaLimitsUp = up(M, Ntrunc + 1);
    double** lambdaLimitsDown = down(M, Ntrunc + 1);
    std::cout << "ac";
    Lambda *l = new Lambda(M, Ntrunc, lambdaLimitsUp, lambdaLimitsDown);
    std::cout << "Sum is:" << std::endl;
    
    std::cout << "Product is:" << std::endl;
   
    std::cout << "Epsilon is " << (*l).epsilon() << std::endl;
    
    ObsParam obs(*l, FIXED);
    
    std::cout << "Generating flow"; 
    Flow flow(M, Ntrunc, *l);
    flow.generate();
    std::cout << "Printing data to file" << std::endl;

    int count = 0;
    std::string fileName = "data/evolution";
    std::string index = "0.out";
    std::ostringstream countString;
    countString << fileName << index;
    const char* fileN;
    
    while (fileExists(countString.str()))
    {
        count++;
        countString.str("");
        countString.clear();
        countString << fileName << count << ".out";
    }
   
    std::ofstream myfile;
    std::string fileNom = countString.str();
    std::cout << "Printed to file: " << fileNom << std::endl;
    fileN = fileNom.c_str();
    myfile.open (fileN);
    flow.printFlow(myfile);
    myfile.close();
    return 0;
}
