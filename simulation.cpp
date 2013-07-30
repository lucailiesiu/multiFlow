


void param::initialize(double **initialParam, const double **lambdaLimitsUp, const double **lambdaLimitsDown)
{
    
    for (int i = 0; i < M; i++)
    {
        for (int j = 0; j < pow(M, i+1); j++)
        {
            initialParam[i][j] = (lambdaLimitsUp[i][j] - lambdaLimitsDown[i][j]) * ((float)rand()/(float)RAND_MAX) + lambdaLimitsDown[i][j];   
    }
}
