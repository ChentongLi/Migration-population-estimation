#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_sf.h>
#include "header/init.h"
#include "header/likelihood.h"
#include "header/MCMC.h"

int main() 
{
    init();
    MCMC();
    return 0;
}
