double my_log(double x)
{
    double a=(1.0+x)/2.0;
    double b=sqrt(x);
    while(1){
        a=(a+b)/2.0;
        b=sqrt(a*b);
        if(fabs(a-b)<0.001) break;
    }
    return 2.0*(x-1.0)/(a+b);
}

double Priors(Par *p)
{
    double vp=gsl_ran_beta_pdf(p->p,8,2);
    double va=gsl_ran_lognormal_pdf(p->alpha,log(10),0.2);
    double vg=gsl_ran_lognormal_pdf(p->gamma,log(2),0.2);
    return my_log(vp)+my_log(va)+my_log(vg);
}

double likelihood(Par *p)
{
    int i,j;
    double llhood=0.0,tmp=0.0,t1;
    double *mu,*It;
    mu=malloc(datal*sizeof(double));
    It=malloc(datal*sizeof(double));
    for (i=0;i<Ndelay;i++)
    {
        mu[i]=p->Mun[i];
        It[i]=p->In[i];
    }
    for(i=Ndelay;i<datal;i++)
    {
        tmp=0.0;
        for(j=1;j<=Ndelay;j++) tmp+=mu[i-j]*(exp(-(j-1)*p->gamma)-exp(-j*p->gamma))/p->gamma+It[i-j]*(exp(-(j-1)*p->delta)-exp(-j*p->delta))/p->delta;
         for(j=1;j<=Ndelay;j++){
             t1=p->alpha/(p->delta-p->alpha)*(  (exp(-(j-1)*p->alpha)-exp(-j*p->alpha))/p->alpha  -  (exp(-(j-1)*p->delta)-exp(-j*p->delta))/p->delta );
            int nl=(i-j)/IntervalL ;
            tmp+=p->lambda[nl]*t1;
         }        
        if (datas[i].sfh){
            It[i]=p->beta*p->p*datas[i].At*tmp;
        }
        else{
            It[i]=p->beta*datas[i].At*tmp;
        }
        tmp=0.0;
        for (j=1;j<=Ndelay;j++) tmp+=It[i-j]*(1-(exp(-(j-1)*p->delta)-exp(-j*p->delta))/p->delta);
        for (j=1;j<=Ndelay;j++) {
             t1=1.0+(p->delta*(exp(-(j-1)*p->alpha)-exp(-j*p->alpha))/p->alpha - p->alpha*(exp(-(j-1)*p->delta)-exp(-j*p->delta))/p->delta )/(p->alpha-p->delta);
            int nl=(i-j)/IntervalL ;
            tmp+=p->lambda[nl]*t1; 
        }
        mu[i]=tmp;
    }
    for(i=0;i<datal;i++)
    {
        llhood+=datas[i].Dt*log(mu[i])-mu[i];
    }
    llhood-=SumLogGamma;
    llhood+=Priors(p);
    free(mu);
    free(It);
    return llhood;
}

void  Printlikelihood(Par *p)
{
    int i,j;
    double tmp=0.0,t1;
    double *mu,*It;
    mu=malloc(datal*sizeof(double));
    It=malloc(datal*sizeof(double));
    for (i=0;i<Ndelay;i++)
    {
        mu[i]=p->Mun[i];
        It[i]=p->In[i];
    }
    for(i=Ndelay;i<datal;i++)
    {
        tmp=0.0;
        for(j=1;j<=Ndelay;j++) tmp+=mu[i-j]*(exp(-(j-1)*p->gamma)-exp(-j*p->gamma))/p->gamma+It[i-j]*(exp(-(j-1)*p->delta)-exp(-j*p->delta))/p->delta;
         for(j=1;j<=Ndelay;j++){
             t1=p->alpha/(p->delta-p->alpha)*(  (exp(-(j-1)*p->alpha)-exp(-j*p->alpha))/p->alpha  -  (exp(-(j-1)*p->delta)-exp(-j*p->delta))/p->delta );
            int nl=(i-j)/IntervalL ;
            tmp+=p->lambda[nl]*t1;
         }        
        if (datas[i].sfh){
            It[i]=p->beta*p->p*datas[i].At*tmp;
        }
        else{
            It[i]=p->beta*datas[i].At*tmp;
        }
        tmp=0.0;
        for (j=1;j<=Ndelay;j++) tmp+=It[i-j]*(1-(exp(-(j-1)*p->delta)-exp(-j*p->delta))/p->delta);
        for (j=1;j<=Ndelay;j++) {
             t1=1.0+(p->delta*(exp(-(j-1)*p->alpha)-exp(-j*p->alpha))/p->alpha - p->alpha*(exp(-(j-1)*p->delta)-exp(-j*p->delta))/p->delta )/(p->alpha-p->delta);
            int nl=(i-j)/IntervalL ;
            tmp+=p->lambda[nl]*t1; 
        }
        mu[i]=tmp;
    }
    for(i=0;i<datal;i++)
    {
    fprintf(Fmu,"%lf \t", mu[i]);
    fprintf(Fit,"%lf \t", It[i]);
    }
    fprintf(Fmu,"\n" );
    fprintf(Fit,"\n" );
    free(mu);
    free(It);
}