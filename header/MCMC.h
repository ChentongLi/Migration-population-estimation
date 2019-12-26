#define MIN(x,y)  (((x)<(y)) ?  (x) : (y))
#define MAX(x,y)  (((x)>(y)) ?  (x) : (y))
#define MAXMCTIMES 1e8

void SetValue (Par *out, Par *in)
{
    int i;
    for (i=0;i< Nlambda;i++) out->lambda[i]=in->lambda[i];
    for (i=0;i<Ndelay;i++){
        out->In[i]=in->In[i];
        out->Mun[i]=in->Mun[i];
    }
    out->alpha=in->alpha;
    out->beta=in->beta;
    out->gamma=in->gamma;
    out->delta=in->delta;
    out->p=in->p;
}

void RandPar(Par *p, gsl_rng *r)
{
    // gsl_ran_lognormal(r,log(mu),10.0/mu);
    // gsl_ran_beta(r,x*n,(1-x)*n);
    int i;
    double mu;
    for (i=0;i<Nlambda;i++) {
        mu=p->lambda[i];
        if (mu<100)
            p->lambda[i]=gsl_ran_lognormal(r,log(mu),1.0/mu);
        else
            p->lambda[i]=gsl_ran_lognormal(r,log(mu),15.0/mu);
    }
    for (i=0;i<Ndelay;i++) {
        mu=p->In[i];
        if (mu<100)
            p->In[i]=gsl_ran_lognormal(r,log(mu),1.0/mu);
        else
            p->In[i]=gsl_ran_lognormal(r,log(mu),10.0/mu);
    }
  /*  
    for (i=0;i<Ndelay;i++) {
        mu=p->Mun[i];
        if (mu<100)
            p->Mun[i]=gsl_ran_lognormal(r,log(mu),1.0/mu);
        else
            p->Mun[i]=gsl_ran_lognormal(r,log(mu),10.0/mu);
    }
*/
    double u;
    u=gsl_rng_uniform(r)-0.5;
    p->alpha*=exp(u*0.2);
    u=gsl_rng_uniform(r)-0.5;
    p->beta*=exp(u);
    u=gsl_rng_uniform(r)-0.5;
    p->gamma*=exp(u*0.2);
    u=gsl_rng_uniform(r)-0.5;
    p->delta*=exp(u);
    p->p= gsl_ran_beta(r,p->p*200.0,(1-p->p)*200.0);
    if (p->p > 0.999) p->p=0.999;
}

void MCMC()
{
    const gsl_rng_type *T;
    gsl_rng *r;   
    T=gsl_rng_ranlxs0; 
    gsl_rng_default_seed = ((unsigned long)(time(NULL))); 
    r=gsl_rng_alloc(T); 

    int ncount=0,i;
    int id=getpid();
    double ALPHA,llhood,R;
    double pllhood=-1e50;
    char fname[30],fname1[30],fname2[30],fname3[30];

    Par *para,*tp;
    para=malloc(sizeof(Par));
    tp=malloc(sizeof(Par));
    InitPar(para);

    sprintf(fname,"result/%d_par.txt",id); 
    sprintf(fname1,"result/%d_lamb.txt",id);
    sprintf(fname2,"result/%d_Mu.txt",id);
    sprintf(fname3,"result/%d_It.txt",id);
    FILE *fp,*fl;
    fp=fopen(fname,"w+");
    fl=fopen(fname1,"w+");
    Fmu=fopen(fname2,"w+");
    Fit=fopen(fname3,"w+");
    printf("The initilization is over, then go to the calculation process, %d!\n",id);
    
    while (ncount<MAXMCTIMES){

        SetValue(tp,para);
        RandPar(para,r);
       llhood=likelihood(para);

       ALPHA=MIN(0,llhood-pllhood);
       R=log(gsl_rng_uniform(r));

      //  printf("%e, %e,%lf\n",llhood,pllhood,R);
        if (R<ALPHA){
            pllhood=llhood;
            printf("%lf\n",llhood);
            for (i=0;i<Nlambda;i++)  fprintf(fl,"%lf \t", para->lambda[i]);
            fprintf(fl,"\n");
            fprintf(fp,"%e\t%e\t%e\t%e\t%e\t%lf\n",para->alpha,para->beta,para->gamma,para->delta,para->p,pllhood);
            Printlikelihood(para);
        }
        else{
            SetValue(para,tp);
        }
        ncount++;
    }
    fclose(fp);
    fclose(fl);
    fclose(Fmu);
    fclose(Fit);
    gsl_rng_free(r); 
}

