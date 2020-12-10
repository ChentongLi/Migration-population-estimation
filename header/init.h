#define Ndelay 12
#define IntervalL 3
#define Nlambda 44
FILE *Fmu,*Fit;

typedef struct Par
{
    double lambda[Nlambda];
    double In[Ndelay],Mun[Ndelay]; 
    double alpha,beta,gamma,delta,p;
}Par;

typedef struct  DataSet
{
    int index;
    int year;
    short  month;
    double At;
    double Dt;
    short  sfh;
}DataSet;
DataSet *datas;
int datal;
double SumLogGamma;

void init()
{
	FILE *fp;
	fp=fopen("data/data.csv","r");
    if (fp==NULL) {
        printf("No that files!\n");
        exit(1);
    }
    int len=128,i=0;
    datas=malloc(sizeof(DataSet)*len);
    while(!feof(fp)){
        fscanf(fp,"%d ,%d,%hd,%lf,%lf,%hd\n",&datas[i].index,&datas[i].year,&datas[i].month,
        &datas[i].Dt,&datas[i].At,&datas[i].sfh);
        datas[i].At*=1e-7;
        i++;
        if (i>=len) {
            len+=128;
            datas=realloc(datas,sizeof(DataSet)*len);
        }
    }
    datal=i;
    SumLogGamma=0;
    for (i=0;i<datal;i++){
        SumLogGamma+=gsl_sf_lngamma(datas[i].Dt+1);
    }
	fclose(fp);
}

void InitPar(Par *p)
{
    int i;
    for (i=0;i<Nlambda;i++) p->lambda[i]=200.0;
    for(i=0;i<Ndelay;i++){
        if (i==0) p->In[i]=datas[i].Dt;
        else p->In[i]=datas[i-1].Dt;
        p->Mun[i]=datas[i].Dt;
    } 
    p->alpha=2.0;
    p->beta=0.4;
    p->delta=4.0/3.0;
    p->gamma=1;
    p->p=0.6;
}

