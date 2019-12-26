function prcc

runs=1000;

beta_LHS=LHS_Call(4.48,4.59,4.7,0.02,runs,'unif')*1e-9;
delta_LHS=LHS_Call(7.7,8.0,8.3,0.04,runs,'unif')*1e-2;
gamma_LHS=LHS_Call(3,4.08,5.43,0.04,runs,'unif');
alpha_LHS=LHS_Call(7.9,11.5,18.4,0.1,runs,'unif');

y=[];
for i=1:runs
    y(i)=MeanInf(beta_LHS(i),delta_LHS(i),gamma_LHS(i),alpha_LHS(i));
end

MI=y';
PM=[beta_LHS,delta_LHS,gamma_LHS,alpha_LHS];
rank_PM=ranking(PM);
rank_MI=ranking(MI);
N=size(PM,2);

PRCC=[];
for j=1:N
    PMt=rank_PM;
    rank_x=PMt(:,j);
    PMt(:,j)=[];
    s1=regstats(rank_MI,PMt,'linear','r');
    y_res=s1.r;
    s2=regstats(rank_x,PMt,'linear','r');
    x_res=s2.r;
    PRCC=[PRCC,corr(y_res,x_res)];
end
PRCC
bar(PRCC)
title('Sensitivity of MI with respect to model parameters')
set(gca,'xTick',1:4,'xTicklabel',{'\beta','\delta','\gamma','\alpha'})

function r=ranking(x)

[a,b]=size(x);
for j=1:b
   [s,i]=sort(x(:,j));
   r(i,j)=[1:a]';
end

function MI=MeanInf(beta,delta,gamma,alpha)
    datal=130;
    Ndelay=6;
    IntervalL=4;
    lambda=[273.0372835    92.91015964    761.1723267    1424.629764    44.28150733    1538.622325    1219.587228    3.684478042    517.2987882    1167.082123    124.7830456    687.6729476    1004.239344    69.79360306    803.1518263    1237.625027    81.09670832    502.7889744    1455.48373    84.17674003    664.9920722    1288.855116    378.7569767    548.1260017    1514.938163    78.34451308    615.5629799    1298.868112    76.8835498    612.0238194    1425.070885    710.108228    603.0603769];
    Mun=[2449    2530    2747    2964    2971    3245];
    In=[2529.625982    2627.855645    2207.895709    2383.901297    2564.407134    2821.792687];
    At=[ 31243457.67    31248494.33    31253531    31258567.67    31263604.33    31268641    31273677.67    31278714.33    31283751    31288787.67    31293824.33    31298861    31305492    31312123    31318754    31325385    31332016    31338647    31345278    31351909    31358540    31365171    31371802    31378433    31383904.58    31389376.17    31394847.75    31400319.33    31405790.92    31411262.5    31416734.08    31422205.67    31427677.25    31433148.83    31438620.42    31444092    31433959.5    31423827    31413694.5    31403562    31393429.5    31383297    31373164.5    31363032    31352899.5    31342767    31332634.5    31322502    31315090    31307678    31300266    31292854    31285442    31278030    31270618    31263206    31255794    31248382    31240970    31233558    31217119    31200680    31184241    31167802    31151363    31134924    31118485    31102046    31085607    31069168    31052729    31036290    31018919.92    31001549.83    30984179.75    30966809.67    30949439.58    30932069.5    30914699.42    30897329.33    30879959.25    30862589.17    30845219.08    30827849    30808926.83    30790004.67    30771082.5    30752160.33    30733238.17    30714316    30695393.83    30676471.67    30657549.5    30638627.33    30619705.17    30600783    30578800.83    30556818.67    30534836.5    30512854.33    30490872.17    30468890    30446907.83    30424925.67    30402943.5    30380961.33    30358979.17    30336997    30327109.25    30317221.5    30307333.75    30297446    30287558.25    30277670.5    30267782.75    30257895    30248007.25    30238119.5    30228231.75    30218344    30208344    30198344    30188344    30178344    30168344    30158344    30148344    30138344    30128344    30118344];
    mu=zeros(1,datal);
    It=zeros(1,datal);
    for i=1:Ndelay
        mu(i)= Mun(i);
        It(i)= In(i);
    end
    for i=Ndelay+1:datal
        tmp=0;
        for j=1:Ndelay   tmp=tmp+mu(i-j)*(exp(-(j-1)*gamma)-exp(-j*gamma))/gamma+It(i-j)*(exp(-(j-1)*delta)-exp(-j*delta))/delta;
        end
        for j=1:Ndelay
            t1=alpha/(delta-alpha)*( (exp(-(j-1)* alpha)-exp(-j* alpha))/ alpha  -  (exp(-(j-1)* delta)-exp(-j* delta))/ delta );
            nl=floor((i-j)/IntervalL)+1;
            tmp=tmp+lambda(nl)*t1;
        end
        It(i)= beta*At(i)*tmp;
        tmp=0;
        for j=1:Ndelay
            tmp+=It(i-j)*(1-(exp(-(j-1)* delta)-exp(-j* delta))/ delta);
        end
        for j=1:Ndelay
             t1=1.0+(delta*(exp(-(j-1)* alpha)-exp(-j* alpha))/ alpha -  alpha*(exp(-(j-1)* delta)-exp(-j* delta))/ delta )/( alpha- delta);
            nl=floor((i-j)/IntervalL)+1;
            tmp=tmp + lambda(nl)*t1;
        end
        mu(i)=tmp;
    end
    MI=sum(It)/datal;

