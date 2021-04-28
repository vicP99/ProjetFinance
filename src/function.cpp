#include "../include/header.hpp"


using namespace std;

vecteur operator +(const vecteur & v1,const vecteur & v2)
{
    int taille=v1.v.size();
    vecteur res(taille);
    for(int i=0;i<taille;i++)
    {
        res[i]=v1[i]+v2[i];
    }
    return res;
}

vecteur operator +(const double a,const vecteur & v)
{
    int taille=v.v.size();
    vecteur res(taille);
    for(int i=0;i<taille;i++)
    {
        res[i]=v[i]+a;
    }
    return res;
}

vecteur operator +(const vecteur & v, const double a)
{
    return a+v;
}
vecteur operator -(const vecteur & v1,const vecteur & v2)
{
    int taille=v1.v.size();
    vecteur res(taille);
    for(int i=0;i<taille;i++)
    {
        res[i]=v1[i]-v2[i];
    }
    return res;
}

vecteur operator -(const double a,const vecteur & v)
{
    int taille=v.v.size();
    vecteur res(taille);
    for(int i=0;i<taille;i++)
    {
        res[i]=v[i]-a;
    }
    return res;
}

vecteur operator -(const vecteur & v, const double a)
{
    return a-v;
}

vecteur operator *(double a ,const vecteur & v)
{
    int taille=v.v.size();
    vecteur res(taille);
    for(int i=0;i<taille;i++)
    {
        res[i]=a*v[i];
    }
    return res;
}

vecteur operator *(const vecteur & v, double a)
{
    return a*v;
}

ostream & operator <<( ostream & flux,const vecteur & v)
{
    int taille=v.v.size();
    flux<<"(";
    for(int i=0;i<taille-1;i++)
    {
        flux<<v[i]<<",";
    }
    flux<<v[taille-1]<<")\n";
    return flux;
}
ostream & operator <<( ostream & flux,const simulation & v)
{
    flux<<"La valeur est: "<<v.val<<" et l'interval de confiance à "<<v.niv<<"\% est: ["<<v.ICinf<<","<<v.ICsup<<"] avec une erreure de "<<v.val-v.ICinf<<endl;
    return flux;
}
vecteur exp(const vecteur& v){
    vecteur res(v.v.size());
    for(uint i=0;i<v.v.size();i++){
        res[i]=exp(v[i]);
    }
    return res;
}
vecteur Plus(const vecteur& v){
    vecteur res(v.v.size());
    for(uint i=0;i<v.v.size();i++){
        res[i]=max(v[i],0.);
    }
    return res;
}
vecteur max(const vecteur& v1,const vecteur& v2){
    vecteur res(v1.v.size());
    for(uint i=0;i<v1.v.size();i++){
        res[i]=max(v1[i],v2[i]);
    }
    return res;
}
vecteur f(const vecteur& v, double K){
    return Plus(v-K) - Plus(v);
} 





vecteur loi_unif(double a, double b, int nb_points){
    vecteur res(nb_points);
    for (int i=0;i<nb_points;i++){
        res[i]  = a+(b-a)*((double)(rand()))/((double)RAND_MAX);
    }
    return res;
}

 
vector<vecteur> normal_indep(int nbpoint,double sig1,double sig2,double mu1,double mu2)
{
    vecteur U=loi_unif(0,1,nbpoint);
    vecteur V=loi_unif(0,1,nbpoint);
    vecteur X(nbpoint);
    vecteur Y(nbpoint);
    for(int i=0;i<nbpoint;i++)
    {
        X[i]=mu1+sig1*sqrt(-2*log(U[i]))*cos(2*M_PI*V[i]);
        Y[i]=mu2+sig2*sqrt(-2*log(U[i]))*sin(2*M_PI*V[i]);
    }
    vector<vecteur> res(2);
    res[0]=X;
    res[1]=Y;
    return res;
}

double E(vecteur v)
{
    int taille=v.v.size();
    double sum=0;
    for (int i=0;i<taille;i++)
    {
        sum+=v[i];
    }
    return sum/((double)taille);
}

double V(vecteur v)
{
    double mu=E(v);
    int taille=v.v.size();
    double sum=0;
    for (int i=0;i<taille;i++)
    {
        sum+=pow((v[i]-mu),2);
    }
    return sum/((double)(taille-1));
}

double Cov(vecteur v1,vecteur v2)
{
    double mu1=E(v1);
    double mu2=E(v2);
    double sum=0;
    int taille=min(v1.v.size(),v2.v.size());
    for(int i=0;i<taille;i++)
    {
        sum+=(v1[i]-mu1)*(v2[i]-mu2);
    }
    return sum/((double)(taille-1));
}

double I(double x)
{
    double t;
    double res;
    if(x>0)
    {
        t=1/(1+0.2316419*x);
        res=(1/sqrt(2*M_PI))*exp(-x*x/2)*(0.319381530*t-0.356563782*pow(t,2)+1.781477937*pow(t,3)-1.821255978*pow(t,4)+1.330274429*pow(t,5));
        return 1.-res;
    }
    else
    {
        t=1/(1-0.2316419*x);
        res=(1/sqrt(2*M_PI))*exp(-x*x/2)*(0.319381530*t-0.356563782*pow(t,2)+1.781477937*pow(t,3)-1.821255978*pow(t,4)+1.330274429*pow(t,5));
        return res;
    }
}

double psi(double y,parametre par){
    double alphaModifie=par.alpha*par.S1*exp((pow(par.sig1,2))/2*-par.T+par.sig1*y);
    double betaModifie=par.beta*par.S2*exp((pow(par.sig2*par.rho,2))*(-par.T/2)+par.sig2*par.rho*y);
    double sigMod=par.sig2*sqrt((1-pow(par.rho,2))*par.T);
    return alphaModifie*I((log(alphaModifie/betaModifie)/sigMod+sigMod/2)) -betaModifie*I((log(alphaModifie/betaModifie)/sigMod-sigMod/2));
}
vecteur psi(const vecteur& v,parametre par){
    vecteur res(v.v.size());
    for (uint i=0;i<v.v.size();i++){
        res[i]=psi(v[i],par);
    }
    return res;
}
double psi_spread(double y,parametre par){
    double alphaModifie=par.alpha*par.S1*exp((pow(par.sig1,2))/2*-par.T+par.sig1*y)- par.K*exp(-par.r*par.T);
    double betaModifie=par.beta*par.S2*exp((pow(par.sig2*par.rho,2))*(-par.T/2)+par.sig2*par.rho*y);
    double sigMod=par.sig2*sqrt((1-pow(par.rho,2))*par.T);
    return alphaModifie*I((log(alphaModifie/betaModifie)/sigMod+sigMod/2)) -betaModifie*I((log(alphaModifie/betaModifie)/sigMod-sigMod/2));
}
double psi_controle(double y,parametre par){
    double betaModifie=par.alpha*par.S1*exp((pow(par.sig1,2))/2*-par.T+par.sig1*y);
    double alphaModifie=par.beta*par.S2*exp((pow(par.sig2*par.rho,2))*(-par.T/2)+par.sig2*par.rho*y)+par.K*exp(-par.r*par.T);
    double sigMod=par.sig1*sqrt((1-pow(par.rho,2))*par.T);
    return alphaModifie*I((log(alphaModifie/betaModifie)/sigMod+sigMod/2)) -betaModifie*I((log(alphaModifie/betaModifie)/sigMod-sigMod/2));
}
vecteur psi_spread(const vecteur& v,parametre par){
    vecteur res(v.v.size());
    for (uint i=0;i<v.v.size();i++){
        res[i]=psi_spread(v[i],par);
    }
    return res;
}
vecteur psi_controle(const vecteur& v,parametre par){
    vecteur res(v.v.size());
    for (uint i=0;i<v.v.size();i++){
        res[i]=psi_controle(v[i],par);
    }
    return res;
}
simulation echange_MC(parametre par){
    simulation res;
    vecteur S1,S2;
    vector<vecteur> inde=normal_indep(par.nb_simul,sqrt(par.T),sqrt(par.T));
    vecteur W2=par.rho*inde[0]+sqrt(1-par.rho*par.rho)*inde[1];
    vecteur W1=inde[0];
    S1=par.alpha*par.S1*exp((par.r-par.sig1*par.sig1/2.)*par.T+par.sig1*W1);
    S2=par.beta*par.S2*exp((par.r-par.sig2*par.sig2/2.)*par.T+par.sig2*W2);
    vecteur Splus=Plus(S1-S2);
    res.val=exp(-par.r*par.T)*E(Splus);
    double var=exp(-2*par.r*par.T)*V(Splus);
    //cout<<"V(X)="<<var<<endl;
    res.ICinf=res.val-par.quant*sqrt(var/(double)par.nb_simul);
    res.ICsup=res.val+par.quant*sqrt(var/(double)par.nb_simul);
    res.niv=par.niv;
    return res;
}
simulation echange_MC_conditionner(parametre par){
    simulation res;
    vecteur S1,S2;
    vector<vecteur> inde=normal_indep(par.nb_simul,sqrt(par.T),sqrt(par.T));
    vecteur W1=inde[0];
    vecteur vec=psi(W1,par);
    res.val=E(vec);
    double var=V(vec);
    //cout<<"V(X|Y)="<<var<<endl;
    res.ICinf=res.val - par.quant*sqrt(var/(double)par.nb_simul);
    res.ICsup=res.val + par.quant*sqrt(var/(double)par.nb_simul);
    res.niv=par.niv;
    return res;
}

simulation spread_MC(parametre par){
    simulation res;
    vecteur S1,S2;
    vector<vecteur> inde=normal_indep(par.nb_simul,sqrt(par.T),sqrt(par.T));
    vecteur W2=par.rho*inde[0]+sqrt(1-par.rho*par.rho)*inde[1];
    vecteur W1=inde[0];
    S1=par.alpha*par.S1*exp((par.r-par.sig1*par.sig1/2.)*par.T+par.sig1*W1);
    S2=par.beta*par.S2*exp((par.r-par.sig2*par.sig2/2.)*par.T+par.sig2*W2);
    vecteur Splus=Plus(S1-S2-exp(-par.r*par.T)*par.K);
    res.val=exp(-par.r*par.T)*E(Splus);
    double var=exp(-2*par.r*par.T)*V(Splus);
    res.ICinf=res.val-par.quant*sqrt(var/(double)par.nb_simul);
    res.ICsup=res.val+par.quant*sqrt(var/(double)par.nb_simul);
    res.niv=par.niv;
    return res;
}
simulation spread_MC_conditionner(parametre par){
    simulation res;
    vector<vecteur> inde=normal_indep(par.nb_simul,sqrt(par.T),sqrt(par.T));
    vecteur W1=inde[0];
    vecteur vec=psi_spread(W1,par);
    res.val=E(vec);
    double var=V(vec);
    res.ICinf=res.val - par.quant*sqrt(var/(double)par.nb_simul);
    res.ICsup=res.val + par.quant*sqrt(var/(double)par.nb_simul);
    res.niv=par.niv;
    return res;
}
simulation variableControle(parametre par){
    simulation res;
    vector<vecteur> inde=normal_indep(par.nb_simul,sqrt(par.T),sqrt(par.T));
    vecteur W2=inde[1];
    vecteur vec=psi_controle(W2,par);
    res.val=E(vec)+par.alpha*par.S1-par.beta*par.S2 - par.K*exp(-par.r*par.T);
    double var=V(vec);
    res.ICinf=res.val - par.quant*sqrt(var/(double)par.nb_simul);
    res.ICsup=res.val + par.quant*sqrt(var/(double)par.nb_simul);
    res.niv=par.niv;
    return res;
}
double forwardBestof(parametre par){
    double sig=sqrt(par.T*(par.sig1*par.sig1+par.sig2*par.sig2-2*par.rho*par.sig1*par.sig2));
    double P_1=par.alpha*par.S1*I((1./sig)*(log(par.alpha*par.S1/(par.beta*par.S2))+(1./2.)*pow(sig,2)))-par.beta*par.S2*I((1./sig)*(log(par.alpha*par.S1/(par.beta*par.S2))-(1./2.)*pow(sig,2)));
    double P_2=par.beta*par.S2*I((1./sig)*(log(par.beta*par.S2/(par.alpha*par.S1))+(1./2.)*pow(sig,2)))-par.alpha*par.S1*I((1./sig)*(log(par.beta*par.S2/(par.alpha*par.S1))-(1./2.)*pow(sig,2)));
    double res=1./2.*(par.alpha*par.S1+par.beta*par.S2+P_1+P_2) - par.K*exp(-par.r*par.T);
    return res;
}
simulation MC_forwardBestof(parametre par){
    simulation res;
    vecteur S1,S2;
    vector<vecteur> inde=normal_indep(par.nb_simul,sqrt(par.T),sqrt(par.T));
    vecteur W2=par.rho*inde[0]+sqrt(1-par.rho*par.rho)*inde[1];
    vecteur W1=inde[0];
    S1=par.alpha*par.S1*exp((par.r-par.sig1*par.sig1/2.)*par.T+par.sig1*W1);
    S2=par.beta*par.S2*exp((par.r-par.sig2*par.sig2/2.)*par.T+par.sig2*W2);
    vecteur Splus=(max(S1,S2)-par.K);
    res.val=exp(-par.r*par.T)*E(Splus);
    double var=exp(-2*par.r*par.T)*V(Splus);
    res.ICinf=res.val-par.quant*sqrt(var/(double)par.nb_simul);
    res.ICsup=res.val+par.quant*sqrt(var/(double)par.nb_simul);
    res.niv=par.niv;
    return res;
}

simulation spread_Controle_2(parametre par){
    simulation res;
    vecteur S1,S2;
    vector<vecteur> inde=normal_indep(par.nb_simul,sqrt(par.T),sqrt(par.T));
    vecteur W2=par.rho*inde[0]+sqrt(1-par.rho*par.rho)*inde[1];
    vecteur W1=inde[0];
    S1=par.alpha*par.S1*exp((par.r-par.sig1*par.sig1/2.)*par.T+par.sig1*W1);
    S2=par.beta*par.S2*exp((par.r-par.sig2*par.sig2/2.)*par.T+par.sig2*W2);
    vecteur Splus=f(S1-S2,par.K);
    double sig=sqrt(par.T*(par.sig1*par.sig1+par.sig2*par.sig2-2*par.rho*par.sig1*par.sig2));
    double P_a=par.alpha*par.S1*I((1./sig)*(log(par.alpha*par.S1/(par.beta*par.S2))+(1./2.)*pow(sig,2)))-par.beta*par.S2*I((1./sig)*(log(par.alpha*par.S1/(par.beta*par.S2))-(1./2.)*pow(sig,2)));
    res.val=exp(-par.r*par.T)*E(Splus)+P_a;
    double var=exp(-2*par.r*par.T)*V(Splus);
    res.ICinf=res.val-par.quant*sqrt(var/(double)par.nb_simul);
    res.ICsup=res.val+par.quant*sqrt(var/(double)par.nb_simul);
    res.niv=par.niv;
    return res;

}
