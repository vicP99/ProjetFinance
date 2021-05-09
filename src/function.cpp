#include "../include/header.hpp"


using namespace std;



//définition des opérations sur les vecteurs
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
    flux<<"La valeur est: "<<v.val<<" et l'interval de confiance à "<<v.niv<<"\% est: ["<<v.ICinf<<","<<v.ICsup<<"] avec une erreur de "<<v.val-v.ICinf<<endl;
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
//fonction utile pour la question 15
vecteur f(const vecteur& v, double K){
    return Plus(v-K) - Plus(v);
} 

//--------------------------


//Simulation d'une loi uniforme
vecteur loi_unif(double a, double b, int nb_points){
    vecteur res(nb_points);
    for (int i=0;i<nb_points;i++){
        res[i]  = a+(b-a)*((double)(rand()))/((double)RAND_MAX);
        while(res[i]==0){
            res[i]  = a+(b-a)*((double)(rand()))/((double)RAND_MAX); //On veut une loi uniforme sur ]0,1] pour appliquer Box-Muller
        }
    }
    return res;
}

//Simulation de loi normales indépendantes
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

//Calcul de moyenne empirique
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

//Calcul de variance empirique
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

//Calcul de covariance empirique(inutile dans tout le projet)
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

//Fonction de répartition d'une loi normale
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
//Sert pour simuler le prix d'une option d'échange par conditionnement. (Donne le résultat analytique E(...|W1=y))
double psi(double y,parametre par){
    double alphaModifie=par.alpha*par.S1*exp((pow(par.sig1,2))/2*-par.T+par.sig1*y);
    double betaModifie=par.beta*par.S2*exp((pow(par.sig2*par.rho,2))*(-par.T/2)+par.sig2*par.rho*y);
    double sigMod=par.sig2*sqrt((1-pow(par.rho,2))*par.T);
    return alphaModifie*I((log(alphaModifie/betaModifie)/sigMod+sigMod/2)) -betaModifie*I((log(alphaModifie/betaModifie)/sigMod-sigMod/2));
}
//Operateur sur un vecteur
vecteur psi(const vecteur& v,parametre par){
    vecteur res(v.v.size());
    for (uint i=0;i<v.v.size();i++){
        res[i]=psi(v[i],par);
    }
    return res;
}
//Même chose mais pour le conditionnement de l'option spread
double psi_spread(double y,parametre par){
    double alphaModifie=par.alpha*par.S1*exp((pow(par.sig1,2))/2*-par.T+par.sig1*y)- par.K*exp(-par.r*par.T);
    if(alphaModifie <= 0){
        return 0;
    }
    double betaModifie=par.beta*par.S2*exp((pow(par.sig2*par.rho,2))*(-par.T/2)+par.sig2*par.rho*y);
    double sigMod=par.sig2*sqrt((1-pow(par.rho,2))*par.T);
    return alphaModifie*I((log(alphaModifie/betaModifie)/sigMod+sigMod/2)) -betaModifie*I((log(alphaModifie/betaModifie)/sigMod-sigMod/2));
}
//Pour simuler par conditionnement la variable de controle
double psi_controle(double y,parametre par){
    double betaModifie=par.alpha*par.S1*exp((pow(par.sig1,2))/2*-par.T+par.sig1*y);
    double alphaModifie=par.beta*par.S2*exp((pow(par.sig2*par.rho,2))*(-par.T/2)+par.sig2*par.rho*y)+par.K*exp(-par.r*par.T);
    if(alphaModifie <=0){
        return 0;
    }
    double sigMod=par.sig1*sqrt((1-pow(par.rho,2))*par.T);
    return alphaModifie*I((log(alphaModifie/betaModifie)/sigMod+sigMod/2)) -betaModifie*I((log(alphaModifie/betaModifie)/sigMod-sigMod/2));
}
//operateur sur les vecteurs
vecteur psi_spread(const vecteur& v,parametre par){
    vecteur res(v.v.size());
    for (uint i=0;i<v.v.size();i++){
        res[i]=psi_spread(v[i],par);
    }
    return res;
}
//operateur sur les vecteurs
vecteur psi_controle(const vecteur& v,parametre par){
    vecteur res(v.v.size());
    for (uint i=0;i<v.v.size();i++){
        res[i]=psi_controle(v[i],par);
    }
    return res;
}
//Prix de l'option d'échange par Monte Carlo classique
simulation echange_MC(parametre par){
//simulation de S1 S2
    simulation res;
    vecteur S1,S2;
    vector<vecteur> inde=normal_indep(par.nb_simul,sqrt(par.T),sqrt(par.T));
    vecteur W2=par.rho*inde[0]+sqrt(1-par.rho*par.rho)*inde[1];
    vecteur W1=inde[0];
    S1=par.alpha*par.S1*exp((par.r-par.sig1*par.sig1/2.)*par.T+par.sig1*W1);
    S2=par.beta*par.S2*exp((par.r-par.sig2*par.sig2/2.)*par.T+par.sig2*W2);
//Simulation de (alpha*S1 -beta*S2)+
    vecteur Splus=Plus(S1-S2);
//resulats de la simulation par monte carlo
    res.val=exp(-par.r*par.T)*E(Splus);
    double var=exp(-2*par.r*par.T)*V(Splus);
    //cout<<"V(X)="<<var<<endl;
    res.ICinf=res.val-par.quant*sqrt(var/(double)par.nb_simul);
    res.ICsup=res.val+par.quant*sqrt(var/(double)par.nb_simul);
    res.niv=par.niv;
    return res;
}
//Prix de l'option d'échange par conditionnement
simulation echange_MC_conditionner(parametre par){
    simulation res;
    vector<vecteur> inde=normal_indep(par.nb_simul,sqrt(par.T),sqrt(par.T));
    vecteur W1=inde[0];
    vecteur vec=psi(W1,par);//E(X|Y=y) pour n_simulation de Y

    //Resultat de la simulation
    res.val=E(vec);
    double var=V(vec);
    //cout<<"V(X|Y)="<<var<<endl;
    res.ICinf=res.val - par.quant*sqrt(var/(double)par.nb_simul);
    res.ICsup=res.val + par.quant*sqrt(var/(double)par.nb_simul);
    res.niv=par.niv;
    return res;
}
//Prix de l'option spread par monte carlo classique (pareil que l'option d'échange)
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
//spread par conditionnement
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
//spread par cariable de controle + conditionnement
simulation variableControle(parametre par){
    srand(0);
    simulation res;
    vector<vecteur> inde=normal_indep(par.nb_simul,sqrt(par.T),sqrt(par.T));
    vecteur W2=inde[1];
    vecteur vec=psi_controle(W2,par); //conditionnement sur la variable de contrôle
    res.val=E(vec)+par.alpha*par.S1-par.beta*par.S2 - par.K*exp(-par.r*par.T); //ajout de la partie analytique
    double var=V(vec);
    res.ICinf=res.val - par.quant*sqrt(var/(double)par.nb_simul);
    res.ICsup=res.val + par.quant*sqrt(var/(double)par.nb_simul);
    res.niv=par.niv;
    return res;
}
//Calcul analytique du BestOf
double forwardBestof(parametre par){
    double sig=sqrt(par.T*(par.sig1*par.sig1+par.sig2*par.sig2-2*par.rho*par.sig1*par.sig2));
    double P_1=par.alpha*par.S1*I((1./sig)*(log(par.alpha*par.S1/(par.beta*par.S2))+(1./2.)*pow(sig,2)))-par.beta*par.S2*I((1./sig)*(log(par.alpha*par.S1/(par.beta*par.S2))-(1./2.)*pow(sig,2)));
    double P_2=par.beta*par.S2*I((1./sig)*(log(par.beta*par.S2/(par.alpha*par.S1))+(1./2.)*pow(sig,2)))-par.alpha*par.S1*I((1./sig)*(log(par.beta*par.S2/(par.alpha*par.S1))-(1./2.)*pow(sig,2)));
    double res=1./2.*(par.alpha*par.S1+par.beta*par.S2+P_1+P_2) - par.K*exp(-par.r*par.T);
    return res;
}
//SImulation du bestOf par Monte Carlo classique
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
//méthode avec une autre variable de contrôle sur Spread
simulation spread_Controle_2(parametre par){
    simulation res;
    vecteur S1,S2;
    vector<vecteur> inde=normal_indep(par.nb_simul,sqrt(par.T),sqrt(par.T));
    vecteur W2=par.rho*inde[0]+sqrt(1-par.rho*par.rho)*inde[1];
    vecteur W1=inde[0];
    S1=par.alpha*par.S1*exp((par.r-par.sig1*par.sig1/2.)*par.T+par.sig1*W1);
    S2=par.beta*par.S2*exp((par.r-par.sig2*par.sig2/2.)*par.T+par.sig2*W2);
    vecteur Splus=f(S1-S2,par.K); //variable de controle
    double sig=sqrt(par.T*(par.sig1*par.sig1+par.sig2*par.sig2-2*par.rho*par.sig1*par.sig2));
    double P_a=par.alpha*par.S1*I((1./sig)*(log(par.alpha*par.S1/(par.beta*par.S2))+(1./2.)*pow(sig,2)))-par.beta*par.S2*I((1./sig)*(log(par.alpha*par.S1/(par.beta*par.S2))-(1./2.)*pow(sig,2))); //prix analytique de l'option d'échange
    res.val=exp(-par.r*par.T)*E(Splus)+P_a;
    double var=exp(-2*par.r*par.T)*V(Splus);
    res.ICinf=res.val-par.quant*sqrt(var/(double)par.nb_simul);
    res.ICsup=res.val+par.quant*sqrt(var/(double)par.nb_simul);
    res.niv=par.niv;
    return res;

}
