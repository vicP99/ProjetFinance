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

double I(double x)
{
    double t=1/(1+0.2316419*x);
    double res=(1/sqrt(2*M_PI))*exp(-x*x/2)*(0.319381530*t-0.356563782*pow(t,2)+1.781477937*pow(t,3)-1.821255978*pow(t,4)+1.330274429*pow(t,5));
    if(x>0)
    {
        return 1.-res;
    }
    else
    {
        return res;
    }
}