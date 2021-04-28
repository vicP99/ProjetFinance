#include <iostream>
#include <fstream>
#include <stdlib.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <vector>
#include <time.h>
using namespace std;

#ifndef MAIN_H
#define MAIN_H

class vecteur{
    public:
    vector<double> v;

   vecteur(int taille=0,double x=0.): v(vector<double>(taille,x)) {}  
   vecteur(const vecteur & u): v(u.v) {}  
   double & operator[](int i){return v[i];}
   double operator[](int i)const {return v[i];}
};
ostream & operator <<( ostream & flux,const vecteur & v);
vecteur operator +(const vecteur & v1,const vecteur & v2);
vecteur operator +(const double a,const vecteur & v);
vecteur operator +(const vecteur & v,const double a);
vecteur operator -(const vecteur & v1,const vecteur & v2);
vecteur operator -(const double a,const vecteur & v);
vecteur operator -(vecteur & v,const double a);
vecteur operator *(double a,const vecteur & v);
vecteur operator *(const vecteur & v,double a);
vecteur exp(const vecteur&);
vecteur Plus(const vecteur&);
vecteur max(const vecteur&, const vecteur&);
vecteur f(const vecteur& v, double K);



vecteur loi_unif(double a, double b, int nb_points);
vector<vecteur> normal_indep(int nbpoint,double sig1=1,double sig2=1,double mu1=0,double mu2=0);


double E(vecteur v);

double V(vecteur v);
double I(double x);
double Cov(vecteur v1,vecteur v2);



class parametre{
    public:
    double alpha;
    double beta;
    double S1;
    double S2;
    double T;
    double sig1;
    double sig2;
    double rho;
    double r;
    double K;
    double quant;
    double niv;
    int nb_simul;
};

class simulation{
    public:
    double val;
    double ICinf;
    double ICsup;
    double niv;
};
double psi(double y, parametre par);
vecteur psi(const vecteur& v,parametre par);
double psi_spread(double y, parametre par);
vecteur psi_spread(const vecteur& v,parametre par);
double psi_controle(double y, parametre par);
vecteur psi_controle(const vecteur& v,parametre par);

simulation echange_MC(parametre par);
simulation echange_MC_conditionner(parametre par);
simulation spread_MC(parametre par);
simulation spread_MC_conditionner(parametre par);
simulation variableControle(parametre par);
ostream & operator <<( ostream & flux,const simulation & v);
double forwardBestof(parametre par);
simulation MC_forwardBestof(parametre par);
simulation spread_Controle_2(parametre par);
#endif