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
    //classe de vecteur pour faciliter l'implémentation des méthodes de monte carlo
    public:
    vector<double> v;

   vecteur(int taille=0,double x=0.): v(vector<double>(taille,x)) {}  
   vecteur(const vecteur & u): v(u.v) {}  
   double & operator[](int i){return v[i];}
   double operator[](int i)const {return v[i];}
};

//affichage des vecteurs
ostream & operator <<( ostream & flux,const vecteur & v);

//opérations classique sur les vecteurs (addition,soustraction,multiplication par un scalaire)
vecteur operator +(const vecteur & v1,const vecteur & v2);
vecteur operator +(const double a,const vecteur & v);
vecteur operator +(const vecteur & v,const double a);
vecteur operator -(const vecteur & v1,const vecteur & v2);
vecteur operator -(const double a,const vecteur & v);
vecteur operator -(vecteur & v,const double a);
vecteur operator *(double a,const vecteur & v);
vecteur operator *(const vecteur & v,double a);

//quelques fonction vectorisés utile dans le cadre du projet
vecteur exp(const vecteur&);
vecteur Plus(const vecteur&); // max(vec,0)
vecteur max(const vecteur&, const vecteur&); //maximum element par element
vecteur f(const vecteur& v, double K); //Q15


//simulation de loi uniforme
vecteur loi_unif(double a, double b, int nb_points);

//Méthode de box muller pour simuler 2 loi normale indépendante
vector<vecteur> normal_indep(int nbpoint,double sig1=1,double sig2=1,double mu1=0,double mu2=0);

//espérance & variance empirique (sans biais)
double E(vecteur v);

double V(vecteur v);
double I(double x); //fonction de répartition de la loi normale (0,1)
double Cov(vecteur v1,vecteur v2); //covariance empirique de 2 vecteurs



class parametre{
    //la classe paramètre contient tout les paramètres du problème
    //son rôle est de faciliter l'implémentation des méthodes et la lecture du code
    // tout en permettant de changer facilement les valeurs des paramètres d'une simulation à l'autre
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
    double quant; //quantile au niveau niv
    double niv; //niveau de l'intervalle de confiance (généralement 0.9)
    int nb_simul;
};

class simulation{
    //classe contenant les résultats de simulation
    //valeur obtenue par monte carlo avec la borne inf et sup de l'intervalle de confiance au niveau niv
    public:
    double val;
    double ICinf;
    double ICsup;
    double niv;
};
//pour cette première fonction psi on a psi= E(f(y,W3)) voir Q4
double psi(double y, parametre par);

//fonction psi vectorisée
vecteur psi(const vecteur& v,parametre par);

//psi_spread=E(f_spread(y,W3)) voir Q7
double psi_spread(double y, parametre par);
vecteur psi_spread(const vecteur& v,parametre par);

//méthode de conditionnement combinée avec la méthode de la variable de contrôle    
double psi_controle(double y, parametre par);
vecteur psi_controle(const vecteur& v,parametre par);


//méthodes de monte carlo, voir function.cpp pour le détail
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