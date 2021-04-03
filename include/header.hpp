#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <vector>
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
vecteur operator *(double a,const vecteur & v);
vecteur operator *(const vecteur & v,double a);


vecteur loi_unif(double a, double b, int nb_points);
vector<vecteur> normal_indep(int nbpoint,double sig1=1,double sig2=1,double mu1=0,double mu2=0);


double E(vecteur v);

double I(double x);
#endif