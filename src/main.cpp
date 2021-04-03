#include "../include/header.hpp"

ofstream out("out/sortie.txt");

int main(){
    int nb_simul=10000;
    double T=2;
    double rho=0.5;
    double r=0.01;
    double sig1=0.25;
    double sig2=0.3;
    double beta=1,alpha=1;
    double S1,S2;
    vector<vecteur> inde=normal_indep(nb_simul,sqrt(T),sqrt(T));
    vecteur W2=rho*inde[0]+sqrt(1-rho*rho)*inde[1];
    vecteur W1=inde[0];
    vecteur Splus(nb_simul);
    for(int i=0;i<nb_simul;i++)
    {
        S1=alpha*exp((r-sig1*sig1/2.)*T+sig1*W1[i]);
        S2=beta*exp((r-sig2*sig2/2.)*T+sig2*W2[i]);
        if(S1-S2<0)
        {
            Splus[i]=0;
        }
        else
        {
            Splus[i]=S1-S2;
        }
        
    }
    double P=exp(-r*T)*E(Splus);
    cout<<"P="<<P<<"\n";

    return 0;
}