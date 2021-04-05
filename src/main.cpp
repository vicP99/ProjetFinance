#include "../include/header.hpp"

ofstream out("out/sortie.txt");

int main(){
    //srand (time(NULL));
    int nb_simul=100;
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
    double var=exp(-2*r*T)*V(Splus);
    double sig=sqrt(T*(sig1*sig1+sig2*sig2-2*rho*sig1*sig2));
    double P_a=I(1/sig+sig/2)-I(1/sig-sig/2);
    cout<<"P="<<P<<"\n";
    cout<<"P_analytique="<<P_a<<"\n";
    cout<<"intervalle de confiance ["<<P-1.645*sqrt(var)/(double)nb_simul<<","<<P+1.645*sqrt(var)/(double)nb_simul<<"]\n";

    return 1;
}