#include "../include/header.hpp"

ofstream out("out/sortie.txt");

int main(){
    parametre par;
    par.nb_simul=100000;
    par.T=2;
    par.rho=0.5;
    par.r=0.01;
    par.sig1=0.25;
    par.sig2=0.3;
    par.beta=1;
    par.alpha=12;
    par.S2=1;
    par.S1=1;
    par.K=0.2;
    
    double sig=sqrt(par.T*(par.sig1*par.sig1+par.sig2*par.sig2-2*par.rho*par.sig1*par.sig2));
    double P_a=par.alpha*par.S1*I((1./sig)*(log(par.alpha*par.S1/(par.beta*par.S2))+(1./2.)*pow(sig,2)))-par.beta*par.S2*I((1./sig)*(log(par.alpha*par.S1/(par.beta*par.S2))-(1./2.)*pow(sig,2)));



    cout<<"P="<<P_a<<endl;
    cout<<echange_MC(par);
    cout<<echange_MC_conditionner(par);
    cout<<spread_MC(par);
    cout<<spread_MC_conditionner(par);

     int n=100;
    for (int i=1;i<n;i++){
        par.rho=-1+2./(double)n*i;
        out<<spread_MC_conditionner(par).val<<"\n";
    } 
    return 1;
}