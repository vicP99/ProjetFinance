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
    double S1,S2,S1_anti,S2_anti;
    vector<vecteur> inde=normal_indep(nb_simul,sqrt(T),sqrt(T));
    vecteur W2=rho*inde[0]+sqrt(1-rho*rho)*inde[1];
    vecteur W1=inde[0];
    vecteur Splus(nb_simul);
    vecteur Splus_anti(nb_simul);
    for(int i=0;i<nb_simul;i++)
    {
        S1=alpha*exp((r-sig1*sig1/2.)*T+sig1*W1[i]);
        S2=beta*exp((r-sig2*sig2/2.)*T+sig2*W2[i]);
        S1_anti=alpha*exp((r-sig1*sig1/2.)*T+sig1*(-W1[i]));
        S2_anti=beta*exp((r-sig2*sig2/2.)*T+sig2*(-W2[i]));
        if(S1-S2<0)
        {
            Splus[i]=0;
        }
        else
        {
            Splus[i]=S1-S2;
        }
        if(S1_anti-S2_anti<0)
        {
            Splus_anti[i]=0;
        }
        else
        {
            Splus_anti[i]=S1_anti-S2_anti;
        }
        
    }
    double P=exp(-r*T)*E(Splus);
    double P_anti=exp(-r*T)*(E(Splus)+E(Splus_anti))/2.;
    double var=exp(-2*r*T)*V(Splus);
    double var_anti=exp(-2*r*T)/4.*(V(Splus)+V(Splus_anti)+(1/((double) nb_simul))*Cov(Splus,Splus_anti));
    double sig=sqrt(T*(sig1*sig1+sig2*sig2-2*rho*sig1*sig2));
    double P_a=I(sig/2)-I(-sig/2);
    cout<<"P_monte-carlo="<<P<<"\n";
    cout<<"P_analytique="<<P_a<<"\n";
    cout<<"erreur="<<1.96*sqrt(var/(double)nb_simul)<<"\n";
    cout<<"intervalle de confiance ["<<P-1.96*sqrt(var/(double)nb_simul)<<","<<P+1.96*sqrt(var/(double)nb_simul)<<"]\n";
    cout<<"\n";
    cout<<"P_antithetique="<<P_anti<<"\n";
    cout<<"erreur="<<1.96*sqrt(var_anti/(2*(double)nb_simul))<<"\n";
    cout<<"intervalle de confiance ["<<P_anti-1.96*sqrt(var_anti/(2*(double)nb_simul))<<","<<P_anti+1.96*sqrt(var_anti/(2*(double)nb_simul))<<"]\n";
    return 1;
}