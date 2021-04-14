#include "../include/header.hpp"

ofstream out("out/sortie.txt");

int main(){
    parametre par;
    par.nb_simul=10000;
    par.T=2;
    par.rho=0.5;
    par.r=0.01;
    par.sig1=0.25;
    par.sig2=0.3;
    par.beta=1;
    par.alpha=1;
    par.S2=1;
    par.S1=1;
    double S1,S2,S1_anti,S2_anti;
    vector<vecteur> inde=normal_indep(par.nb_simul,sqrt(par.T),sqrt(par.T));
    vecteur W2=par.rho*inde[0]+sqrt(1-par.rho*par.rho)*inde[1];
    vecteur W1=inde[0];
    vecteur Splus(par.nb_simul);
    vecteur Splus_anti(par.nb_simul);
    for(int i=0;i<par.nb_simul;i++)
    {
        S1=par.alpha*exp((par.r-par.sig1*par.sig1/2.)*par.T+par.sig1*W1[i]);
        S2=par.beta*exp((par.r-par.sig2*par.sig2/2.)*par.T+par.sig2*W2[i]);
        S1_anti=par.alpha*exp((par.r-par.sig1*par.sig1/2.)*par.T+par.sig1*(-W1[i]));
        S2_anti=par.beta*exp((par.r-par.sig2*par.sig2/2.)*par.T+par.sig2*(-W2[i]));
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
    double P=exp(-par.r*par.T)*E(Splus);
    double P_anti=exp(-par.r*par.T)*(E(Splus)+E(Splus_anti))/2.;
    double var=exp(-2*par.r*par.T)*V(Splus);
    double var_anti=exp(-2*par.r*par.T)/4.*(V(Splus)+V(Splus_anti)+(1/((double) par.nb_simul))*Cov(Splus,Splus_anti));
    double sig=sqrt(par.T*(par.sig1*par.sig1+par.sig2*par.sig2-2*par.rho*par.sig1*par.sig2));
    double P_a=I(sig/2)-I(-sig/2);


    double MCcond=0;
    for (int i=0;i<par.nb_simul;i++){
        MCcond+=psi(inde[1][i],par);
    }
    y=0;
    double K=par.alpha*par.S1*exp((-pow(par.sig1,2))/2*par.T+par.sig1*y);

    double MCcond=0;
    for (int i=0;i<par.nb_simul;i++){
        ST=par.beta*par.S2*exp(()*par.T + par.sg2*par.rho*W1[i])
    }





    cout<<"P_monte-carlo="<<P<<"\n";
    cout<<"var="<<var<<"\n";
    cout<<"P_analytique="<<P_a<<"\n";
    cout<<"P_condtionnement="<<MCcond<<"\n";
    cout<<"erreur="<<1.96*sqrt(var/(double)par.nb_simul)<<"\n";
    cout<<"intervalle de confiance ["<<P-1.96*sqrt(var/(double)par.nb_simul)<<","<<P+1.96*sqrt(var/(double)par.nb_simul)<<"]\n";
    cout<<"\n";
    cout<<"P_antithetique="<<P_anti<<"\n";
    cout<<"erreur="<<1.96*sqrt(var_anti/(2*(double)par.nb_simul))<<"\n";
    cout<<"intervalle de confiance ["<<P_anti-1.96*sqrt(var_anti/(2*(double)par.nb_simul))<<","<<P_anti+1.96*sqrt(var_anti/(2*(double)par.nb_simul))<<"]\n";
    return 1;
}