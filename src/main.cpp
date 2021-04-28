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
    par.alpha=1.2;
    par.S2=1;
    par.S1=1;
    par.K=0.2;
    par.quant=1.644854;
    par.niv=90;
    
    double sig=sqrt(par.T*(par.sig1*par.sig1+par.sig2*par.sig2-2*par.rho*par.sig1*par.sig2));
    double P_a=par.alpha*par.S1*I((1./sig)*(log(par.alpha*par.S1/(par.beta*par.S2))+(1./2.)*pow(sig,2)))-par.beta*par.S2*I((1./sig)*(log(par.alpha*par.S1/(par.beta*par.S2))-(1./2.)*pow(sig,2)));



   // cout<<"P="<<P_a<<endl;
    //cout<<echange_MC(par);
  //  cout<<echange_MC_conditionner(par);
    cout<<spread_MC(par);
    cout<<spread_MC_conditionner(par);
    cout<<variableControle(par);
    cout<<spread_Controle_2(par);
    cout<<MC_forwardBestof(par);
    cout<<forwardBestof(par);

//affichage de la fonction d'echange avec rédution de variance pour rho dans ]-1,1[

/* int n=100;
for (int i=1;i<n;i++){
    par.nb_simul= 10000; 
    par.T=2;
    par.r=0.01;
    par.sig1=0.25;
    par.sig2=0.3;
    par.beta=1;
    par.alpha=1;
    par.S2=1;
    par.K=0.2;
    par.rho= -1.+(double)i*2./(double)n; 

    simulation SMC=echange_MC_conditionner(par);
    out<<SMC.val<<"\n";
    out<<SMC.val-SMC.ICinf<<"\n";
    }   */

//Estimateur de spread par conditionnement tracer de trajectoire
/*  int n=100;
    for (int i=1;i<=n;i++){
        par.nb_simul= 1000*i; 
        par.rho=0.5;
        par.T=2;
        par.r=0.01;
        par.sig1=0.25;
        par.sig2=0.3;
        par.beta=1;
        par.alpha=1.2;
        par.S2=1;
        par.K=0.2;
        simulation SMC=spread_MC_conditionner(par);
        out<<SMC.val<<"\n";
        out<<SMC.val-SMC.ICinf<<"\n";
    }  */

//affichage de la fonction spread avec rédution de variance pour rho dans ]-1,1[

/* int n=100;
for (int i=1;i<n;i++){
    par.nb_simul= 10000; 
    par.T=2;
    par.r=0.01;
    par.sig1=0.25;
    par.sig2=0.3;
    par.beta=1;
    par.alpha=1;
    par.S2=1;
    par.K=0.2;
    par.rho= -1.+(double)i*2./(double)n; 

    simulation SMC=spread_MC_conditionner(par);
    out<<SMC.val<<"\n";
    out<<SMC.val-SMC.ICinf<<"\n";
    }   */
//affichage de la fonction spread avec rédution de variance pour K dans ]-1,1[

/* int n=100;
for (int i=0;i<=n;i++){
    par.nb_simul= 10000; 
    par.T=2;
    par.r=0.01;
    par.sig1=0.25;
    par.sig2=0.3;
    par.beta=1;
    par.alpha=1.2;
    par.S2=1;
    par.rho=0.5;
    par.K= -1.+(double)i*2./(double)n; 

    simulation SMC=spread_MC_conditionner(par);
    out<<SMC.val<<"\n";
    out<<SMC.val-SMC.ICinf<<"\n";
    }   */
//affichage de la fonction spread avec rédution de variance controle pour sigma dans [0,0.]^2[

/* int n=100;
for (int i=0;i<=n;i++){
    for(int j=0;j<=n;j++){
        par.nb_simul= 10000; 
        par.T=2;
        par.r=0.01;
        par.sig1=(double)i*0.8/(double)n;
        par.sig2=-1.+(double)j*0.8/(double)n;
        par.beta=1;
        par.alpha=1.2;
        par.S2=1;
        par.rho=0.5;
        par.K=0.2; 

        simulation SMC=variableControle(par);
        out<<SMC.val<<"\n";
        out<<SMC.val-SMC.ICinf<<"\n";
    }

    }   */

//comparaison option échange conditionner et sans conditionnement

/*     int n=100;
    for (int i=1;i<=n;i++){
        par.nb_simul= 1000*i; 
        simulation EMC=echange_MC(par);
        simulation SMC=echange_MC_conditionner(par);
        out<<EMC.val<<"\n";
        out<<SMC.val<<"\n";
        out<<EMC.val-EMC.ICinf <<"\n";
        out<<SMC.val-SMC.ICinf<<"\n";
    }  */

    return 1;
}