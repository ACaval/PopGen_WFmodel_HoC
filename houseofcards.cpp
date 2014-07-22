//_____________________________________LIBRERIE______________________________________//
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>
#include <sstream>
#include <time.h>
#include <iomanip>
using namespace std;

#include "bnldev_mod.c"
#include "poidev_mod.c"
#include "gammln_mod.c"
#include "gammln.c"
#include "gammp.c"
#include "mt19937-64.c"

//__________________________________________________________________________________________//
//______________________________________*** MAIN ***________________________________________//
//__________________________________________________________________________________________//

int main() {

  init_genrand64 ( time(NULL) );
  ofstream wmeantotFile;   wmeantotFile.open("hoc_wmeantot.dat", ios::trunc);
  ofstream kmeantotFile;   kmeantotFile.open("hoc_kmeantot.dat", ios::trunc); 
  ofstream advtotFile;     advtotFile.open("hoc_advtot.dat", ios::trunc);
  ofstream advtot2File;    advtot2File.open("hoc_advtot2.dat", ios::trunc); 
  ofstream vktotFile;      vktotFile.open("hoc_vktot.dat", ios::trunc);
  ofstream vstotFile;      vstotFile.open("hoc_vstot.dat", ios::trunc);
  ofstream histoFile;      histoFile.open("hoc_histo.dat", ios::trunc);
  ofstream histowFile;     histowFile.open("hoc_histow.dat", ios::trunc);
  ofstream numclassFile;   numclassFile.open("hoc_numclass.dat", ios::trunc);
  ofstream benrateFile;    benrateFile.open("hoc_ub.dat", ios::trunc);
 
  cout.setf(ios::fixed,ios::floatfield);
  cout.precision(16); 


 //____Dichiarazione delle variabili globali (iterazioni ensemble)____// 
  int VarPar = 1;
  int iter = 20;
  double N0 = 3.3E+7, N = N0;
  double U0 = 3.91E-04;
  int t0 = 0;
  int tmax = 50001;
  int cell = 20;
  int step = floor(tmax/cell); //floor(tmax / cell);
  int cc = 0;
  int counter0 = 2*step; int counter = counter0;
 
  cout << "t= " << tmax << "; " << "step-width= " << step << "; " << "#points= " << cell << "; " << flush << endl;

  //double w0 = 0.28;   //EXPONENTIAL
  double pot = 5.70;  // POWERLAW //pot > 2 altrimenti l'integrale non converge
  //double beta = 10.8543; double nu = 139.437; double w0 = 0.707293;   // POWERLAW+EXP_TAIL

  double w_in = 1.;  
 
  //____Dichiarazione delle variabili di popolazione (iterazioni temporali)____//
  int k = 1;
  int kmax = 2*tmax;
  int MaxMut=0;
  double wmean = 0.;
  double wmean2 = 0.;
  double kmean = 0.;
  double kmean2 = 0.;
  double benrate = 0.;
  double adv = 0.;
  int numclass = 0;

   //____Dichiarazione delle variabili di classe (iterazioni sulle frequenze)____//
  double nr, Nres = N;
  double pr=0., qr=0., sumpr=0.;
  
  double *Vwmeantot = new double[tmax];
  double *Vwmean2tot = new double [tmax];
  double *Vkmeantot = new double[tmax+1];
  double *Vkmean2tot = new double[tmax];
  double *Vwmeanvar = new double[tmax];
  double *Vkmeanvar = new double[tmax];

  double *Vadvtot = new double[tmax];
  double *Vvel_k = new double[cell]; 
  double *Vvel_s = new double[cell]; 
  double *Vnumclasstot = new double[tmax];
  double *Vbenratetot = new double[tmax];


  //int *Vmutation = new int[kmax]; //semplice vettore del # di mutazioni {0,1,2,3...}
  //double *Vprob = new double[kmax]; //vettore delle pr probabilita'  delle classi r {p0,p1,p2,p3...}
  //double *Vqrob = new double[kmax]; //vettore delle qr prob ridotte delle classi r {q0,q1,q2,q3...}
  //double *Vfreq = new double[kmax]; // vettore delle fr frequenze delle classi r {f0,f1,f2,f3...} 
  double *Vfitness = new double[kmax]; //vettore delle wr fitness delle classi r {w0,w1,w2,w3...}
  double *Vadvantage = new double[kmax]; //vettore delle s(k) delle classi r {s0(k),s1(k),s2(k)...}
  double *Vprobmut = new double [kmax]; //vettore delle U(k)

  //double *Vprob_new = new double[kmax]; //vettore delle pr delle classi r {p0,p1,p2,p3...}
  //double *Vqrob_new = new double[kmax]; //vettore delle qr delle classi r {q0,q1,q2,q3...}
  //double *Vfreq_new = new double[kmax]; //vettore delle fr delle classi r {f0,f1,f2,f3...}

 
  //_____________________________ ||||| FOR SUI PARAMETRI ||||| _____________________________//
  for (int m = 1; m<=VarPar; m++) {

    cout << "Variazione di parametro #" << m << "/" << VarPar << endl;

    // Vfitness[0]=1.;
   
    //___________Inizializzazione della fitness_____________//
    for (int i=0; i<kmax; i++) {

      if(i==0) {Vfitness[i]= w_in;}
      else {
	//Vfitness[i] = Vfitness[i-1] + exp(-Vfitness[i-1]/w0)*(w0+Vfitness[i-1]);  //EXP
	Vfitness[i] = Vfitness[i-1] + (pot*Vfitness[i-1] + 1)/((pow((Vfitness[i-1]+1), pot)*(pot-1))); //POWLAW
	//Vfitness[i] = Vfitness[i-1] + w0*(1.-gammp((nu+2)/beta,pow((Vfitness[i-1]/w0),beta)))*(exp(gammln_mod((nu+2)/beta)))/(exp(gammln_mod((nu+1)/beta))); //POWLAW EXP TAIL
                
            };
 
     Vadvantage[i] = log(Vfitness[i]);                                  

     //Vprobmut[i] = U0 * exp(-Vfitness[i]/w0);          //Distribuzione esponenziale g(w) = 1/w0 * exp(-w/w0)
     Vprobmut[i] = U0 / pow((Vfitness[i]+1),pot);      //Distribuzione power law
     //Vprobmut[i] = U0*(1.-gammp((nu+1)/beta,pow((Vfitness[i]/w0),beta)))/(exp(gammln_mod((nu+1)/beta))); ///POWLAW

     //cout << "k= " << i << "; " << "s= " << Vadvantage[i] << "; " << "w= " << Vfitness[i] << "; " << "U= " << Vprobmut[i] << endl;

   };

     
       //___________________________ **** FOR SUGLI ENSEMBLE ****_____________________________//
       for (int a=1; a<=iter; a++) {
	 cerr << "Iterazione #" << a << "/" << iter << "\r" << flush;

	 //___________Inizializzazione dei vettori_____________//
	 
  vector<int> Vmutation (kmax,0); //semplice vettore del # di mutazioni {0,1,2,3...}	
     for (int i=0; i<kmax; i++) {
	   Vmutation[i]=i;
	   //Vprob[i]=0.;
	   //Vqrob[i]=0.;
	   //Vfreq[i]=0.;
	   //Vprob_new[i]=0.;
	   //Vqrob_new[i]=0.;
	   //Vfreq_new[i]=0.;
	 };

  vector <double> Vprob (kmax,0); //vettore delle pr probabilita'  delle classi r {p0,p1,p2,p3...}
  vector <double> Vqrob (kmax,0); //vettore delle qr prob ridotte delle classi r {q0,q1,q2,q3...}
  vector <double> Vfreq (kmax,0); // vettore delle fr frequenze delle classi r {f0,f1,f2,f3...}
  vector <double> Vprob_new (kmax,0); //vettore delle fr delle classi r {f0,f1,f2,f3...}
  vector <double> Vqrob_new (kmax,0); //vettore delle qr delle classi r {q0,q1,q2,q3...}
  vector <double> Vfreq_new (kmax,0); //vettore delle fr delle classi r {f0,f1,f2,f3...}


	  Vfreq[0] = 1.;
      k = 1;
      counter = counter0; 
      //	int kminimo=0;
      //	int kmin=kminimo;
     //___________________________ °°° FOR TEMPORALE °°° ____________________________//
      for(int t=t0; t<tmax; t++) {
	    k++;
        //cerr << " qui 1" << endl;

	    wmean = 0.;
	    wmean2 = 0.;
	    kmean = 0.;
	    kmean2 = 0.;
	    benrate = 0.;
	    adv = 0.;
	    for (int i=0; i<k; i++) { 
            if(Vfitness[i]<5000){
	            wmean += (Vfreq[i]*Vfitness[i]);
	            wmean2 += (Vfreq[i]*Vfitness[i]*Vfitness[i]);
                };
            kmean += (Vfreq[i]*Vmutation[i]);
            kmean2 += (Vfreq[i]*Vmutation[i]*Vmutation[i]);
            adv += (Vfreq[i]*Vadvantage[i]);
            benrate += (Vfreq[i]*Vprobmut[i]);	
         };
	
	    //benrate = U0 * exp(-wmean/w0);
	    benrate  = U0 / pow((wmean+1),pot);
	    //benrate = U0* (1.-(gammp((nu+1)/beta,pow((wmean/w0),beta))))/(exp(gammln_mod((nu+1)/beta)));


	//_______________ ^^ FOR SULLE CLASSI DI FREQUENZA ^^ _________________//
	    for (int r=0; r<k; r++) {
	  	   //cerr << " qui 2" << endl;

	        if(Vfitness[r]>10) {break;};

	        if (Vfreq[r] != 0.) {numclass++;};

	        if (r!=0) {
	            pr  = ((Vfreq[r-1] * Vprobmut[r] * Vfitness[r-1]) + (Vfreq[r] * (1-Vprobmut[r]) * Vfitness[r]))/wmean;
	        }
	        else{ 
                pr = (Vfreq[r] * (1-Vprobmut[r]) * Vfitness[r])/wmean;};
	            sumpr += pr;
	            if (pr!=0.)
	                { qr= pr/sumpr; 
	                MaxMut=r+1;
	                //  if(kmin==0){ kmin=r;};
	                }
	            else {
	                qr = 0.;
	            };
	  //if (pr==0.) continue;
	  
	        Vprob_new[r] = pr;
	        Vqrob_new[r] = qr;
	        Vprob[r] = Vprob_new[r];
	        Vqrob[r] = Vqrob_new[r];
	
	        if (pr!=0.){
	            // cout << "t= " << t << "; k=" << Vmutation[r] << "; fr=" << Vfreq[r] <<  "; pr= " << Vprob_new[r]  << "; sumpr=" << sumpr << "; qr=" << Vqrob_new[r]  << "; wr=" << Vfitness[r] << endl;
	        };
	  
	        if (r>0   &&   Vprob[r-1] != 0   &&   Vprob[r] == 0){ r=k+1;};	 
	        };

	    qr = 0.;
	  	//cerr << " qui 3" << endl;

	  //MaxMut = k;
	  //	  kminimo=kmin;
	  //cout << "sono kmin " << kmin << " sono kminimo " << kminimo << endl;

	  //___________ + FOR SULLA BINOMIALE + ____________//
	    for (int i = MaxMut; i>=0; i--) {
	        qr = Vqrob_new[i];
	        if((qr*Nres) < 3){	
	            if((qr*Nres)< 1){
		            if((qr*Nres)< 1E-4){ nr=0;
		            }else{ 
		            if((qr*Nres)< 1E-2)
		            {
		                nr= floor(double(poidev_mod(Nres*qr*1000))/1000.);
		            }else{
		                nr= floor(double(poidev_mod(Nres*qr*10))/10.);
		            };
		            };
	            }else{
		            nr= poidev_mod(Nres*qr);
	            };
	        }else{
	            nr= bnldev_mod(qr, Nres);
	        };
	    /* if((qr*Nres)< 1E-2)
	       nr=0.; 
	       else
	       nr = bnldev_mod(qr, Nres);*/
	  
    	    double freq =nr/N;	  
	        Vfreq_new[i] = freq;
	        Nres -= nr;
	        nr = -1.;
	        Vfreq[i] = Vfreq_new[i];

	        if (t == counter) {
	      	   //cerr << " qui 4" << endl;

	            histoFile << i << " " << Vfreq[i] << endl;
	            histowFile << Vfitness[i] << " " << Vfreq[i] << endl;
	        };

	  }; //__ + Chiusura del for sulla binomiale + __ 



	  //Vfreq[r] = pr; //Scommentare questa e commentare il for sulla binomiale per l'algoritmo deterministico
	  pr = 0.;
	  qr = 0.;
	  Nres = N;


	  //	}; // __ ^^ Chiusura del for sulle classi di frequenza ^^ __
 
	  if (t == counter) {
	    // cout << t << endl;
	    histoFile << endl;
	    histowFile << endl;
	  // counter += counter0;
	  };

	  if (t == counter){
          histoFile << endl;
	      counter += counter0;
	    };

	   Vwmeantot[t] += wmean;
	   Vwmean2tot[t] += wmean2;
	   Vkmeantot[t] += kmean;
	   Vkmean2tot[t] += kmean2;
	   Vadvtot[t] += adv;
	   Vbenratetot[t] += benrate;
	   Vnumclasstot[t] += /*(double)*/ numclass;
	   sumpr = 0.;
	   numclass = 0; 


	//			if (t % 500 == 0) {
	//	  cout << t << " " << kmean << endl;
	//	}
	//cout << endl;

      }; // __ °°° Chiusura del for temporale °°° __
    
      //cerr << " qui 5 "  << endl;

    };  // __ **** Chiusura della singola iterazione **** __

    //int conto = 0; 
   cc=0;
    for (int n=0; n<tmax; n++) {
      
      Vwmeantot[n] /= iter;
      Vwmean2tot[n] /= iter;
      Vkmeantot[n] /= iter;
      Vkmean2tot[n] /= iter;
      Vadvtot[n] /= iter;
      Vnumclasstot[n] /= iter;
      Vwmeanvar[n] = sqrt(abs(Vwmean2tot[n]-(Vwmeantot[n]*Vwmeantot[n])));
      Vkmeanvar[n] = sqrt(abs(Vkmean2tot[n]-(Vkmeantot[n]*Vkmeantot[n])));
       
      if ((n % 500==0 && n<=5000) || (n % 1000 == 0 && n<=20000) || (n % 2000 == 0))  {
	 //cerr << n << " " << Vwmeantot[n] << " " << Vwmeanvar[n] << endl;
	 wmeantotFile << n << " " << Vwmeantot[n] << " " << Vwmeanvar[n] << endl;
	 kmeantotFile << n << " " << Vkmeantot[n] << " " << Vkmeanvar[n] << endl;
	 advtotFile << n << " " << Vadvtot[n] << endl;
	 advtot2File << Vkmeantot[n] << " " << Vadvtot[n] << endl;
	 numclassFile << n << " " << Vnumclasstot[n] << endl;
	 benrateFile << n << " " << Vbenratetot[n] << endl;
    
	 if (n>0) {
        //cerr << "step " << step << endl;
        //cerr << "n " << n << endl;
        //cerr << (n-step) <<endl;
	    //Vvel_k[cc] = (Vkmeantot[n] - Vkmeantot[n-step])/step;
	    
        //vktotFile << n << " " << Vvel_k[cc] << endl;
        vktotFile << n << " " << ((Vkmeantot[n] - Vkmeantot[n-step])/double(step)) << endl;
	    //Vvel_s[cc] = (Vadvtot[n] - Vadvtot[n-step])/double(step);
	    //vstotFile << n << " " << Vvel_s[cc] << endl;
        vstotFile << n << " " << (Vadvtot[n] - Vadvtot[n-step])/double(step) << endl;
        
	    if (cc == cell){
	        vktotFile << endl;
	        vstotFile << endl;
	    };
        
	   };
	 //cout << "qui no 2" << endl;
	 cc++;
     };

      if (n == (tmax-1)) {
	wmeantotFile << endl;
	kmeantotFile << endl;
	advtotFile << endl;
	numclassFile << endl;
      }
    }

  
  vktotFile << endl;
  vstotFile << endl;


  U0 += 1E-04;

  }; 

  // __ ||||| Chiusura del for sui parametri ||||| __
 
 
  wmeantotFile.close();
  kmeantotFile.close();
  histoFile.close();
  histowFile.close();
  advtotFile.close();
  advtot2File.close();
  vktotFile.close();
  vstotFile.close();
  numclassFile.close();
  benrateFile.close();
 
  
  delete Vwmeantot;
  delete Vwmean2tot;
  delete Vkmeantot;
  delete Vkmean2tot;
  delete Vkmeanvar;
  
  delete Vadvtot;
  delete Vvel_k;
  delete Vvel_s;
  delete Vnumclasstot;
  delete Vbenratetot;

  /*
  delete Vmutation;
  delete Vprob;
  delete Vqrob;
  delete Vfreq;
    */
  delete Vfitness;
  delete Vadvantage;
  delete Vprobmut;
  //delete Vprob_new;
  //delete Vqrob_new;
  //delete Vfreq_new;
 
 
  return 0;
  
}//Chiusura del main

