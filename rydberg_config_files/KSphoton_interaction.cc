#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "RydbergFerenc.h"
#include "KConst.h"
#include "KRandom.h"
#include "KSInteractionsMessage.h"
#include "QuadGaussLegendre.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include "random.c"
#include "Clebsch-Gordon_look_up_table.cc"

using katrin::KRandom;
using namespace Kassiopeia {


double PLS(double f,double FWHM, double I,int n,int l,int np, int sign);// LS=Light Source
float h_cross = 6.5821e-16;  //eV/s

////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////

double RadInt2BoundFreeBurgess(int n, int l, double E, int sign)
{

    int lp = l + sign;
    if (lp < 0)
        return 0.;

    if (n < 1 || l < 0 || l > n - 1) {
        puts("Message from function RadInt2BoundFreeBurgess:");
        puts("n has to be larger than 0!");
        puts("l has to be larger than -1 and smaller than n !");
        puts("Program running is stopped !!! ");
        exit(1);
    }

    if (n > nMAX) {
        puts("Message from function RadInt2BoundFreeBurgess");
        puts("n cannot be larger than nMAX!");
        puts("Program running is stopped !!! ");
        exit(1);
    }

    if (E < 0.) {
        puts("Message from function RadInt2BoundFreeComplex");
        puts("The outgoing electron energy E has to be positive!");
        puts("Program running is stopped !!! ");
        exit(1);
    }

    if (sign != -1 && sign != +1) {
        puts("Message from function RadInt2BoundFreeBurgess:");
        puts("sign has to be either -1 or +1!");
        puts("Program running is stopped !!! ");
        exit(1);
    }

    double k = sqrt(2. * fabs(E));  // free electron momentum (velocity) in atomic units

    if (k < 1.e-10)  // k=0 would cause division by zero error
    {
        k = 1.e-10;
    }
    double k2 = k * k;
    double D = 1. + n * n * k2;

    // Eq. 30 of Burgess 1965:
    double logg0 = 0.5 * (logpi - LogN[2] - LogNFactor[2 * n - 1]) + LogN[4] + n * (LogN[4] + LogN[n]) - 2 * n;
    //  Eq. 31 of Burgess 1965:
    double logPs = 0.;
    for (int s = 1; s <= n; s++) {
        logPs += log(1. + s * s * k2);
    }
    double logg1 =
        0.5 * (logPs - log(1. - exp(-2. * M_PI / k))) + 2 * n - 2. / k * atan(n * k) - (n + 2) * log(D) + logg0;

    double g1 = 1.;  // we use here g1=1 instead of exp(logg1), to avoid overflow or underflow
    //  Eq. 32 of Burgess 1965:
    double g2 = 0.5 * sqrt((2. * n - 1.) * D) * g1;
    //  Eq. 33 of Burgess 1965:
    double g3 = 1. / (2 * n) * sqrt(D / (1. + (n - 1.) * (n - 1.) * k2)) * g1;
    //  Eq. 34 of Burgess 1965:
    double g4 = (4. + (n - 1.) * D) / (2. * n) * sqrt((2. * n - 1.) / (1. + (n - 2.) * (n - 2.) * k2)) * g3;

    //    printf("logg0,logg1= %12.3f %12.3f  \t\n",logg0,logg1);

    double Exponent = logg1;

    double g = 0., RadInt;

    // We have to multiply the output by 2/M_PI, due to the different wave func. normalizations;
    // Burgess uses a normalization of pi*delta(k^2-k'^2), and I use delta(E-E');
    // therefore my free wave function is sqrt(2./M_PI) times larger than the free wave func. of Burgess.

    double Cfac = 2. / M_PI;

    if (l == n - 1 && sign == 1) {
        g = g1 * exp(Exponent);
        RadInt = n * n * g;
        return RadInt * RadInt * Cfac;
    }
    else if (l == n - 2 && sign == 1) {
        g = g2 * exp(Exponent);
        RadInt = n * n * g;
        return RadInt * RadInt * Cfac;
    }
    else if (l == n - 1 && sign == -1) {
        g = g3 * exp(Exponent);
        RadInt = n * n * g;
        return RadInt * RadInt * Cfac;
    }
    else if (l == n - 2 && sign == -1) {
        g = g4 * exp(Exponent);
        RadInt = n * n * g;
        return RadInt * RadInt * Cfac;
    }

    double A, B, C;
    int j;


    if (sign == 1) {
        for (int L = n - 3; L >= l; L--) {
            j = L + 2;
            A = 2. * n * sqrt((n * n - (j - 1) * (j - 1)) * (1. + j * j * k2));
            B = 4. * n * n - 4. * j * j + j * (2 * j - 1) * D;
            C = -2. * n * sqrt((n * n - j * j) * (1. + (j + 1) * (j + 1) * k2));
            g = (B * g2 + C * g1) / A;
            g1 = g2;
            g2 = g;
            if (fabs(g) > Chigh) {
                g *= Clow;
                g1 *= Clow;
                g2 *= Clow;
                Exponent += Ehigh;
            }
            else if (fabs(g) < Clow) {
                g *= Chigh;
                g1 *= Chigh;
                g2 *= Chigh;
                Exponent -= Ehigh;
            }
            //              printf("L,g,Exponent= %12i  %12.2e  %12.3f  \t\n",L,g,Exponent);
        }
    }
    else {
        for (int L = n - 3; L >= l; L--) {
            j = L + 1;
            A = 2. * n * sqrt((n * n - j * j) * (1. + (j - 1) * (j - 1) * k2));
            B = 4. * n * n - 4. * j * j + j * (2 * j + 1) * D;
            C = -2. * n * sqrt((n * n - (j + 1) * (j + 1)) * (1. + j * j * k2));
            g = (B * g4 + C * g3) / A;
            g3 = g4;
            g4 = g;
            if (fabs(g) > Chigh) {
                g *= Clow;
                g3 *= Clow;
                g4 *= Clow;
                Exponent += Ehigh;
            }
            else if (fabs(g) < Clow) {
                g *= Chigh;
                g3 *= Chigh;
                g4 *= Chigh;
                Exponent -= Ehigh;
            }
            //              printf("L,g,Exponent= %12i  %12.2e  %12.3f  \t\n",L,g,Exponent);
        }
    }

    g *= exp(Exponent);
    RadInt = n * n * g;

    // We have to multiply the output by 2/M_PI, due to the different wave func. normalizations;
    // Burgess uses a normalization of pi*delta(k^2-k'^2), and I use delta(E-E');
    // therefore my free wave function is sqrt(2./M_PI) times larger than the free wave func. of Burgess.

    return RadInt * RadInt * Cfac;
}

//////////////////////////////////////////////////////////////////////

  double SigmaPhotoionization(int n, int l, double omega, double E)
{

    // Ionization energy of the (n,l) state in atomic units:
    double Enl = 1. / (2. * n * n);

    // E is the energy of the outgoing electron in atomic units
     double omega = Enl + E;

    // Zero cross section below the ionization limit:
    if (omega < Enl)
        return 0.;

    // Outgoing electron energy in atomic units:
    //double E = omega - Enl;

    double Sum = 0.;
    for (int sign = -1; sign <= +1; sign += 2) {
        int lp = l + sign;  // angular momentum quantum number of final state
        if (lp < 0)
            continue;

        double lplmax;
        if (sign < 0)
            lplmax = l;
        else
            lplmax = l + 1;

        double M2 = RadInt2BoundFreeBurgess(n, l, E, sign);
        Sum += lplmax * M2;
    }

    return fC * omega / (2. * l + 1.) * Sum;
}

  /////////////////////////////////////////////////////////////////////


  /////////////////////////////////////////////////////////////////////
  double RadInt2Gordon(int n, int l, int np, int sign)
  {
  // Square of integral of RadialH(n,l,r)*RadialH(np,lp,r)*r*r*r from zero to infinity.
  // n,np,  l, lp=l+sign: principal and angular momentum quantum numbers.
  // lp=l+sign; sign has to be +1 or -1
  // If lp<0 or lp>np-1:  returns 0.  !!!
  // We use eq. 63.2 (page 262) in Bethe-Salpeter:Quantum mechanics of one and two electron atoms,
  //  which is for sign=-1.
  // This is known in the literature by Gordon formula (published by W. Gordon in 1929).
  // Integral RadialH(n,l,r)*RadialH(np,l+sign,r) =Integral RadialH(np,l+sign,r)*RadialH(n,l,r),
  // therefore for sign=+1 we change: np--> n, n--> np,  l+1 --> l
  // Maximal values of n and np: nMAX
     static double LogN[2*nMAX+10],LogNFactor[2*nMAX+10];
     static int start=0;
     if(start==0)
     {
  	for(int k=1;k<=2*nMAX+9;k++)
  	    LogN[k]=log((double)k);
  	LogN[0]=0.;
  	LogNFactor[0]=0;
  	LogNFactor[1]=0;
  	for(int k=2;k<=2*nMAX+9;k++)
  	    LogNFactor[k]=LogNFactor[k-1]+LogN[k];
  	start=1;
     }

     int lp=l+sign;
     if(lp<0 || lp>np-1)
       return 0.;

     if(n<1 || np<1 || l<0 || l>n-1)
     {
        puts("Message from function RadInt2Gordon:");
        puts("n and np have to be larger than 0!");
        puts("l has to be larger than -1 and smaller than n !");
        puts("Program running is stopped !!! ");
        exit(1);
     }

     if(n>nMAX || np>nMAX )
     {
        puts("Message from function RadInt2Gordon");
        puts("n and np cannot be larger than nMAX!");
        puts("Program running is stopped !!! ");
        exit(1);
     }

     if(sign!=-1 && sign!=+1)
     {
        puts("Message from function RadInt2Gordon:");
        puts("sign has to be either -1 or +1!");
        puts("Program running is stopped !!! ");
        exit(1);
     }

     if(sign==1)
     {
       int N,NP,L;
        N=np; NP=n; L=l+1;
        n=N; np=NP; l=L;
     }

     if(n==np)  // exception!
     {
         return 9./4.*n*n*(n*n-l*l);
     }

     int nr=n-l-1;
     int npr=np-l;

  // A, B, C, D, E0:
     double A=-LogN[4] - LogNFactor[2*l-1];
     double B=0.5*(LogNFactor[n+l] + LogNFactor[np+l-1] - LogNFactor[n-l-1] - LogNFactor[np-l]);
     double C=(l+1) * (LogN[4] + LogN[n] + LogN[np]);
     double D=(n+np-2*l-2) * LogN[abs(n-np)] - (n+np) * LogN[n+np];
     double E0=A+B+C+D;

  // x, del2, sum2:
     double del2=(n-np) * (n-np);
     double sum2=(n+np) * (n+np);
     double x=-4.*n*np/del2;

  // F1, F2 :
     double E1, E2, F1, F2;

     if(n>np)
     {
  	F1=HypergeometricFHoangBinh(-nr, -npr, 2*l, x, E1);
          F2=HypergeometricFHoangBinh(-nr-2, -npr, 2*l, x, E2);
     }
     else
     {
  	F1=HypergeometricFHoangBinh(-npr, -nr, 2*l, x, E1);
  	F2=HypergeometricFHoangBinh(-npr,-nr-2, 2*l, x, E2);
     }

  // RGordon:
     double RGordon=exp(E0+E1) * F1 - (del2/sum2) * exp(E0+E2) * F2;

  //  printf("E0,E1,E2,F1,F2= %10.2f %10.2f %10.2f  %10.2e %10.2e   \t\n",E0,E1,E2,F1,F2);


    return RGordon * RGordon;
  }

  transition_rate(double f,double FWHM, double I,int n,int l,int np, int sign)
  {
      // Rate of light source induced transition from state (n,l) to state (np,l+sign) in s^-1 .
      // sign has to be either -1 or +1.
      // n and np have to be different integers!
      // Light source with mean frequency f, Intensity I (W/m^2) and certain width (FWHM).
      const double BohrRadius=5.29177210903e-11, e_charge=1.602176634e-19, epsilon_zero=8.8541878128e-12, hbar=1.054571817e-34, light_speed=299792458, c=137.036,TimeAtomicUnit=2.4189e-17,k=3.1668e-6, h_eV=4.135667696e-15, h=6.62607015e-34,kb=1.380649e-23;

      double C=M_PI*pow(BohrRadius,2)*pow(e_charge,2)*(2*l+1)/pow(hbar,2)/epsilon_zero/light_speed;
      double En=-1./(n*n);
      double Enp=-1./(np*np);
      double fE=fabs(En-Enp)*13.598433/h_eV;

      double sigma=FWHM/2.3548;
      double nphoton=I/sqrt(2.*M_PI*pow(sigma,2))*exp(-((pow((f-fE),2))/(2*pow(sigma,2))));
      double CGC=getCGC(l,sign);

      //std::cout<<"Clebsch is: "<<CGC<<'\n';
      double P=nphoton*C*pow(CGC,2)*RadInt2Gordon(n,l,np,sign);
      //std::cout<<"P PSL is: "<<P<<'\n';
      return P;
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////////////////
  double transition_frequency(double I, double f, w_n1, w_n, FWHM)
   {
      double sigma=FWHM/2.3548;
      double freq_trans = w_n1 - w_n;
      double P_trans = (I/(sigma *sqrt(2*3.14)))*exp(-((freq_trans-f)**pow(2))/(2*sigma*sigma));
      return P_trans;
    }
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////
  double ionisation_light(double I, double FWHM, double f, int n, int l, double pos)
    {
      double P = (I * SigmaPhotoionization(n,l,omega))/(h_cross * f);
      return P;
    }
  ///////////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////////////////////////////////////////////////////////

}
