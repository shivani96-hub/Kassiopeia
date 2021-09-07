//Calculating rates for Photon interaction
transition_rate(double f,double FWHM, double I,int n,int l,int np, int sign)
  {

double transition_frequency(double I, double f, w_n1, w_n, FWHM)
   {
      double sigma=FWHM/2.3548;
      double freq_trans = w_n1 - w_n;
      double P_trans = (I/(sigma *sqrt(2*3.14)))*exp(-((freq_trans-f)**pow(2))/(2*sigma*sigma));
      return P_trans;
    }
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
double ionisation_light(double I, double FWHM, double f, int n, int l, double pos)
    {
      double P = (I * SigmaPhotoionization(n,l,omega))/(h_cross * f);
      return P;
    }
////////////////////////////////////////////////////////////////////////
