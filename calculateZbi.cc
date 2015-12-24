double calculateZbi(double signal=1.1, double bkg=5.5, double unc=0.3){

  double n_on = signal+bkg;
  double mu_b_hat=bkg;
  double sigma_b=unc*bkg;
  double tau = mu_b_hat/(sigma_b*sigma_b);
  double n_off = tau*mu_b_hat;
  double P_Bi = TMath::BetaIncomplete(1./(1.+tau),n_on,n_off+1);
  double Z_Bi = sqrt(2)*TMath::ErfInverse(1 - 2*P_Bi);
  return Z_Bi;
  //  std::cout<<"The calculated Zbi for a signal of "<<signal<<" events and background of "<<bkg<<" events with a systematic uncertainty of "<<unc*100<<"% is "<<Z_Bi<<std::endl;
 
}
