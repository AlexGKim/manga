functions{

  // Linear interpolator  THIS IS NOT USED

  // The algorithm assumes that the target wavelengths are in order
  // and that the model wavelengths are in order

  vector lerp(vector x0, vector x1, vector y1){
    vector[dims(x0)[1]] ans;
    int ind1;

    ind1=1;
    for (i in 1:dims(ans)[1]){
      while (x1[ind1+1] < x0[i]) {
        ind1=ind1+1;
      }
      ans[i] = y1[ind1] + (y1[ind1+1]-y1[ind1])/(x1[ind1+1]-x1[ind1])*(x0[i]-x1[ind1]);
    }
    return ans;
  }

  // y = Kz
  // To make STAN faster, I don't do the full matrix multiplication, but rather only include terms
  // within 5 resolution elements of the target wavelength

  // The algorithm assumes that the target wavelengths are in order
  // and that the model wavelengths are in order

  vector kerp(vector lam0, vector lam_m, vector flux, real R){ 
    vector[dims(lam0)[1]] ans;
    real diff;
    real R2over2;
    real Rinv;
    int ind;
    int mark;
    real range;

    range=5.;
    R2over2=R*R/2;
    Rinv=1./R;

    mark=1; // marks the minimum possible model wavelength

    // Loop over all data wavelengths

    for (i in 1:dims(lam0)[1]) {
      ans[i]=0;

      // go to the lowest model index that contributes
      ind=mark;
      while(lam_m[ind+1] < lam0[i]-range*Rinv){
        ind = ind+1;
      }
      mark=ind; //since the data wavelengths are ordered 

      // do the matrix multiplication for the model terms that contribute
      while(lam_m[ind] < lam0[i]+range*Rinv){
        diff = lam_m[ind]-lam0[i];
        diff = diff * diff*R2over2;
        ans[i]=ans[i]+ exp(-diff)*flux[ind];
        ind = ind+1;
      }
    }
    return ans;
  } 
}

data {
  int nlam_m;
  int nlam0;
  int nlam1;

  real arange;
  real R;
  vector[nlam0] lam0;
  vector[nlam0] flux0;
  vector[nlam0] flux0_var;
  vector[nlam1] lam1;
  vector[nlam1] flux1;
  vector[nlam1] flux1_var;

  vector[nlam_m] lam_m;
}

parameters {
  vector<lower=-.2, upper=15>[nlam_m] flux; // The "z"'s in the model y = Kz
  simplex[2] norm;  // The normalization
  real<lower=-arange,upper=arange> scale_a;  // log-wavelength shift
  real<lower=0> sigma2; // variance
}

//  Model for data from spaxel i is
//  y_i ~ n_i Cauchy (K_i(ln lambda+scale_a_i) z, sigma)

model {
  vector[nlam0] flux0_model;
  vector[nlam1] flux1_model;

  flux0_model = 2*norm[1]*kerp(-0.5*scale_a+lam0,lam_m,flux,R);
  flux1_model = 2*norm[2]*kerp(0.5*scale_a+lam1,lam_m,flux,R);

  target += cauchy_lpdf(flux0| flux0_model, sqrt(flux0_var+sigma2));
  target += cauchy_lpdf(flux1| flux1_model, sqrt(flux1_var+sigma2));

  target += cauchy_lpdf(sigma2|0,.0009); // this is based on the dispersion from zero shift
}