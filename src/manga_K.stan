functions{
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

  vector kerp(vector lam0, vector lam_m, vector flux, real R){ 
    vector[dims(lam0)[1]] ans;
    // matrix[dims (lam0)[1],dims(lam_m)[1]] K_0;
    real diff;
    real R2over2;
    real Rinv;
    int ind;
    int mark;
    // vector[dims(lam_m)[1]] diff;

    R2over2=R*R/2;
    Rinv=1./R;
    mark=1;
    for (i in 1:dims(lam0)[1]) {
      ans[i]=0;
      ind=mark;
      // print(lam_m[ind]," ",lam0[i]);
      while(lam_m[ind+1] < lam0[i]-5*Rinv){
        // print(lam_m[ind]," ",lam0[i]);
        ind = ind+1;
      }

      mark=ind;
      while(lam_m[ind] < lam0[i]+5*Rinv){
        // print(lam_m[ind]," ",lam0[i]);
        diff = lam_m[ind]-lam0[i];
        diff = diff * diff*R2over2;
        ans[i]=ans[i]+ exp(-diff)*flux[ind];
        ind = ind+1;
      }
      // print(lam_m[ind]," ",lam0[i]);
      // diff = lam_m-lam0[i];
      // diff = diff .* diff*R2over2;
      // for (j in 1:dims(lam_m)[1]) {
      //     if (fabs(diff[j]) < 50){
      //       ans[i]=ans[i]+ exp(-diff[j])*flux[j];
      //     } 
      // }
    }
    // K_0 = exp(-K_0 .* K_0/2.*R2);
    // return K_0*flux;
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
  vector[nlam1] lam1;
  vector[nlam1] flux1;

  vector[nlam_m] lam_m;

}

parameters {
  vector<lower=-.2, upper=15>[nlam_m] flux;
  simplex[2] norm;
  // simplex[2] a_raw;
  real<lower=-arange,upper=arange> scale_a;
  real<lower=0> sigma;
}

model {
  vector[nlam0] flux0_model;
  vector[nlam1] flux1_model;

  // matrix[nlam0,nlam_m] K_0;
  // matrix[nlam1,nlam_m] K_1;

  // for (i in 1:nlam0) {
  //   for (j in 1:nlam_m) {
  //         K_0[i,j] = -0.5*scale_a+lam0[i]-lam_m[j];
  //   }
  // }
  // K_0 = exp(-K_0 .* K_0/2.*R2);
   
  // for (i in 1:nlam1) {
  //   for (j in 1:nlam_m) {
  //         K_1[i,j] = 0.5*scale_a+lam1[i]-lam_m[j];
  //   }
  // }

  // K_1 = exp(-K_1 .* K_1/2.*R2);

  flux0_model = 2*norm[1]*kerp(-0.5*scale_a+lam0,lam_m,flux,R);
  flux1_model = 2*norm[2]*kerp(0.5*scale_a+lam1,lam_m,flux,R);

  target += cauchy_lpdf(flux0| flux0_model, sigma);
  target += cauchy_lpdf(flux1| flux1_model, sigma);

  target += cauchy_lpdf(sigma|0,.03); // this is based on the dispersion from zero shift
}