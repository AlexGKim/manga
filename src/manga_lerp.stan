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
}

data {
  int nlam_m;
  int nlam0;
  int nlam1;

  real amin;
  real amax;

  vector[nlam0] lam0;
  vector[nlam0] flux0;
  vector[nlam1] lam1;
  vector[nlam1] flux1;

  vector[nlam_m] lam_m;
}

parameters {
  vector<lower=0, upper=100>[nlam_m] flux;
  simplex[2] norm;
  // simplex[2] a_raw;
  real<lower=-1e-2,upper=1e-2> scale_a;
  real<lower=0> sigma;
}

model {
  vector[nlam0] flux0_model;
  vector[nlam1] flux1_model;


  flux0_model = 2*norm[1]*lerp(lam0, -0.5*scale_a+lam_m, flux);
  flux1_model = 2*norm[2]*lerp(lam1, 0.5*scale_a+lam_m, flux);


  target += normal_lpdf(flux0| flux0_model, sigma);
  target += normal_lpdf(flux1| flux1_model, sigma);

  target += cauchy_lpdf(sigma|0,.5);
}