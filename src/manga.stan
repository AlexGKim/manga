functions{
  vector lerp(vector x0, vector x1, vector y1){
    vector[size(x0)] ans;
    int ind1;

    ind1=0
    for (i in 1:size(ans)){
      while (ind1<size(x0) and x0[i]< x1[ind1+1]) {
        ind1=ind1+1
      }
      ans[i] = 
    }
    return sum((v1-mean(v1)) .* (v2-mean(v2)))/num_elements(v1);
  }
}

data {
  int nlam_m;
  real dlam;
  int nlam0;
  int nlam1;

  vector[nlam0] lam0;
  vector[nlam0] flux0;
  vector[nlam1] lam1;
  vector[nlam1] flux1;


  vector[nlam_m] lam_m;
  matrix[nlam0,nlam_m] K_0;
}

parameters {
  vector<lower=0.01,upper=100>[nlam_m] flux;

  real<lower=0.95,upper=1.05> a_1;
  real<lower=0> n_1;
  real<lower=0> sigma;
}


model {
  vector[nlam0] flux0_model;
  vector[nlam1] flux1_model;

  matrix[nlam1,nlam_m] K_1;
  real temp;

  for (i in 1:nlam1) {
    for (j in 1:nlam_m) {
        // temp = a_1*lam1[i]-lam_m[j];
        // if (fabs(temp) < 5*dlam){
          K_1[i,j] = a_1*lam1[i]-lam_m[j];
        // } else{
          // K_1[i,j] =0;
        // }
    }
  }
  K_1 = exp(-K_1 .* K_1/2./dlam^2);

  flux0_model = K_0 * flux;
  // print(K_0);
  // print(a_1);
  // print("flux0");
  // print(flux0_model);
  // print("flux");
  // print(flux);
  flux1_model = n_1*(K_1 * flux);
  // print("flux1")
  // print(flux1_model);

  target += normal_lpdf(flux0| flux0_model, sigma);
  target += normal_lpdf(flux1| flux1_model, sigma);

  target += cauchy_lpdf(sigma|0,1);
}