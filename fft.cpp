#include <bits/stdc++.h>
using namespace std;

typedef complex<double> zz;

//fast fourier transform of a polynomial, (invert = true) gives the inverse transform
vector<zz> fft(vector<zz> const &pol, bool invert = false){
  int n = pol.size();
  vector<zz> prd (n, 0);
  
  if (n==1) {
    prd[0]=pol[0];
    return prd;
  }
  vector<zz> evpol(n/2, 0), oddpol(n/2, 0);
  for (int i = 0; i < n; i++){
    if (i % 2) oddpol[i/2] = pol[i];
    else evpol[i/2] = pol[i];
  }
  
  //classic divide and conquer
  evpol = fft(evpol, invert);
  oddpol = fft(oddpol, invert);
  
  
  double arg = 2*acos(-1) *(1-2*(!invert))/n;
  zz omega(cos(arg), sin(arg)), mtp(1, 0);
  
  for (int i = 0; i < n/2; i++, mtp*=omega){
    prd[i] = evpol[i] + mtp *oddpol[i];
  }
  mtp = 1;
  for (int i = 0; i < n/2; i++, mtp*=omega){
    prd[i+n/2] = evpol[i] - mtp*oddpol[i];
  }
  if (invert)
    for (int i = 0; i < n; i++)
      prd[i]/=2;
    
   // for (zz x : prd)
    // printf("%lf %lf \n", real(x), imag(x));
  // printf("\n");
  return prd;
}
    
vector<zz> multiply(vector<zz> const &f, vector<zz>const &g){//fast polynomial multiplication using fft
  int n;
  vector<zz> p(f.begin(), f.end()), q(g.begin(), g.end());
  for (n = 1; n < p.size() + q.size(); n<<=1);
  p.resize(n);
  q.resize(n);
  p=fft(p); q=fft(q);
  for (zz x : p)
    printf("%lf %lf\n", real(x), imag(x));
  printf("\n\n");
   for (zz x : q)
    printf("%lf %lf\n", real(x), imag(x));
  printf("\n\n");
  for (int i = 0; i < n; i++)
    p[i] *= q[i];
  p = fft(p, true);
  return p;
}
  
int main(){
  int df, dg;  double x;
  scanf("%d %d", &df, &dg);
  // scanf("%d", &df);
  vector<zz> f(df+1, 0), g(dg+1, 0);
  for (int i = 0; i <=df; i++){
    scanf("%lf", &x);
    f[i] = x;
  }
  for (int i = 0; i <=dg; i++){
    scanf("%lf", &x);
    g[i] = x;
  }
  
  vector<zz> h = multiply(f, g);
  
  for (zz x : h)
    printf("%lf %lf \n", real(x), imag(x));
  return 0;
}
  
  