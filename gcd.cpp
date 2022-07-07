
//Extended Euclidean GCD Algorithm, 

#include <bits/stdc++.h>
using namespace std;

pair<int, int> divi(int a, int b){ //return quotient and remainder on dividing a by b, without sign issues (remainder is always in [-b/2..b/2])
  int mtp1 = 1, mtp2 = 1;
  if (a < 0){
    a = -a; mtp1 *= -1;
  }
  if (b < 0){
    b = -b; mtp2 *= -1;
  }
  int q = a/b;
  int r = a % b;
  if (2*r > b && r!= 0){
    r -=b; q++;
  }
  // printf("%d %d\n", q, r);
  return {q*mtp1 *mtp2, r*mtp1};
}

pair<int, pair<int, int>> gcd (int a, int b){ // returns {d, {u, v}} such that au + bv = d = (gcd(a, b)) and |u| < b and |v| < a
  int m = a, n = b; 
  if (a == 0) return {b, {0, 1}};
  int u = 1, v = 0, x = 0, y = 1, u_, v_, x_, y_, q, r, tmp; pair<int, int> dvsn;
  while (b!= 0){
    dvsn = divi(a, b);
    q = dvsn.first; r = dvsn.second;
    a = b; b = r;
    // printf("%d\n", r);
    u_ = x;        
    v_ = y;
    x_ = u - q*x;  
    y_ = v - q*y;
    u = u_; v = v_; x = x_; y = y_;
    // printf("%d %d : %d %d %d %d\n", a, b, u, v, x, y);
  }
  if (u >= n) {
    int q_ = u/n;
    u %= n;
    v += q_ * m;
  }
  if (v >= m) {
    int q_ = v/m;
    v %= m;
    u += q_ * n;
  }
  if (a < 0){
    a = -a;
    u = -u;
    v = -v;
  }
  return {a, {u, v}};
}

int inv(int a, int b){ //Compute the Multiplicative Inverse of a in Z/bZ
  auto x = gcd(a, b);
  assert(x.first == 1);
  int ainv = x.second.first;
  if (ainv < 0) ainv += b;
  return ainv;
}

// int main(){
  // int a, b;
  // scanf("%d %d", &a, &b);
  // pair<int, pair<int, int>> eq = gcd(a, b);
  // int d, u, v;
  // d = eq.first; pair<int, int> coeffs = eq.second;
  // u = coeffs.first; v = coeffs.second;
  // printf("The gcd of %d and %d is %d [= %d * (%d) + %d * (%d)]", a, b, d, a, u, b, v);
  // printf("The inverse of %d modulo %d is %d.", a, b, inv(a, b));
  // auto x = divi(a, b);
  // printf("%d %d", x.first, x.second);
// }