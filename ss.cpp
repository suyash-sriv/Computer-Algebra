#include "gcd.cpp"
using namespace std;
long long rse(long long base, long long expt, bool mod = false, long long modbase = 0){ //repeated squaring exponentiation
  if (mod) base%=modbase;
  long long mtp = base;
  long long ans = 1;
  while (expt){
    if (expt % 2) {ans *= mtp; if (mod) ans %= modbase;}
    expt /= 2; mtp *= mtp; if (mod) mtp %= modbase;
  }
  return ans;
}
long long jacobi(long long a, long long n){ //returns the jacobi symbol (a/n) defined as \Product{p^e totally divides n}{((a/p)^e)} where for prime p (a/p) is the Legendre symbol = {+-1, 0}, = a^(p-1)/2 (modulo(p))
  if (gcd(a, n).first != 1)
    return 0;
  long long mtp = 1;
  
  while (true){
    a = a % n; a = (2*a > n)?(a-n):a;
    if (a < 0){
      if (n%2 == 0){
        mtp *= 1;
        n/=2;
        a=-a;
      }
      else{
        mtp*=((n-1)/2 % 2)?-1:1;
        a = -a;
      }
    }
    while (a % 2 == 0 && a!= 0){
      a/=2;
      mtp *= (((n*n - 1)/8)%2)?-1:1;
    }
    if (a== 1){
      break;
    }
    if (a == 0) break;
    mtp *= ((((a-1)/2) * ((n-1)/2))%2)?-1:1;
    swap(a, n);
  }
  return mtp;
}
bool binarysearch(long long i, long long n){
  long long lo = 1, hi = round(pow(2, log(n)/(log(2)*i))) + 1;
  while (lo <= hi){
    //printf("(%d, %d) ", lo, hi);
    long long mid = (1ll*lo + hi)/2;
    long long pw = rse(mid, i);
    //printf("%lf;\n", pw);
    if ( pw > n)
      hi = mid - 1;
    else if (pw<n)
      lo = mid + 1;
    else return true;
  }
  return false;
}
bool isperfectpower(long long n){ //checks if n is a perfect power
  if (n <= 1) return true;
  double D = log(n)/log(2) + 1;
  //printf("%lf\n", D);
  for (long long i = 2; i <= D; i++)
    if (binarysearch(i, n)) return true;
  return false;
}

bool ss(long long n, int prec = 1){//solovay strassen test
  if (n == 2) return true;
  if (n%2== 0 || isperfectpower(n)) return false;
  long long a;
  
  while (prec--){
  while (a = rand()%(n-1) + 1)
    if (gcd(a, n).first == 1) break;
  int jc = (jacobi(a, n)+n)%n;
  //printf("%d %d %d \n", a, n, jc);
  if (jc != rse(a, (n-1)/2, true, n))
    return false;
  }
  return true;
}
pair<int, int> l2pd(int n){ // largest 2-power decomposition, i.e. returns k, m such that 2^k is the largest 2-power dividing n, and n = m*2**k
  int p = 0;
  while (!(n%2)){
    n/=2;
    p++;
  }
  return {p, n};
}
bool mr(long long n, int prec = 1){ //miller rabin test
  if (n == 2) return true;
  if (n%2 == 0 || isperfectpower(n))
    return false;
  long long a;
  while (prec--){
  while (a=rand()%(n-1) + 1)
    if (gcd(a, n).first == 1)
      break;
  if (rse(a, n-1, true, n) != 1) return false;
  auto x = l2pd(n-1);
  int k, m; k = x.first; m = x.second;
  int prev, curr;
  for (int i = 0; i < k; i++){
    curr = rse(a, m*(1<<i), true, n);
    prev = curr;
    if (curr == 1 && (prev != 1 && prev != -1))
      return false;
  }
  }
  return true;
}


//aks test was implemented in pari/gp

int main(){
  srand(time(0));
  long long x;
   long long a, n; scanf("%lld %lld", &a, &n);
  // printf("%d\n--\n", jacobi(a, n));
  // for (long long i = 1; i <= 2*a; i+=2, printf("\n"))
    // for (long long j = 1; j <= n; printf("%c ", (x = jacobi(j++, i), (x == 0)?'0':((x==1)?'+':'-' ))));
  
  // for (long long i = 0; i < 20; i++, printf("\n"))
    // for (long long j = 0; j < 20; j++)
      // printf("%lld ", rse(i, j));
  // printf("%lld", isperfectpower(128));
  // for (long long i = 0; i < n; i++)
    // if (isperfectpower(i))
      // printf("%d, ", i);
  
  for (int i = 0; i < 100; i++)
    if (mr(i, 5))
      printf("%d ", i);
    
}