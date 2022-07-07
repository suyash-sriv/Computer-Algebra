lg(n)=
{
  return(  log(n)/log(2));
}


aks(n) =
{ 
my(R, k, r, a);
  if (n <= 1, return(0);,
  if (n==2, return (1);,
  if (n%2 == 0, return(0);,
    if (ispower(n), return(0);, \\if n is a perfect power, return 0, implemented in ss.cpp in (log(n))^2 time
    k = 4 * (lg(n))^2;
    R = ceil(16 *(lg(n) )^5); \\guaranteed to find a suitable r in 1..R (such that multiplicative order of n in Z/rZ is > k)
    for (tmp = 1, R, 
      r = tmp;
      \\print(r);
      a = iferr (znorder(Mod(n, r)), ERR, -1, errname(ERR) == "e_COPRIME"); //a is the multiplicative order (ord_r (n)),  ignore this r if n and r are not coprime
      if (a  > k, 
         break();););
  
  for (a = 1, min(r, n-1), 
    if (gcd(a, n) > 1, return(0);););
  t = ceil(2*sqrt(r)*lg(n));
  for (a = 1, t, 
    if (Mod((Mod(1, n)*x+Mod(a, n))^n - (Mod(1, n)*x^n +Mod(a, n)), x^r -1), return(0);););
  return(1);
    ););););
}
