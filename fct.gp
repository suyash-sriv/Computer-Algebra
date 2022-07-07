ddf(pol, polvar, fieldvar) = \\distinct degree factorization, factorize f into Product_{1 <= i <= d} f_i where each f_i is the product of all  degree i irreducible polyomials in the finite field that divide f
  {
    L = List();
    d = poldegree(pol, polvar);
    flag = 1;
    i = 1; 
    while (i < d, 
      a = subst(pol, polvar, fieldvar);
      b = fieldvar^(d^i) - fieldvar;
      k = gcd(subst(pol, polvar, fieldvar), fieldvar^(d^i) - fieldvar);
      if (k != 1,flag = 0;listput(L, k); pol = pol/k;, i=i+1));
   print(L); 
    flag;
  }
\\ddf(x^4+ x^3+x^2+x+1, x, ffgen(4, 'x)); \\example run of ddf factorization in F_4

auxedf(pol, dg, polvar, fieldvar)=  \\helper function for edf , splits pol into two smaller factors with probability > 1/2, is called repeatedly by edf till a split is found, upon which divide and conquer does the job for us
{
  n = poldegree(pol);
  a =Pol(vector(2, i, random(fieldvar)), 'x);
  if (poldegree(a) == 0, -1;);
  g = gcd(a, pol);
  if (g!=1, g;);
  b = a^(((fieldvar.p*fieldvar.f)^dg -1)/2) ;
  b = b%pol;
  g2 = gcd(b-1, pol);
  if (g2!=1 && g2!=pol, g2;,-1);
}

edf(pol, dg, polvar, fieldvar)= \\equal degree factorization of a polynomial with all factors of known equal degree dg
{
  n = poldegree(pol);
  if (n==dg, vector(1,i,pol);,
  a = -1;
  while (a == -1,
  a = auxedf(pol, dg, polvar, fieldvar););
  L = edf(a, dg, polvar, fieldvar); //recursively factorize the factors
  K = edf(pol/a, dg, polvar, fieldvar);
  v = vector(#L+#K, i, if (i <= #L, L[i], K[i-#L]););//merge the results
  v;);
}
\\print(edf(x^4+ x^3 + x - 1, 2, x, ffgen(3, 'x))); \\example run of edf in F_3


squarefreepart(pol, polvar, fieldvar) = \\returns squarefreepart of the polynomial
{
  my(g, d, v, k, t, v2);
  if (poldegree(pol) < 1, pol;,
  g=Mod(deriv(pol,  polvar), fieldvar.p);
  d = Mod(gcd(pol, g), fieldvar.p);  
  print(pol);
  print(g);  print(d);
  if (g==0, ppol = takeroot(pol, polvar, fieldvar); squarefreepart(ppol, polvar, fieldvar);,
  v = pol/d;
print(v);
  k = gcd(d, v^(poldegree(pol)));
print(k);
  t = takeroot(d/k, polvar, fieldvar);
  print(t);
  v2 = squarefreepart(t, polvar, fieldvar);
  v*v2;
  ););
  
}
   
   

\\in a finite field with characteristic p (>0, ofcourse), derivative of a non-constant polynomial may be zero (!), (iff the polynomial is a p^th power(!!)). function takeroot, takes the p^th root of such a polynomial.

takeroot(pol, polvar, fieldvar)=
{
  d = poldegree(pol, polvar);
  p = fieldvar.p; f = fieldvar.f; \\print(d, p ,f);
  Polrev(vector(d/p+1, i, polcoef(pol, p*(i-1) , polvar)^(p^(f - 1))));
}



\\some example runs
\\print(squarefreepart(x^4, x, ffgen(2, 'x)));

\\print(squarefreepart(x^4*(x+1)^2*(x-1)^2*(x^2+1)^2*(x^2+x+1), x, ffgen(2, 'x)));

\\takeroot(x^4 +x ^2, x, ffgen(4, 'x));

