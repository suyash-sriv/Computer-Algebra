\\check the irreducibility of a polynomial in a finite field

irredp(pol, polvar, fieldvar) = 
  {
    d = poldegree(pol, polvar);
    flag = 1; 
    for (i = 1, d/2, 
      a = subst(pol, polvar, fieldvar);
      b = fieldvar^(d^i) - fieldvar;
      k = gcd(subst(pol, polvar, fieldvar), fieldvar^(d^i) - fieldvar);
      if (k != 1,flag = 0;));
    flag;
  }
  print(irredp(x^4+ x^3+x^2+x+1, x, ffgen(4, 'x)));
