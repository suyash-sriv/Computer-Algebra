//Matrix Multiplication using Strassen's algorithm, in O(n^lg(7)) time


#include<bits/stdc++.h>
using namespace std;

typedef vector<vector<int>> vvi;

void ppprint(vvi x){
  printf("%d \n", x[0][0]);
}


vvi add(vvi const &a, vvi const &b, int s1 = 1, int s2 = 1){
  int n = a.size();
  vvi sm (n, vector<int>(n,0));
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      sm[i][j] = a[i][j]*s1 + b[i][j]*s2;
  return sm;
}


vvi prod(vvi const &m1, vvi const &m2){
  int n = m1.size();
  vvi prd(n, vector<int> (n, 0));
  if (n==1){
    prd[0][0] = m1[0][0]*m2[0][0];
    return prd;
  }
  int k = n/2;
  vvi x11(k, vector<int> (k, 0)); vvi x12 = x11, x21 = x11, x22 = x11;
  vvi y11(k, vector<int> (k, 0)); vvi y12 = y11, y21 = y11, y22 = y11;
  for (int i = 0; i < k; i++)
    for (int j = 0; j < k; j++){
      x11[i][j] = m1[i][j];
      y11[i][j] = m2[i][j];
    }
  
  for (int i = 0; i < k; i++)
    for (int j = k; j < n; j++){
      x12[i][j-k] = m1[i][j];
      y12[i][j-k] = m2[i][j];
    }
  
  for (int i = k; i < n; i++)
    for (int j = 0; j < k; j++){
      x21[i-k][j] = m1[i][j];
      y21[i-k][j] = m2[i][j];
    }
  
  for (int i = k; i < n; i++)
    for (int j = k; j < n; j++){
      x22[i-k][j-k] = m1[i][j];
      y22[i-k][j-k] = m2[i][j];
    }
  
  vvi p1 = prod(add(x11, x22), add(y11, y22)), 
  p2 = prod(add(x21, x22), y11), 
  p3 = prod(x11, add(y12, y22, 1, -1)), 
  p4 = prod (x22, add(y11, y21, -1, 1)), 
  p5 = prod(add(x11, x12), y22), 
  p6 = prod(add(x11, x21, -1, 1), add(y11, y12)), 
  p7 = prod(add(x12, x22, 1, -1), add(y21, y22));
  
  vvi A = add(add(p1, p4), add(p5, p7, -1, 1)), 
  B = add(p3, p5), 
  C = add(p2, p4), 
  D = add(add(p1, p3), add(p2, p6, -1, 1));
  
  for (int i = 0; i < k; i++)
    for (int j = 0; j < k; j++)
      prd[i][j] = A[i][j];
  
  for (int i = 0; i < k; i++)
    for (int j = k; j < n; j++)
      prd[i][j]= B[i][j-k];
  
  for (int i = k; i < n; i++)
    for (int j = 0; j < k; j++)
      prd[i][j] = C[i-k][j];
  
  for (int i = k; i < n; i++)
    for (int j = k; j < n; j++)
      prd[i][j] = D[i-k][j-k];
    
  return prd;
}

int main(){
  
  // freopen("input.txt", "r", stdin);
  
  int n;
  scanf("%d", &n);
  vector<vector<int>> m1(n, vector<int> (n, 0)), m2(n,vector<int> (n, 0)), prd(n, vector<int> (n, 0));
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      scanf("%d", &m1[i][j]);
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      scanf("%d", &m2[i][j]);
  
  printf("Hello");
  prd = prod(m1, m2);
  
  for (int i = 0; i < n; i++){
    for (int j = 0; j < n; j++)
      printf("%d ", prd[i][j]);
    printf("\n");
  }
}