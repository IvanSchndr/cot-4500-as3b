import numpy as np


def IVP(f, r, y0, h, e):
  t = np.arange(r[0], r[1] + h, h)
  y = np.zeros(len(t))
  y[0] = y0
  x = len(t) - 1
  for i in range(x):
    y[i + 1] = e(f, t[i], y[i], h)

  return y



def rk(f, tn, yn, h):
  k1 = f(tn, yn)
  k2 = f(tn + h / 2, yn + h * k1 / 2)  
  k3 = f(tn + h / 2, yn + h * k2 / 2)
  k4 = f(tn + h, yn + h * k3)
  return yn+h*(k1+2*k2+2*k3+k4)/6

def euler(f, tn, yn, h):
  return yn + h * f(tn, yn)
  

def f(t,y):
  return t-y*y


def gauss(a, b):
  n = len(b)
  m = n-1
  i = 0
  j = i-1
  x = np.zeros(n)

  new = np.concatenate((a, b), axis=1, dtype=float)


  while (i < n):
    if (new[i][i] != 0):
      for j in range(i + 1, n):
        scaling =( new[j][i] )/ (new[i][i])
        new[j] = new[j] - (scaling * new[i])
   
    i = i + 1
  
  x[m] = new[m][n] / new[m][m]

  for k in range(n - 2, -1, -1):
    x[k] = new[k][n]

    for j in range(k + 1, n):
      x[k] = x[k] - new[k][j] * x[j]

    x[k] = x[k] / new[k][k]
  print("\n[", end = " ")
  for ans in range(n):
    z = x[ans]
    print(f"{z:.2f}", end=" ")
  print("]", )

def l_u(a):
  n = len(a)
  i = 0
  j = i-1
  

  U = np.zeros((n,n))
  U = a.astype(float)
  L = np.zeros((n,n))
  L = L.astype(float)

  while (i < n):
    if (U[i][i] != 0):
      for j in range(i + 1, n):
        scaling =( U[j][i] )/ (U[i][i])
        L[j][i] = scaling
        U[j] = U[j] - (scaling * U[i])

    i = i + 1
  for k in range(n):
    L[k][k] = 1

  d = np.linalg.det(L) * np.linalg.det(U)
    
  print("\n", d, "\n\n", L, "\n\n", U)


def isdidominant(a):
  for i in range(len(a)):
    sum = 0
    for j in range(len(a)):
      if (i != j):
        sum += abs(a[i][j])
    if (abs(a[i][i]) <= sum):
      return False

  return True


def isposdef(a):
  for i in range(len(a)):
    for j in range(i, len(a)):
      if (a[i][j] != a[j][i]):
        return False
  ig = np.linalg.eigvals(a)
  for i in range(len(ig)):
    if (ig[i] < 0):
      return False
  return True

def main():
  r = [0, 2]
  iter = 10
  h = (r[1] - r[0]) / iter
  f0 = 1
  y = IVP(f, r, f0, h, euler)
  print(y[10])
  print("\n")
  y = IVP(f, r, f0, h, rk)
  print(y[10])
  print("")

  a = np.array([[2, -1, 1], [1, 3, 1], [-1, 5, 4]])
  at = np.array([[1, 1, 0, 3], [2, 1, -1, 1], [3, -1, -1, 2], [-1, 2, 3, -1]])
  b = np.array([[6], [0], [-3]])
 
  gauss(a, b)
  print("")

  l_u(at)

  dd = np.array([[9, 0, 5, 2, 1], [3, 9, 1, 2, 1], [0, 1, 7, 2, 3],[4, 2, 3, 12, 2],[3, 2, 4, 0, 8]])
  if isdidominant(dd):
    isdd = "YES"
  else:
    isdd = "NO"
  print("\n\n3. Is the matrix diagonally dominant:", isdd)
  pd = np.array([[2,2,1],[2,3,0],[1,0,2]])
  if isposdef(pd):
    ispd = "YES"
  else:
    ispd = "NO"
  print("\n4. Is the matrix a positive definite:", ispd)

if __name__ == "__main__":
  main()
