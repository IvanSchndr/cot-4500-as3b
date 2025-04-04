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
  
