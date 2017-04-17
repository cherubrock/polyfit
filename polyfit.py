from sympy import *

def make_poly_even(proto,x_arr,dx_arr):
  # generates a polynomial fitting given curve
  # proto = target expression/function eg. sin(x), the curve wil be fitted against function of x argument (x symbol is mandatory)
  # x_arr = list of x arguments values for which there will be exact match between proto and resulting polynomial
  # dx_arr = list of x arguments values for which there will be exact match between proto first derivative and resulting polynomial first derivative
  proto_diff = proto.diff()
  a = symbols(make_coeff_symbols(len(x_arr+dx_arr)))
  # insert zeros for odd coefficients
  coeffs = [val for pair in zip([0]*len(a),list(reversed(a))) for val in pair]
  poly = Poly(coeffs,gens=x)
  poly_diff = poly.diff()
  # prepare equations for solver
  eqns = []
  for pts in x_arr:
    cpoly = poly.subs(x,pts)
    cpoly = Add(cpoly,Mul(-1,proto.subs(x,pts)))
    eqns.append(cpoly);
  for pts in dx_arr:
    cpoly = poly_diff.subs(x,pts)
    cpoly = Add(cpoly,Mul(-1,proto_diff.subs(x,pts)))
    eqns.append(cpoly);
  #for eqn in eqns:
  #  print eqn
  solution = poly.subs(solve(eqns,a))
  #print solution
  return solution

def make_coeff_symbols(n,prefix='a'):
  symlist = []
  for i in range(n):
    symlist.append(prefix+str(i))
  return ','.join(symlist)


x = Symbol('x')
print make_poly_even(cos(pi*x),["0","1/4","1/2"],["1/4","1/2"])
print make_poly_even(cos(pi*x/2),["0","1"],["1"])


  