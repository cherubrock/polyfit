from sympy import *

def make_poly(proto,x_arr,dx_arr,*args):
  # generates a polynomial approximatio of given curve
  # proto = target expression/function eg. sin(x), the curve wil be fitted against function of x argument (x symbol is mandatory)
  # x_arr = list of x arguments values for which there will be exact match between proto and resulting polynomial
  # dx_arr = list of x arguments values for which there will be exact match between proto first derivative and resulting polynomial first derivative
  proto_diff = proto.diff()
  a = symbols(make_coeff_symbols(len(x_arr+dx_arr)))
  # insert zeros for odd/even coefficients
  if 'even' in args:
    coeffs = [val for pair in zip([0]*len(a),list(reversed(a))) for val in pair]
  elif 'odd' in args:
    coeffs = [val for pair in zip(list(reversed(a)),[0]*len(a)) for val in pair]
  else:
    coeffs = a
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
  solution = poly.subs(solve(eqns,a))
  if 'debug' in args:
    for eqn in eqns:
      print str(eqn) + " = 0"
    print solution
  return solution

def print_metrics(proto,poly_approx,x_range,*args):
  # proto = target expression/function eg. sin(x), (x symbol is mandatory)
  # poly_approx = polynomial to test against prototype
  # x_range = [min,max] 2 element array
  signal = Integral(poly_approx, (x,x_range[0],x_range[1])).as_sum(10,'trapezoid')
  dev_fun = sqrt((proto - poly_approx)**2)
  deviation = Integral(dev_fun, (x,x_range[0],x_range[1])).as_sum(10,'trapezoid')
  snr = (signal/deviation)
  snrdb = '%0.1f' % (20*log(snr,10))
  enob = '%0.1f' % (log(snr,2))
  print "SNR: " + str(snrdb) + " [dB]"
  print "ENOB: " + str(enob)
  return 0

def make_coeff_symbols(n,first_letter='a'):
  first = ord(first_letter)
  return ','.join(chr(i) for i in range(first,first+n))


x = Symbol('x')
# calculate symbolic
poly_approx = make_poly(sin(pi*x/2),[1],[1,0],'odd','debug')
# show coeffs as float
print poly_approx.evalf()
# show SNR and ENOB
print_metrics(sin(pi*x/2),poly_approx,[0,1])

#print_metrics(cos(pi*x),make_poly(cos(pi*x),[0,"1/4","1/2"],["1/4","1/2"],'even','debug'),[0,"1/2"])


  