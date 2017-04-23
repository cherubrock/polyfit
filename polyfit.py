from sympy import *

def make_poly(proto,x_arr,dx_arr,*args):
  # generates a polynomial approximatio of given curve
  # proto = target expression/function eg. sin(x), the curve wil be fitted against function of x argument (x symbol is mandatory)
  # x_arr = list of x arguments values for which there will be exact match between proto and resulting polynomial
  # dx_arr = list of x arguments values for which there will be exact match between proto first derivative and resulting polynomial first derivative
  # *args => 'odd' = hint to solver that prototype is odd function
  # *args => 'even' = hint to solver that prototype is even function
  # *args => 'debug' = print equations to solve
  # *args => 'zero-average' = make additional equation to satisfy zero error average over x range

  proto_diff = proto.diff()
  poly_order = len(x_arr+dx_arr)
  x_range = [min(x_arr+dx_arr),max(x_arr+dx_arr)]
  if 'zero-average' in args:
    poly_order+=1
  a = symbols(make_coeff_symbols(poly_order))
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
  if 'zero-average' in args:
    cpoly = integrate(poly,(x,x_range[0],x_range[1])).as_expr()
    cpoly = Add(cpoly,Mul(-1,integrate(proto, (x,x_range[0],x_range[1]))))
    eqns.append(cpoly);
  solution = poly.subs(solve(eqns,a))
  if 'debug' in args:
    for eqn in eqns:
      print "Equation: " + str(eqn) + " = 0"
    print "Solution: " + str(solution)
  return solution

def print_metrics(proto,poly_approx,x_range,*args):
  # proto = target expression/function eg. sin(x), (x symbol is mandatory)
  # poly_approx = polynomial to test against prototype
  # x_range = [min,max] 2 element array
  signal = Integral(proto, (x,x_range[0],x_range[1])).as_sum(4,'trapezoid')
  dev_fun = sqrt((proto - poly_approx)**2)
  deviation = Integral(dev_fun, (x,x_range[0],x_range[1])).as_sum(4,'trapezoid')
  snr = signal/deviation
  snrdb = '%0.1f' % (20*log(snr,10))
  enob = '%0.1f' % (log(snr,2))
  print "SNR: " + str(snrdb) + " [dB]"
  print "ENOB: " + str(enob)

def make_coeff_symbols(n,first_letter='a'):
  first = ord(first_letter)
  return ','.join(chr(i) for i in range(first,first+n))

def make_poly_all(proto,x_arr,dx_arr,*args):
  solution = make_poly(proto,x_arr,dx_arr,*args)
  print "Float Eval: " + str(solution.evalf())
  x_range = [min(x_arr+dx_arr),max(x_arr+dx_arr)]
  print_metrics(proto,solution,x_range)

# usage
x = Symbol('x')
#make_poly_all(sin(pi*x/2),[1],[0,1],'odd','zero-average','debug')
#make_poly_all(cos(pi*x),["0","1/2"],["1/2"],'even','zero-average','debug')
make_poly_all(sin(pi*x),["1/2"],[0,"1/2"],'odd','zero-average','debug')
make_poly_all(sin(pi*x),["1/4","1/2"],[0,"1/2"],'odd','debug')