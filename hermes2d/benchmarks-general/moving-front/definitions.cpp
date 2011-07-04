#include "definitions.h"

double CustomExactSolution::value (double x, double y) const 
{
  double t = *t_ptr;
  return 0;
}

void CustomExactSolution::derivatives (double x, double y, scalar& dx, scalar& dy) const 
{
  double t = *t_ptr;
  
  dx = 0;
  dy = 0;
}

Ord CustomExactSolution::ord(Ord x, Ord y) const 
{
  return Ord(20);
}


double CustomFunction::value(double x, double y) const 
{
  double S = s;
  double C = c;
  double t = *t_ptr;
  double r2 = x*x + y*y;
  double r = sqrt(r2);
  double r3 = r2*r;
  double bubble = (x - x0)*(x - x1)*(y - y0)*(y - y1);
  double r_minus_t_squared = (-t + r)*(-t + r);
  double atan_term = (-atan(S*(-t + r)) + M_PI/2)/C;
  double S2 = s*s;
  double S3 = S2 * S;
  double x2 = x*x;
  double y2 = y*y;
  double S2_term = S2*r_minus_t_squared + 1;

  return - 2 * S3 * t * x2 * (-t + t1) * (-t + r) * bubble / (C * r2 * S2_term * S2_term) 
         - 2 * S3 * t * y2 * (-t + t1) * (-t + r) * bubble / (C * r2 * S2_term * S2_term) 
         - S * t * x2 * (-t + t1) * bubble / (C * r3 * S2_term) 
         + 2 * S * t * x * (-t + t1) * (x - x0) * (y - y0) * (y - y1) / (C * r * S2_term) 
         + 2 * S * t * x * (-t + t1) * (x - x1) * (y - y0) * (y - y1) / (C * r * S2_term) 
         - S * t * y2 * (-t + t1) * bubble / (C * r3 * S2_term) 
         + 2 * S * t * y * (-t + t1) * (x - x0) * (x - x1) * (y - y0) / (C * r * S2_term) 
         + 2 * S * t * y * (-t + t1) * (x - x0) * (x - x1) * (y - y1) / (C * r * S2_term) 
         + S * t * (-t + t1) * bubble / (C * S2_term) 
         + 2 * S * t * (-t + t1) * bubble / (C * r * S2_term) 
         - 2 * t * (-t + t1) * (x - x0) * (x - x1) * atan_term 
         - 2 * t * (-t + t1) * (y - y0) * (y - y1) * atan_term 
         - t * bubble * atan_term 
         + (-t + t1) * bubble * atan_term;
}

Ord CustomFunction::value(Ord x, Ord y) const 
{
  return Ord(20);
}
