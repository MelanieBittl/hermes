#include "weakform/weakform.h"
#include "integrals/h1.h"
#include "boundaryconditions/essential_bcs.h"

 // CustomWeakFormConvectionLinear wf(time_step, &u_prev_time);

class CustomWeakFormConvectionLinear : public WeakForm
{
public:
  CustomWeakFormConvectionLinear( double tau, Solution* sln_prev_time) : WeakForm(1) {
    add_matrix_form(new CustomMatrixFormVolConvectionLinear(0, 0,  tau));
    VectorFormVolConvection* vector_form = new VectorFormVolConvection(0,  tau);
    vector_form->ext.push_back(sln_prev_time);
    add_vector_form(vector_form);
  };

private:
  class CustomMatrixFormVolConvectionLinear : public WeakForm::MatrixFormVol
  {
  public:
    // This weak form is custom since it contains a nonlinearity in the diffusion term.
    CustomMatrixFormVolConvectionLinear(int i, int j, double tau) 
      : WeakForm::MatrixFormVol(i, j, HERMES_ANY, HERMES_NONSYM), tau(tau) { }

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, 
                       Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const {

         Scalar result = 0; 
  for (int i = 0; i < n; i++)
    result += wt[i] * (u->val[i] * v->val[i] / tau + 0.5*v->val[i] *(u->dx[i] * (0.5- e->y[i]) + u->dy[i] * (e->x[i]-0.5) ));
  return result;

    }

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, 
                 Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const {
      return matrix_form<double, scalar>(n, wt, u_ext, u, v, e, ext);
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
            Geom<Ord> *e, ExtData<Ord> *ext) const {
      return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
    }


    
    // Members.
  
    double tau;
  };

  // This form (residual) is custom since it contains a nonlinear term.
  class VectorFormVolConvection : public WeakForm::VectorFormVol
  {
  public:
    VectorFormVolConvection(int i,double tau) : WeakForm::VectorFormVol(i), tau(tau) { }

    template<typename Real, typename Scalar>
    Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const {
  Scalar result = 0;
  Func<Scalar>* u_prev_newton = u_ext[0];
  Func<Scalar>* u_prev_time = ext->fn[0];
  for (int i = 0; i < n; i++)
    result += wt[i] * ((u_prev_newton->val[i] - u_prev_time->val[i]) * v->val[i] / tau +
                      0.5*v->val[i] * (u_prev_newton->dx[i] * (0.5- e->y[i]) + u_prev_newton->dy[i] *  (e->x[i]-0.5) +
			u_prev_time->dx[i] *  (0.5- e->y[i]) + u_prev_time->dy[i] *  (e->x[i]-0.5)));
  return result;

    }

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const {
      return vector_form<double, scalar>(n, wt, u_ext, v, e, ext);
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const {
      return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
    }


    // Members.
    
    double tau;
  };
};



/* Initial condition */

class CustomInitialCondition : public ExactSolutionScalar
{
public:
  CustomInitialCondition(Mesh* mesh) : ExactSolutionScalar(mesh) {};

  virtual void derivatives (double x, double y, scalar& dx, scalar& dy) const {
      
    	double radius;
        //hump
	double x_0 =0.25;
	double y_0= 0.5;	
	radius = (1.0/0.15) * sqrt( pow((x-x_0),2.0) + pow((y-y_0),2.0));
	if( radius< 1) {		
		dx = -sin(radius*PI)/4.0*(1.0/(0.15 * sqrt( pow((x-x_0),2.0) + pow((y-y_0),2.0))))*2*x;
		dy = -sin(radius*PI)/4.0*(1.0/(0.15 * sqrt( pow((x-x_0),2.0) + pow((y-y_0),2.0))))*2*y;	
	/*}else{			
	//cone
	x_0 = 0.5;
	y_0 = 0.25;
	radius = 1.0/0.15 * sqrt( pow((x-x_0),2.0) + pow((y-y_0),2.0));
	if((radius< 1)&&(x!=x_0)) { 	
			dx = 1-(1.0/(0.15 * sqrt( pow((x-x_0),2.0) + pow((y-y_0),2.0))))*2*x;
		dy = 1-(1.0/(0.15 * sqrt( pow((x-x_0),2.0) + pow((y-y_0),2.0))))*2*y;	*/
	}else{dx=0.; dy=0.;}	
      // }



  };

  virtual scalar value (double x, double y) const {
           double result = 0.0;
    	double radius;
        //hump
	double x_0 =0.25;
	double y_0= 0.5;	
	radius = (1.0/0.15) * sqrt( pow((x-x_0),2.0) + pow((y-y_0),2.0));
	if( radius<= 1) { 
		 result = (1.0+ cos(PI*radius))/4.0;
		return result;	
	}
	//slotted cylinder
/*	x_0 = 0.5;
	y_0 = 0.75;
	radius = 1.0/0.15 * sqrt( pow((x-x_0),2.0) + pow((y-y_0),2.0));
	if(radius <= 1) { 	
		if(fabs((x-x_0))>= 0.025) return 1.0;
		if(y>=0.85) return 1.0;
	}*/
	
	//cone
/*	x_0 = 0.5;
	y_0 = 0.25;
	radius = 1.0/0.15 * sqrt( pow((x-x_0),2.0) + pow((y-y_0),2.0));
	if(radius<= 1) { 	
		result = 1-radius;
	}*/

       return result;
  }

  virtual Ord ord(Ord x, Ord y) const {
      return x+y;
  }
};
