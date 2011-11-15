#include "definitions.h"

 CustomWeakFormMassmatrix::CustomWeakFormMassmatrix(double time_step,Solution* sln_prev_time) : WeakForm(1) {
		CustomMatrixFormVolMassmatrix* mass_form= new CustomMatrixFormVolMassmatrix(0, 0, time_step);	
		add_matrix_form(mass_form);
		VectorFormVolMass* vector_form = new VectorFormVolMass(0, time_step);
		vector_form->ext.push_back(sln_prev_time);
		add_vector_form(vector_form);
  }
 CustomWeakFormMassmatrix::~CustomWeakFormMassmatrix(){
		delete get_mfvol()[0];
		delete get_vfvol()[0];
		WeakForm::delete_all();

	};


    template<typename Real, typename Scalar>
    Scalar CustomMatrixFormVolMassmatrix::matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, 
                       Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const {
		     Scalar result = 0; 
	  for (int i = 0; i < n; i++)
		result += wt[i] * (u->val[i] * v->val[i])/time_step;
	  return result;

    };

   scalar CustomMatrixFormVolMassmatrix::value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, 
                 Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const {
      return matrix_form<double, scalar>(n, wt, u_ext, u, v, e, ext);
    };

    Ord CustomMatrixFormVolMassmatrix::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
            Geom<Ord> *e, ExtData<Ord> *ext) const {
      return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
    };



    template<typename Real, typename Scalar>
    Scalar VectorFormVolMass::vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const {
	  Scalar result = 0;
	  Func<Scalar>* u_prev_time = ext->fn[0];
	  for (int i = 0; i < n; i++)
		result += wt[i] *  u_prev_time->val[i] * v->val[i]/time_step;
	  return result;

    };

   scalar VectorFormVolMass::value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const {
      return vector_form<double, scalar>(n, wt, u_ext, v, e, ext);
    };

    Ord VectorFormVolMass::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const {
      return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
    };



  CustomWeakFormConvection::CustomWeakFormConvection(Solution* sln_prev_time) : WeakForm(1) {
    add_matrix_form(new CustomMatrixFormVolConvection(0, 0));
   VectorFormVolConvection* vector_form = new VectorFormVolConvection(0);
    vector_form->ext.push_back(sln_prev_time);
    add_vector_form(vector_form);
  };
	CustomWeakFormConvection::~CustomWeakFormConvection(){
		delete get_mfvol()[0];
		delete get_vfvol()[0];
		WeakForm::delete_all();
	};




    template<typename Real, typename Scalar>
    Scalar CustomMatrixFormVolConvection::matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, 
                       Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const {

     Scalar result = 0; 
  for (int i = 0; i < n; i++)
    result += -wt[i] * (v->val[i] *(u->dx[i] * (0.5- e->y[i]) + u->dy[i] * (e->x[i]-0.5) ));
  return result;

    };

scalar CustomMatrixFormVolConvection::value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, 
                 Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const {
      return matrix_form<double, scalar>(n, wt, u_ext, u, v, e, ext);
    };

   Ord CustomMatrixFormVolConvection::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
            Geom<Ord> *e, ExtData<Ord> *ext) const {
      return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
    };

   



    template<typename Real, typename Scalar>
    Scalar VectorFormVolConvection::vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const {
  Scalar result = 0; 
  Func<Scalar>* u_prev_time = ext->fn[0];
  for (int i = 0; i < n; i++)
    result += -wt[i] * ( v->val[i] * (u_prev_time->dx[i] * (0.5- e->y[i]) + u_prev_time->dy[i] *  (e->x[i]-0.5)));
  return result;
    };
     scalar VectorFormVolConvection::value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const {
      return vector_form<double, scalar>(n, wt, u_ext, v, e, ext);
    };

     Ord VectorFormVolConvection::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const {
      return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
    };


//------------------------------------ANFANG fuer HERMES ohne FC-------------------------



 ConvectionForm::ConvectionForm( double tau, Solution* sln_prev_time) : WeakForm(1) {
    add_matrix_form(new ConvectionMatForm(0, 0,  tau));
    VectorConvection* vector_form = new VectorConvection(0,  tau);
    vector_form->ext.push_back(sln_prev_time);
    add_vector_form(vector_form);
  };

	ConvectionForm::~ConvectionForm(){
		delete get_mfvol()[0];
		delete get_vfvol()[0];
		WeakForm::delete_all();
	};

    template<typename Real, typename Scalar>
    Scalar   ConvectionMatForm::matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, 
                       Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const {

		       Scalar result = 0; 
		for (int i = 0; i < n; i++) 
		  result += wt[i] * (u->val[i] * v->val[i] / tau + 0.5*v->val[i] *(u->dx[i] * (0.5- e->y[i]) + u->dy[i] * (e->x[i]-0.5) ));	
		return result;

    };

scalar ConvectionMatForm::value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, 
                 Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const {
      return matrix_form<double, scalar>(n, wt, u_ext, u, v, e, ext);
    };
Ord ConvectionMatForm::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
            Geom<Ord> *e, ExtData<Ord> *ext) const {
      return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
    };


    


    template<typename Real, typename Scalar>
    Scalar VectorConvection::vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const {
			Scalar result = 0;
			Func<Scalar>* u_prev_newton = u_ext[0];
			Func<Scalar>* u_prev_time = ext->fn[0];
			for (int i = 0; i < n; i++) 
				result += wt[i] * ((u_prev_newton->val[i] - u_prev_time->val[i]) * v->val[i] / tau +
				                  0.5*v->val[i] * (u_prev_newton->dx[i] * (0.5- e->y[i]) + u_prev_newton->dy[i] *  (e->x[i]-0.5) +
					u_prev_time->dx[i] *  (0.5- e->y[i]) + u_prev_time->dy[i] *  (e->x[i]-0.5)));
	
			return result;

    };

  scalar VectorConvection::value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const {
      return vector_form<double, scalar>(n, wt, u_ext, v, e, ext);
    };

    Ord VectorConvection::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const {
      return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
    };



//---------------------------ENDE-------------

/* Initial condition */

 void CustomInitialCondition::derivatives(double x, double y, scalar& dx, scalar& dy) const {
      
    	double radius;
        //hump
	double x_0 =0.25;
	double y_0= 0.5;	
	radius = (1.0/0.15) * sqrt( pow((x-x_0),2.0) + pow((y-y_0),2.0));
	if( radius< 1.0) {		
		dx = -sin(radius*PI)/4.0*(1.0/(0.15 * sqrt( pow((x-x_0),2.0) + pow((y-y_0),2.0))))*2*x;
		dy = -sin(radius*PI)/4.0*(1.0/(0.15 * sqrt( pow((x-x_0),2.0) + pow((y-y_0),2.0))))*2*y;	
	}else{		
	//cone
		x_0 = 0.5;
		y_0 = 0.25;
		radius = 1.0/0.15 * sqrt( pow((x-x_0),2.0) + pow((y-y_0),2.0));
		if((radius< 1.0)&&(x!=x_0)) { 	
				dx = 1.0-(1.0/(0.15 * sqrt( pow((x-x_0),2.0) + pow((y-y_0),2.0))))*2*x;
			dy = 1.0-(1.0/(0.15 * sqrt( pow((x-x_0),2.0) + pow((y-y_0),2.0))))*2*y;	
		}else{dx=0.; dy=0.;}	
  }



};

 scalar CustomInitialCondition::value(double x, double y) const {
        scalar result = 0.0;
    	double radius;
        //hump
	double x_0 =0.25;
	double y_0= 0.5;	
	radius = (1.0/0.15) * sqrt( pow((x-x_0),2.0) + pow((y-y_0),2.0));
	if( radius<= 1.0) { 
		 result = (1.0+ cos(PI*radius))/4.0;
		return result;	
	}
	//slotted cylinder
	x_0 = 0.5;
	y_0 = 0.75;
	radius = 1.0/0.15 * sqrt( pow((x-x_0),2.0) + pow((y-y_0),2.0));
	if(radius <= 1) { 	
		if(fabs((x-x_0))>= 0.025) return 1.0;
		if(y>=0.85) return 1.0;
	}	
	//cone
	x_0 = 0.5;
	y_0 = 0.25;
	radius = 1.0/0.15 * sqrt( pow((x-x_0),2.0) + pow((y-y_0),2.0));
	if(radius<= 1.0) { 	
		result = 1.0-radius;
	}
       return result;
};

 Ord CustomInitialCondition::ord(Ord x, Ord y) const {
      return Ord(10);
};





