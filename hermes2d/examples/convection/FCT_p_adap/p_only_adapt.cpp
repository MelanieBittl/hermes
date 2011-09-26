#include "p_only_adapt.h"


bool PonlyAdapt::adapt(int* elements_to_refine)
{  info("p-adapt");
	if(elements_to_refine==NULL) return false;
	  //get meshes
  Mesh* meshes[H2D_MAX_COMPONENTS];
  for (int j = 0; j < this->num; j++) {
    meshes[j] = this->spaces[j]->get_mesh();
  }
	bool changed = false;
 Space* space = this->spaces[0];
int order, v_ord, h_ord;
  //apply refinements
	Element* e;
	for_all_active_elements(e, space->get_mesh()){
		if(elements_to_refine[e->id] == 2) {			//refine
			if(e->is_triangle()==true){
				order = space->get_element_order(e->id); 
				if(order >=1){
					space->set_element_order_internal(e->id, order+1);
					changed = true;
				}
			}else if(e->is_quad()==true){
				v_ord = H2D_GET_V_ORDER(space->get_element_order(e->id));
				h_ord = H2D_GET_H_ORDER(space->get_element_order(e->id));
				order = H2D_MAKE_QUAD_ORDER(h_ord+1, v_ord+1);
					space->set_element_order_internal(e->id, order);
					changed = true;
				}
								
			}else if(elements_to_refine[e->id] == 0) {			//coarse
			if(e->is_triangle()==true){
				order = space->get_element_order(e->id); 
				if(order >1){
					space->set_element_order_internal(e->id, order-1);
					changed = true;
				}
			}else if(e->is_quad()==true){
				v_ord = H2D_GET_V_ORDER(space->get_element_order(e->id));
				h_ord = H2D_GET_H_ORDER(space->get_element_order(e->id));
				order = H2D_MAKE_QUAD_ORDER(h_ord-1, v_ord-1);
					if((v_ord >1)&&(h_ord>1)){
						space->set_element_order_internal(e->id, order);
						changed = true;
					}
				}
								
			}		
		
		}
	  
	if(changed==false){ info("nothing to refine/coarse");return false;}

  // in singlemesh case, impose same orders across meshes
 // homogenize_shared_mesh_orders(meshes);

  // since space changed, assign dofs:
  Space::assign_dofs(this->spaces);

  return true;
}



double ErrorEstimation::calc_err_extrapolation(Solution *sln, Solution *rsln, bool solutions_for_adapt,
                    unsigned int error_flags )
{
return calc_err_internal(sln, rsln, NULL, solutions_for_adapt, error_flags);

};

double ErrorEstimation::calc_l2_norm(Solution* sln, Solution* rsln,
                             bool solutions_for_adapt,
                               unsigned int error_flags){

Hermes::vector<Solution *> slns;
slns.push_back(sln);
Hermes::vector<Solution *> rslns;
rslns.push_back(rsln);
return calc_l2_norm(slns, rslns,  solutions_for_adapt, error_flags);

};

double ErrorEstimation::l2_error(MeshFunction *sln1, MeshFunction *sln2, MeshFunction *rsln1,
                     MeshFunction *rsln2, Space* space,bool solutions_for_adapt, MatrixSolverType solver_type){
if (space == NULL) error("Space == NULL in  ErrorEstimation::l2_error");
	double result =0;
	int ndof = space->get_num_dofs();
	scalar* coeffs= new scalar[ndof];
	OGProjection::project_global(space, sln1, coeffs, solver_type, HERMES_L2_NORM);
	result = this->l2_error(sln1, sln2, rsln1,rsln2, coeffs,space,solutions_for_adapt);
	delete [] coeffs;
	return result;

};


double ErrorEstimation::calc_l2_norm(Hermes::vector<Solution *> slns, Hermes::vector<Solution *> rslns,
                                 bool solutions_for_adapt, unsigned int error_flags)
{

 _F_
  int i, j, k;

  int n = slns.size();
  if (n != this->num) EXIT("Wrong number of solutions.");

  TimePeriod tmr;

  Solution* rslns_original[H2D_MAX_COMPONENTS];
  Solution* slns_original[H2D_MAX_COMPONENTS];

  for (i = 0; i < n; i++) {
    slns_original[i] = this->sln[i];
    this->sln[i] = slns[i];
    sln[i]->set_quad_2d(&g_quad_2d_std);
  }
  for (i = 0; i < n; i++) {
    rslns_original[i] = this->rsln[i];
    this->rsln[i] = rslns[i];
    rsln[i]->set_quad_2d(&g_quad_2d_std);
  }

  have_coarse_solutions = true;
  have_reference_solutions = true;

  // Prepare multi-mesh traversal and error arrays.
  Mesh **meshes = new Mesh *[2 * num];
  Transformable **tr = new Transformable *[2 * num];
  Traverse trav;
  num_act_elems = 0;
  for (i = 0; i < num; i++) {
    meshes[i] = sln[i]->get_mesh();
    meshes[i + num] = rsln[i]->get_mesh();
    tr[i] = sln[i];
    tr[i + num] = rsln[i];

    num_act_elems += sln[i]->get_mesh()->get_num_active_elements();

    int max = meshes[i]->get_max_element_id();
    if(solutions_for_adapt) {
      if (errors[i] != NULL) delete [] errors[i];
      errors[i] = new double[max];
      memset(errors[i], 0, sizeof(double) * max);
    }
  }

  double total_norm = 0.0;
  double *norms = new double[num];
  memset(norms, 0, num * sizeof(double));
 
  if(solutions_for_adapt) this->errors_squared_sum = 0.0;
  double total_error = 0.0;

  // Calculate error.
  Element **ee;
  trav.begin(2 * num, meshes, tr);
  while ((ee = trav.get_next_state(NULL, NULL)) != NULL) {
    for (i = 0; i < num; i++) {
      for (j = 0; j < num; j++) {
        if (error_form[i][j] != NULL) {
          double err, nrm;
          err = eval_error(error_form[i][j], sln[i], sln[j], rsln[i], rsln[j]);
          nrm = eval_error_norm(error_form[i][j], rsln[i], rsln[j]);

          norms[i] += nrm;
          total_norm  += nrm;
          total_error += err;
  
          if(solutions_for_adapt)
            this->errors[i][ee[i]->id] += err;
        }
      }
    }
  }
  trav.finish();



  tmr.tick();
  error_time = tmr.accumulated();

  // Make the error relative if needed.
  if(solutions_for_adapt) {
    if ((error_flags & HERMES_ELEMENT_ERROR_MASK) == HERMES_ELEMENT_ERROR_REL) {
      for (int i = 0; i < this->num; i++) {
        Element* e;
        for_all_active_elements(e, meshes[i])
          errors[i][e->id] /= norms[i];
      }
    }

    this->errors_squared_sum = total_error;

    // Element error mask is used here, because this variable is used in the adapt()
    // function, where the processed error (sum of errors of processed element errors)
    // is matched to this variable.
    if ((error_flags & HERMES_TOTAL_ERROR_MASK) == HERMES_ELEMENT_ERROR_REL)
      errors_squared_sum = errors_squared_sum / total_norm;
  }

  // Prepare an ordered list of elements according to an error.
  if(solutions_for_adapt) {
    fill_regular_queue(meshes);
    have_errors = true;
  }
  else {
    for (i = 0; i < n; i++) {
      this->sln[i] = slns_original[i];
      this->rsln[i] = rslns_original[i];
    }
  }

  delete [] meshes;
  delete [] tr;
  delete [] norms;



  // Return error value.
/*
  if ((error_flags & HERMES_TOTAL_ERROR_MASK) == HERMES_TOTAL_ERROR_ABS)
   return sqrt(total_error);
  else if ((error_flags & HERMES_TOTAL_ERROR_MASK) == HERMES_TOTAL_ERROR_REL)
    return sqrt(total_error / total_norm);
  else {
    error("Unknown total error type (0x%x).", error_flags & HERMES_TOTAL_ERROR_MASK);
    return -1.0;
  }*/
	return sqrt(total_norm);
}







double  ErrorEstimation::l2_error(MeshFunction *sln1, MeshFunction *sln2, MeshFunction *rsln1,
                         MeshFunction *rsln2, scalar* coeffs, Space* space, bool solutions_for_adapt){
	if(elem_ex==NULL) error("elem_ex==NULL in ErrorEstimation::l2_error" );

	Adapt::MatrixFormVolError* form =error_form[0][0];
	Mesh* mesh = space->get_mesh();	
	PrecalcShapeset* pss = new PrecalcShapeset(space->get_shapeset());

	RefMap *rv1 = sln1->get_refmap();
 	 RefMap *rv2 = sln2->get_refmap();
  	RefMap *rrv1 = rsln1->get_refmap();
  	RefMap *rrv2 = rsln2->get_refmap();

	RefMap* refmap_ex = new RefMap;
	refmap_ex->set_active_element(elem_ex);
	 

	Element* el;
	double total_error = 0.0;
	AsmList* al_ex = new AsmList;
	space->get_element_assembly_list(elem_ex, al_ex); 

	int nc = 0;
	int max = mesh->get_max_element_id();
	if(solutions_for_adapt) {
      if (errors[nc] != NULL) delete [] errors[nc];
      errors[nc] = new double[max];
      memset(errors[nc], 0, sizeof(double) * max);
    }

 for_all_active_elements(el, mesh) 
 {	 
	Quad2D* quad_e = sln1->get_quad_2d();   
    quad_e->set_mode(el->get_mode());
   	int o = space->get_element_order(el->id);
	 o = std::max(H2D_GET_H_ORDER(o), H2D_GET_V_ORDER(o));
    int np = quad_e->get_num_points(o);
	
	rv1->set_active_element(el);
	rv1->set_quad_2d(quad_e);
	refmap_ex->set_quad_2d(quad_e);	

	double3* quad_points = rv1->get_quad_2d()->get_points(o);  //Quadraturpunkte bzgl e!
    double* x_phys = rv1->get_phys_x(o);   //physikalische Koordinaten von e
	double* y_phys = rv1->get_phys_y(o);
	double x, y;
	double3* new_quad_points = new double3[np];

	for(int i = 0; i<np; i++){ //Transformation von physikalischen Koordinaten von e, auf Referenzelement von elem_ex
		refmap_ex->untransform(elem_ex,x_phys[i],y_phys[i], x, y ); 
		new_quad_points[i][0]= x;
		new_quad_points[i][1] = y; 
		new_quad_points[i][2] = quad_points[i][2];
	}
		
	double* jac = NULL;

	Geom<double>* geo_e = init_geom_vol(rv1, o);		
		if(!rv1->is_jacobian_const()) {
			jac = rv1->get_jacobian(o);				 
		}
		double* jwt = new double[np];
		for(int i = 0; i < np; i++) {
		  if(rv1->is_jacobian_const())
		    jwt[i] = quad_points[i][2] * rv1->get_const_jacobian();
		  else
		   jwt[i] = quad_points[i][2] * jac[i];
		}	

	Func<scalar>* err1 = new Func<double>(np, nc);	
	err1->val = new scalar[np];
	Func<scalar>* err2 = new Func<double>(np, nc);	
	err2->val = new scalar[np];

	rsln1->set_active_element(el);
	rsln2->set_active_element(el);
	Func<scalar>* v1 = init_fn(rsln1, o);
	Func<scalar>* v2 = init_fn(rsln2, o);

	double* u = new scalar[np];
	for(int k = 0; k <np; k++)	u[k]= 0;

    for (unsigned int i = 0; i < al_ex->cnt; i++) {
     if (al_ex->dof[i] >= 0) {
		 int dof = al_ex->dof[i];
        scalar coef = al_ex->coef[i] * coeffs[dof];
		for(int k = 0; k <np; k++)					
				u[k] += pss->get_shapeset()->get_fn_value(al_ex->idx[i], new_quad_points[k][0],new_quad_points[k][1],0)*coef;				
			       	
      }
		   
    }
	for(int k = 0; k <np; k++){				
		err1->val[k] = u[k];
		err2->val[k] = u[k];
		
	}  
	  err1->subtract(*v1);
	  err2->subtract(*v2);
		
	 scalar res = form->value(np, jwt, NULL, err1, err2, geo_e, NULL);

	 double err = std::abs(res);

	if(solutions_for_adapt)
            this->errors[nc][el->id] = err;
	total_error += err;	
	

	delete [] new_quad_points;
	delete [] jwt;
	delete [] u;	
	  geo_e->free(); delete geo_e;
	  err1->free_fn(); delete err1;
	  err2->free_fn(); delete err2;
	  v1->free_fn(); delete v1;
	  v2->free_fn(); delete v2;
  }

	if(solutions_for_adapt) {
   	 have_errors = true;
 	 }

	delete refmap_ex;
	delete pss;
	delete al_ex;
	this->errors_squared_sum = total_error;
	return sqrt(total_error);
}













