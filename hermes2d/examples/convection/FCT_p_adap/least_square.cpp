#include "../hermes_common/common.h"
#include "function/solution.h"
#include "function/forms.h"
#include "hermes2d.h"

#include "weakform/weakform.h"
#include "integrals/h1.h"




using namespace WeakFormsH1;


class  CustomWeakFormLS  : public WeakForm
{
public:
  CustomWeakFormLS(Solution* sln) : WeakForm(1) {
    add_matrix_form(new CustomWeakFormLSVol(0, 0));
    VectorFormVolLS* vector_form = new VectorFormVolLS(0);
    vector_form->ext.push_back(sln);
    add_vector_form(vector_form);
  };

private:
  class CustomWeakFormLSVol : public WeakForm::MatrixFormVol
  {
  public:
   
    CustomWeakFormLSVol(int i, int j) 
      : WeakForm::MatrixFormVol(i, j, HERMES_ANY, HERMES_NONSYM) { }

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, 
                       Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const {
		if((u->val==NULL)||(v->val==NULL)||(wt==NULL)) error("CustomWeakFormLSVol- NULL");
		if((u->num_gip != n)||(v->num_gip != n))error("CustomWeakFormLSVol-n");
         Scalar result = 0; 			
		  for (int i = 0; i < n; i++)				
			result += wt[i] * (u->val[i] * v->val[i]);		
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



  };


  class VectorFormVolLS : public WeakForm::VectorFormVol
  {
  public:
    VectorFormVolLS(int i) : WeakForm::VectorFormVol(i) { }

    template<typename Real, typename Scalar>
    Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const {
		  Scalar result = 0;
		  Func<Scalar>* sln = ext->fn[0];
		 if(sln ==NULL) error("VectorFormVolLS: sln ==NULL");
		  for (int i = 0; i < n; i++)
			result += wt[i] *  sln->val[i] * v->val[i];
		  return result;

    }

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const {
      return vector_form<double, scalar>(n, wt, u_ext, v, e, ext);
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const {
      return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
    }


  
  };
};



class HERMES_API DiscreteProblem_LS : public DiscreteProblem
{



public:
  /// Constructors. 
  DiscreteProblem_LS(WeakForm* wf, Space* space, Element* e): DiscreteProblem(wf, space){
			elem_ex = e;
			//refmap_ex = new RefMap_Extrapolation(e);		
			refmap_ex = new RefMap;	
			pss_ex_i = new PrecalcShapeset(space->get_shapeset());
			pss_ex_j = new PrecalcShapeset(space->get_shapeset());
			al_ex = new AsmList();
			space->get_element_assembly_list(elem_ex, al_ex); 
			 };

  /// Destuctor.
 virtual ~DiscreteProblem_LS(){	
		if(refmap_ex!=NULL){ delete refmap_ex; refmap_ex=NULL;}
		if(pss_ex_i!=NULL){delete pss_ex_i; pss_ex_i=NULL;}
		if(pss_ex_j!=NULL){delete pss_ex_j; pss_ex_j =NULL;}
		if(al_ex!=NULL){delete al_ex; al_ex =NULL;}
	};

  void create_sparse_structure(SparseMatrix* mat, Vector* rhs = NULL);
	void assemble_all(SparseMatrix* mat, WeakForm::MatrixFormVol *mfv, Vector* rhs, WeakForm::VectorFormVol *vfv);

	



protected:	
	Element* elem_ex;	//Element e bzgl dessen Basis extrapoliert wird
	//RefMap_Extrapolation* refmap_ex;   //RefMap bzgl elem_ex;
	RefMap* refmap_ex;   //RefMap bzgl elem_ex;
	PrecalcShapeset* pss_ex_i ;       //PrecalcShapeset  bzgl elem_ex
	PrecalcShapeset* pss_ex_j ;       //PrecalcShapeset bzgl elem_ex
	AsmList* al_ex;					//AsmList von elem_ex
	   

};


//Die Freiheitsgrade befinden sich lediglich auf Element e! Erst bei der Berechnung der Integrale
// kommen andere Elemente hinzu!
void DiscreteProblem_LS::create_sparse_structure(SparseMatrix* mat, Vector* rhs)
{
  _F_

  if (is_up_to_date())
  {
    if (mat != NULL)
    {
      verbose("Reusing matrix sparse structure.");
      mat->zero();
    }
    if (rhs != NULL) rhs->zero();
    return;
  }
  
  int ndof = al_ex->cnt;  //Achtung: die Anzahl der Freiheitsgrade entspricht nur noch der Anzahl der Freiheitsgrade von elem_ex!!!

  if (mat != NULL)  
  {
    // Spaces have changed: create the matrix from scratch.
    have_matrix = true;
    mat->free();
    mat->prealloc(ndof);
	if(wf->get_neq()!=1.0) error("too many equations in DiscreteProblem_LS::create_sparse_structure");

    AsmList* am = al_ex;
        AsmList* an = al_ex;

        // Pretend assembling of the element stiffness matrix.
        for (unsigned int i = 0; i < am->cnt; i++) {
          if (am->dof[i] >= 0) {
            for (unsigned int j = 0; j < an->cnt; j++) {
              if (an->dof[j] >= 0) {
                mat->pre_add_ij(am->dof[i], an->dof[j]);
              }
            }
          }
        }
    mat->alloc();
  }

  // WARNING: unlike Matrix::alloc(), Vector::alloc(ndof) frees the memory occupied
  // by previous vector before allocating
  if (rhs != NULL) rhs->alloc(ndof);

  // save space seq numbers and weakform seq number, so we can detect their changes
  for (unsigned int i = 0; i < wf->get_neq(); i++) sp_seq[i] = spaces[i]->get_seq();

  wf_seq = wf->get_seq();

  struct_changed = true;
}

void DiscreteProblem_LS::assemble_all(SparseMatrix* mat, WeakForm::MatrixFormVol *mfv, Vector* rhs, WeakForm::VectorFormVol *vfv){
	Space* space = spaces[0];
	if(wf->get_neq()>1) error(" DiscreteProblem_LS::assemble_mat: kein multi-mesh moeglich!!!");
	if(mfv ==NULL) error(" DiscreteProblem_LS::assemble_mat: kein MatrixFormVol!!!");
	Mesh* mesh = space->get_mesh();	  
	init_cache();

	Quad2D* quad = pss_ex_i->get_quad_2d(); 
	pss_ex_i->set_active_element(elem_ex);
	pss_ex_j->set_active_element(elem_ex);	
    quad->set_mode(elem_ex->get_mode());
	refmap_ex->set_active_element(elem_ex);
	 refmap_ex->set_quad_2d(quad);

	Element* el;
	RefMap* refmap = new RefMap();
	PrecalcShapeset* pss = new PrecalcShapeset(space->get_shapeset());
	Quad2D* quad_e = pss->get_quad_2d(); 
	AsmList* ai = al_ex;
     AsmList* aj = al_ex;
	scalar **local_stiffness_matrix = NULL;

  for_all_active_elements(el, mesh) 
 {	    
    quad_e->set_mode(el->get_mode());
   	int o = space->get_element_order(el->id);
	 o = std::max(H2D_GET_H_ORDER(o), H2D_GET_V_ORDER(o));
    int np = quad_e->get_num_points(o);
  // printf("\n el->id= %i, np = %i, o = %i\n",el->id, np, o);
	
	refmap->set_active_element(el);
	refmap->set_quad_2d(quad_e);

	double3* quad_points = refmap->get_quad_2d()->get_points(o);  //Quadraturpunkte bzgl e!
    double* x_phys = refmap->get_phys_x(o);   //physikalische Koordinaten von e
	double* y_phys = refmap->get_phys_y(o);
	double x, y;
	double3* new_quad_points = new double3[np];

	for(int i = 0; i<np; i++){ //Transformation von physikalischen Koordinaten von e, auf Referenzelement von elem_ex
		refmap_ex->untransform(elem_ex,x_phys[i],y_phys[i], x, y ); 
		new_quad_points[i][0]= x;
		new_quad_points[i][1] = y; 
		new_quad_points[i][2] = quad_points[i][2];
	}

	/*if(el!=elem_ex){
		for(int i = 0; i<np; i++){ 
			printf(" x_ref= %f, y_ref= %f \n",quad_points[i][0], quad_points[i][1]);
			printf(" x_phys= %f, y_phys= %f \n",x_phys[i], y_phys[i]);
			printf(" x_ref_ex= %f, y_ref_ex= %f \n",new_quad_points[i][0], new_quad_points[i][1] );
			printf(" --------------------------------------------" );
		}
	}*/
		
	double* jac = NULL;
		  // Init geometry and jacobian*weights.	
	 if (cache_e[o] == NULL){
		cache_e[o] = init_geom_vol(refmap, o);		
		if(!refmap->is_jacobian_const()) {
			jac = refmap->get_jacobian(o);		 
		}
		cache_jwt[o] = new double[np];
		for(int i = 0; i < np; i++) {
		  if(refmap->is_jacobian_const())
		    cache_jwt[o][i] = quad_points[i][2] * refmap->get_const_jacobian();
		  else
		    cache_jwt[o][i] = quad_points[i][2] * jac[i];
		}
	 } else error("cache NULL sein");
	  Geom<double>* geo_e = cache_e[o];
	  double* jwt = cache_jwt[o];
	if(jwt==NULL) error("jwt =NULL");

	int nc = pss_ex_i->get_num_components();

	  Func<double>* u = new Func<double>(np, nc);
	 Func<double>* v = new Func<double>(np, nc);
	u->val = new double[np];
	v->val = new double[np];

	vfv->ext[0]->set_active_element(el);
	ExtData<scalar>* ext = init_ext_fns(vfv->ext, refmap, o);   //externe Funktion (u_h) setzen

	/*printf("ext: \n");
	Func<scalar>* sln = ext->fn[0];
	for (int i = 0; i < np; i++) printf(" x_phys= %f, y_phys= %f, u(x,y) = %f \n",x_phys[i], y_phys[i], sln->val[i]);
		printf(" \n");*/	
	

	scalar* result = new scalar[ai->cnt];

  	local_stiffness_matrix = get_matrix_buffer(al_ex->cnt);

    for (unsigned int i = 0; i < ai->cnt; i++) {
	//printf("\n ai_coef(idx= %i): %f, dof:%i, \n",ai->idx[i], ai->coef[i],ai->dof[i]);
     if (ai->dof[i] >= 0) {
		for(int k = 0; k <np; k++)				
			u->val[k] = pss_ex_i->get_shapeset()->get_fn_value(ai->idx[i], new_quad_points[k][0],new_quad_points[k][1],0)*ai->coef[i];	
	/*printf("\n u: \n");
	for (int l = 0; l < np; l++) printf(" %f, ", u->val[l]);
		printf(" \n");*/

		result[ai->dof[i]] = vfv->value(np, jwt, NULL, u, geo_e, ext);		

        for (unsigned int j = 0; j < aj->cnt; j++) {
          if (aj->dof[j] >= 0) {		
			for(int k = 0; k <np; k++)
				v->val[k] = pss_ex_j->get_shapeset()->get_fn_value(aj->idx[j], new_quad_points[k][0],new_quad_points[k][1],0)*aj->coef[j];
	     	
  			// The actual calculation takes place here.
 			 scalar res = mfv->value(np, jwt, NULL, u, v, geo_e, NULL) ; 		
			local_stiffness_matrix[i][j] = res;
        	
          }
		   
        }
      }
    }
	/*if(el==elem_ex) printf(" Result(e): \n");
	else	printf(" Result: \n");
	for (unsigned int i = 0; i < ai->cnt; i++) printf(" %f, ", result[i]);
*/

    if (mat != NULL)  mat->add(al_ex->cnt, al_ex->cnt, local_stiffness_matrix, al_ex->dof, al_ex->dof);  //liste mit dof->ordnet eintraege richtig zu!!
	if(rhs!=NULL)  rhs->add_vector(result);


//Clean up
	delete_single_geom_cache(o);
	if( new_quad_points !=NULL)
	{ 	delete [] new_quad_points;
		 new_quad_points =NULL;}
	if (u != NULL) {
		u->free_fn();	
		delete u;
	}
	if (v != NULL) {
		v->free_fn();	
		delete v;
	}
	if(result!=NULL){ delete [] result;  result= NULL;}
	if (ext != NULL) {
    	ext->free(); 
   	 	delete ext;
  	}


  }

/*printf("Matrix A: \n");
	 for (unsigned int i = 0; i < ai->cnt; i++) {
		for (unsigned int j = 0; j < aj->cnt; j++) {
				printf("%f ", mat->get(i,j));
				//if(j<aj->cnt-1) printf(", ");
				//else printf(" ; ");
		} //printf("= %f \n", rhs->get(i));
		printf(" ; \n");
	}
	printf(" \n  = Vektor: \n");
	for (unsigned int i = 0; i < ai->cnt; i++) printf(" %f; ", rhs->get(i));
	printf(" \n  \n");
*/

	delete refmap;
	delete pss;

	// Deinitialize matrix buffer.
	if(matrix_buffer != NULL)
	delete [] matrix_buffer;
	matrix_buffer = NULL;
	matrix_buffer_dim = 0;




}

class HERMES_API LSSolution : public Solution
{
 public:
	  LSSolution(Mesh* mesh, Element* e): Solution(mesh) {
  		//this->refmap = new RefMap_Extrapolation(e);			
		elem_ex = e;
		refmap_ex = new RefMap();
		coeff_basis = NULL;		
		}; 

  
	  ~LSSolution(){ 
			if(coeff_basis!=NULL){ 
				delete [] coeff_basis; 
				coeff_basis =NULL;
			}
			if(refmap_ex!=NULL) delete refmap_ex;
			refmap_ex = NULL;
			elem_ex =NULL;
	};

	void set_elem_ex(Element* e){				
		elem_ex = e;
	};

	void LS_project(Solution* u_h, Space* space, MatrixSolverType matrix_solver = SOLVER_UMFPACK);
	void init(Space* space, scalar* coeffs);  
	
	scalar* get_coeff_basis(){
		return coeff_basis;
	};
	

protected:	
	Element* elem_ex;			//Element e bzgl dessen Basis 
	scalar* coeff_basis;
	RefMap* refmap_ex;


};


void LSSolution::init(Space* space, scalar* coeffs){

 // sanity check
    if (space == NULL) error("Space == NULL in LSSolution::init(Space* space, scalar* coeffs).");
	
    // initialize precalc shapeset using the space's shapeset
    Shapeset *shapeset = space->get_shapeset();
    if (space->get_shapeset() == NULL) error("Space->shapeset == NULL in LSSolution::init(Space* space, scalar* coeffs)");
    PrecalcShapeset *pss = new PrecalcShapeset(shapeset);
    if (pss == NULL) error("PrecalcShapeset could not be allocated in LSSolution::init(Space* space, scalar* coeffs).");
	 int o;
	
  free(); 
  space_type = space->get_type();
  num_components = pss->get_num_components();
  sln_type = HERMES_SLN;

  // copy the mesh  
  mesh = space->get_mesh();

  // allocate the coefficient arrays
  num_elems = mesh->get_max_element_id();
  if(elem_orders != NULL)
    delete [] elem_orders;
  elem_orders = new int[num_elems];
  memset(elem_orders, 0, sizeof(int) * num_elems);
  for (int l = 0; l < num_components; l++) {
    if(elem_coefs[l] != NULL)
      delete [] elem_coefs[l];
    elem_coefs[l] = new int[num_elems];
    memset(elem_coefs[l], 0, sizeof(int) * num_elems);
  }

  // obtain element orders, allocate mono_coefs
  Element* e;
  num_coefs = 0;
  for_all_active_elements(e, mesh)
  {
    mode = e->get_mode();
    o = space->get_element_order(e->id);
    o = std::max(H2D_GET_H_ORDER(o), H2D_GET_V_ORDER(o));
	for (unsigned int k = 0; k < e->nvert; k++) {
      int eo = space->get_edge_order(e, k);
      if (eo > o) o = eo;
    }

    num_coefs += mode ? sqr(o+1) : (o+1)*(o+2)/2;
    elem_orders[e->id] = o;
  }
  num_coefs *= num_components;
  if(mono_coefs != NULL)
    delete [] mono_coefs;
  mono_coefs = new scalar[num_coefs];

  // express the solution on elements as a linear combination of monomials
  Quad2D* quad = &g_quad_2d_cheb;
  pss->set_quad_2d(quad);
  scalar* mono = mono_coefs;
	
	AsmList* al_ex = new AsmList;
	space->get_element_assembly_list(elem_ex, al_ex);    
	pss->set_active_element(elem_ex);
	mode = elem_ex->get_mode();
    quad->set_mode(mode);


	//coeff_basis setzen (Speicherung auch als lin Kombi von basisfunktionen!)
	if(coeff_basis!=NULL) delete [] coeff_basis;
	coeff_basis = new scalar[al_ex->cnt];
	for (unsigned int k = 0; k < al_ex->cnt; k++)   coeff_basis[k]= coeffs[k]; 

/*printf("coeffs:\n ");
 for (unsigned int k = 0; k < al_ex->cnt; k++)   printf("%f ,",coeffs[k]);*/


	Quad2D* quad_e = &g_quad_2d_cheb;

  for_all_active_elements(e, mesh)
  {
    mode = e->get_mode();
    quad_e->set_mode(mode);
    o = elem_orders[e->id];
    int np = quad->get_num_points(o);
	
	refmap->set_active_element(e);
	refmap->set_quad_2d(quad_e);
	double3* quad_points = refmap->get_quad_2d()->get_points(o);  //Quadraturpunkte bzgl e!
    double* x_phys = refmap->get_phys_x(o);   //physikalische Koordinaten von e
	double* y_phys = refmap->get_phys_y(o);
	double x, y;
	double3* new_quad_points = new double3[np];

	refmap_ex->set_active_element(elem_ex);
	 refmap_ex->set_quad_2d(quad);
	for(int i = 0; i<np; i++){ //Transformation von physikalischen Koordinaten von e, auf Referenzelement von elem_ex		
		refmap_ex->untransform(elem_ex,x_phys[i],y_phys[i], x, y ); 
		new_quad_points[i][0]= x;
		new_quad_points[i][1] = y; 
		new_quad_points[i][2] = quad_points[i][2];
	}	

    for (int l = 0; l < num_components; l++)
    {
      // obtain solution values for the current element
      scalar* val = mono;
      elem_coefs[l][e->id] = (int) (mono - mono_coefs);  //zeigt auf Eintrag in mono_coeffs
      memset(val, 0, sizeof(scalar)*np);
      for (unsigned int k = 0; k < al_ex->cnt; k++)     
	 {
        pss->set_active_shape(al_ex->idx[k]);  //shapeset von elem_ex!!
        pss->set_quad_order(o, H2D_FN_VAL);    
        int dof = al_ex->dof[k];    
        scalar coef = al_ex->coef[k] * coeffs[dof];
		for (int i = 0; i < np; i++){
        	double shape = pss->get_shapeset()->get_fn_value(al_ex->idx[k], new_quad_points[i][0],new_quad_points[i][1],0); 
																			//Wert an Quadpoint in refDomain       
          	val[i] += shape * coef;   //Werte an verschiedenen Quadraturpunkten
		}
      }
      mono += np;

      // solve for the monomial coefficients
      if (mono_lu.mat[mode][o] == NULL)
        mono_lu.mat[mode][o] = calc_mono_matrix(o, mono_lu.perm[mode][o]);
      lubksb(mono_lu.mat[mode][o], np, mono_lu.perm[mode][o], val);
    }
	delete [] new_quad_points;
  }

  if(mesh == NULL) error("mesh == NULL.\n");
  init_dxdy_buffer();
  element = NULL;
  delete al_ex;
  delete pss;

	
}



  

void LSSolution::LS_project(Solution* u_h, Space* space, MatrixSolverType matrix_solver){

	if(elem_ex==NULL)error("LSSolution::LS_project(): elem_ex =NULL");
	CustomWeakFormLS* wf = new CustomWeakFormLS(u_h); 
	 DiscreteProblem_LS* dp = new DiscreteProblem_LS(wf, space, elem_ex);
	SparseMatrix* mat = create_matrix(matrix_solver);  
	Vector* rhs = create_vector(matrix_solver); 	

    dp->create_sparse_structure(mat,rhs);	
	dp->assemble_all(mat, wf->get_mfvol()[0], rhs, wf->get_vfvol()[0]);

	Solver* solver = create_linear_solver(matrix_solver, mat, rhs);
	scalar* coeff_vec =NULL;
	if(solver->solve()){ 
			coeff_vec = solver->get_solution();		
		}else error("solver doesn't work  in LS_project(Solution* u_h, Space* space, MatrixSolverType matrix_solver)");
	
	/*int ndof = mat->get_size();
	for(int i = 0; i<ndof; i++) printf("LS-coeff: %f, ", coeff_vec[i]);
	printf("\n");*/

/*
//neu: Newton-Verfahren
	//rhs(i) = u_h*phi_i 
	//residual = rhs_const + jac*Y^n
	//jac(i,j)  = phi_i * phi_j
	//Jac(Y^n) \deltaY^{n+1} = -F(Y^n) (=rhs).

	int ndof = mat->get_size();
	UMFPackVector* residual = new UMFPackVector(ndof);	
	residual->zero();
	Solver* solver_newton = create_linear_solver(matrix_solver, mat, residual);
	scalar* sln_newton = new scalar[ndof];
	for(int i = 0;i<ndof;i++) sln_newton[i] = coeff_vec[i];
	
	double max_allowed_residual_norm = 1e6;
	double damping_coeff = 1.0; 
	Hermes2D hermes2d;

	scalar* rhs_var = new scalar[ndof];
	for(int i=0; i<ndof;i++) rhs_var[i]=0.0;

	  // The Newton's loop.
	  double residual_norm;
	  int it = 1;
	  while (1)
	  {	
		residual->zero(); residual->add_vector(rhs);
	    	mat->multiply_with_vector(sln_newton, rhs_var);
		for(int i = 0; i<ndof;i++) rhs_var[i] *= (-1.0);	
		residual->add_vector(rhs_var); // residual= -F(Y^n);
	   
	      // Calculate the l2-norm of residual vector
	      //residual_norm = Hermes2D::get_l2_norm(rhs);
	   residual_norm = hermes2d.get_l2_norm(residual);

	 // If maximum allowed residual norm is exceeded, fail.
	    if (residual_norm > max_allowed_residual_norm) {	 
		info("Current residual norm: %g", residual_norm);
		info("Maximum allowed residual norm: %g", max_allowed_residual_norm);
		info("Newton solve not successful, returning false.");	      		 
	    }

	    // If residual norm is within tolerance, or the maximum number
	    // of iteration has been reached, then quit.
	    if ((residual_norm < NEWTON_TOL || it > NEWTON_MAX_ITER) && it > 1) break; 
	
	    // Solve the linear system.
	    if(!solver_newton->solve()) error ("Matrix solver failed.\n");

	    // Add \deltaY^{n+1} to Y^n.
		for (int i = 0; i < ndof; i++) sln_newton[i] += damping_coeff * solver_newton->get_solution()[i];

	    it++;
	  }
	  if (it >= NEWTON_MAX_ITER) {
	    info("Maximum allowed number of Newton iterations (%i) exceeded, returning false.", it);
		info("Current residual norm: %g", residual_norm);	   
	  }
	delete [] rhs_var;
	this->init(space,sln_newton);
	delete solver_newton;
	delete [] sln_newton;
*/

	this->init(space,coeff_vec);



	wf->delete_all();
	delete wf;
	delete dp;
	delete mat;
	delete rhs;
	delete solver;
	
	

}
