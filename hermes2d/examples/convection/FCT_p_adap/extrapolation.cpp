#include "../src/function/solution.h"
#include "../src/function/solution.cpp"
#include "../src/mesh/refmap.h"
#include "../src/h2d_common.h"

/*

class HERMES_API RefMap_Extrapolation : public RefMap
{
 public:
  	RefMap_Extrapolation(Element* e) : RefMap() {
		set_active_element(e);			
		};
 	 ~RefMap_Extrapolation() { };	
	 void jac_at_point(double xi1, double xi2, double& x, double& y, double jac, Shapeset* shapeset);
};
//Berechnung Jacobi/Funktionaldeterminante an Referenzpunkt (xi1,xi2)=>jac, (x,y) physikalische Koord.
void RefMap_Extrapolation::jac_at_point(double xi1, double xi2, double& x, double& y, double jac, Shapeset* shapeset)
{
  double2x2 tmp;
  memset(tmp, 0, sizeof(double2x2));
  x = y = 0;
  for (int i = 0; i < nc; i++)
  {
    double val = shapeset->get_fn_value(indices[i], xi1, xi2, 0);
   x += coeffs[i][0] * val;
    y += coeffs[i][1] * val;

    double dx =  shapeset->get_dx_value(indices[i], xi1, xi2, 0);
    double dy =  shapeset->get_dy_value(indices[i], xi1, xi2, 0);
    tmp[0][0] += coeffs[i][0] * dx;
    tmp[0][1] += coeffs[i][0] * dy;
    tmp[1][0] += coeffs[i][1] * dx;
    tmp[1][1] += coeffs[i][1] * dy;
  }  
  jac = tmp[0][0] * tmp[1][1] - tmp[0][1] * tmp[1][0];

}

*/


class HERMES_API ExtrapolationSolution : public Solution
{
 public:
	  ExtrapolationSolution(Mesh* mesh, Element* e): Solution(mesh) {
  		//this->refmap = new RefMap_Extrapolation(e);	
		refmap_ex= new RefMap();		
		elem_ex = e;
		coeff_basis = NULL;		
		};   
	  ~ExtrapolationSolution(){
			if(coeff_basis!=NULL){ 
				delete [] coeff_basis; 
				coeff_basis=NULL;
			}
			elem_ex =NULL;
			u_h =NULL;
			if(refmap_ex!=NULL) delete refmap_ex;
			refmap_ex = NULL;
		};

	void set_elem_ex(Element* e){				
		elem_ex = e;
	};
	
	void set_fn(Solution* u){
		u_h = u;
	
	};	
	void init_extrapolation(Space* space);

	void init_zero(Space* space);

	scalar* get_coeff_basis(){
		return coeff_basis;
	};
	


protected:	
	Element* elem_ex;			//Element e bzgl dessen Basis extrapoliert wird
	Solution* u_h;   				//zu extrapolierende Funktion
	scalar* coeff_basis;
	RefMap* refmap_ex;

};


void ExtrapolationSolution::init_extrapolation(Space* space){   	
  // sanity check
    if (space == NULL) error("Space == NULL in ExtrapolationSolution::init_extrapolation(Space* space).");
	  if (space->get_mesh() == NULL) error("Mesh == NULL in ExtrapolationSolution::init_extrapolation(Space* space)");
	if(u_h == NULL) error("uh == NULL in ExtrapolationSolution:init_extrapolation(Space* space).");

    // initialize precalc shapeset using the space's shapeset
    Shapeset *shapeset = space->get_shapeset();
    if (space->get_shapeset() == NULL) error("Space->shapeset == NULL in ExtrapolationSolution::init_extrapolation(Space* space)");
    PrecalcShapeset *pss = new PrecalcShapeset(shapeset);
    if (pss == NULL) error("PrecalcShapeset could not be allocated in ExtrapolationSolution::init_extrapolation(Space* space).");	 

  free();

	int o;
 
  space_type = space->get_type();

  num_components = pss->get_num_components();
  sln_type = HERMES_SLN;
  num_dofs = space->get_num_dofs();

  //bool add_dir_lift = false;
  scalar* coeffs= new scalar[num_dofs];
  for(int i = 0; i< num_dofs; i++) coeffs[i]= 0.0;

	OGProjection::project_global(space,u_h, coeffs, matrix_solver, HERMES_L2_NORM);

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
      elem_coefs[l][e->id] = (int) (mono - mono_coefs);
      memset(val, 0, sizeof(scalar)*np);
      for (unsigned int k = 0; k < al_ex->cnt; k++)     
	 {
        pss->set_active_shape(al_ex->idx[k]);  //shapeset von elem_ex!!
        pss->set_quad_order(o, H2D_FN_VAL);    
        int dof = al_ex->dof[k];
        //double dir_lift_coeff = add_dir_lift ? 1.0 : 0.0;
        //scalar coef = al_ex->coef[k] * (dof >= 0 ? coeffs[dof] : dir_lift_coeff);
		scalar coef = al_ex->coef[k] * coeffs[dof] ;
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


	delete [] coeffs;
	delete pss;
	delete al_ex;

	
}



void ExtrapolationSolution::init_zero(Space* space){  //erstmal mit 0 initialisieren
 
  // sanity check
    if (space == NULL) error("Space == NULL in ExtrapolationSolution::init_zero(Space* space).");

    // initialize precalc shapeset using the space's shapeset
    Shapeset *shapeset = space->get_shapeset();
    if (space->get_shapeset() == NULL) error("Space->shapeset == NULL in ExtrapolationSolution::init_zero(Space* space)");
    PrecalcShapeset *pss = new PrecalcShapeset(shapeset);
    if (pss == NULL) error("PrecalcShapeset could not be allocated in ExtrapolationSolution::init_zero(Space* space).");

	 int o;

  free();
 
  space_type = space->get_type();

  num_components = pss->get_num_components();
  sln_type = HERMES_SLN;
  num_dofs = space->get_num_dofs();

  bool add_dir_lift = false;
  scalar* coeffs= new scalar[num_dofs];
  for(int i = 0; i< num_dofs; i++) coeffs[i]= 0.0;

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
	
	AsmList al_ex;
	space->get_element_assembly_list(elem_ex, &al_ex);    
	pss->set_active_element(elem_ex);
	mode = elem_ex->get_mode();
    quad->set_mode(mode);

	//coeff_basis setzen (Speicherung auch als lin Kombi von basisfunktionen!)
	if(coeff_basis!=NULL) delete [] coeff_basis;
	coeff_basis = new scalar[al_ex.cnt];
	for (unsigned int k = 0; k < al_ex.cnt; k++)   coeff_basis[k]= coeffs[k]; 

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
    double* x_phys = refmap->get_phys_x(o);
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
      elem_coefs[l][e->id] = (int) (mono - mono_coefs);
      memset(val, 0, sizeof(scalar)*np);
      for (unsigned int k = 0; k < al_ex.cnt; k++)     
	 {
        pss->set_active_shape(al_ex.idx[k]);  //shapeset von elem_ex!!
        pss->set_quad_order(o, H2D_FN_VAL);    
        int dof = al_ex.dof[k];
        double dir_lift_coeff = add_dir_lift ? 1.0 : 0.0;
        scalar coef = al_ex.coef[k] * (dof >= 0 ? coeffs[dof] : dir_lift_coeff);
		for (int i = 0; i < np; i++){
        	double shape = pss->get_shapeset()->get_fn_value(al_ex.idx[k], new_quad_points[i][0],new_quad_points[i][1],0); 
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

	delete [] coeffs;
	delete pss;
}











