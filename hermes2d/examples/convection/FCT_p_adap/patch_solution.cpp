#include "function/solution.h"



class HERMES_API PatchSolution : public Solution
{
 public:
		PatchSolution(Mesh* mesh): Solution(mesh) {
		};  


	void init(Space* patch_space, Space* space, Solution* sln,scalar* coeffs, int* list,  int elems);

	void init(Space* patch_space, Space* space, Solution* sln, int* list,  int elems, MatrixSolverType solver_type = SOLVER_UMFPACK){
		if ((space == NULL)||(patch_space == NULL)) error("Space == NULL in PatchSolution::init(Space* patch_space, Space* space, Solution* sln).");
		int ndof = space->get_num_dofs();
		scalar* coeffs= new scalar[ndof];
		OGProjection::project_global(space, sln, coeffs, solver_type, HERMES_L2_NORM);
		this->init( patch_space, space,  sln, coeffs,  list,  elems);
		delete [] coeffs;

	};

};



void PatchSolution::init(Space* patch_space, Space* space, Solution* sln,scalar* coeffs,  int* list, int elems){

 // sanity check
    if ((space == NULL)||(patch_space == NULL)) error("Space == NULL in PatchSolution::init(Space* patch_space, Space* space, Solution* sln).");
	  if ((space->get_mesh() == NULL)||  (patch_space->get_mesh() == NULL))  error("Mesh == NULL in PatchSolution::init(Space* patch_space, Space* space, Solution* sln).");
	  if (!space->is_up_to_date())   error("Provided 'space' is not up to date.");

	if((mesh!=NULL)&&(patch_space->get_mesh()!=mesh))error("Mesh !=patch_mesh in PatchSolution::init(Space* patch_space, Space* space, Solution* sln).");

    // initialize precalc shapeset using the space's shapeset
    Shapeset *shapeset = space->get_shapeset();
    if (space->get_shapeset() == NULL) error("Space->shapeset == NULL in PatchSolution::init(Space* patch_space, Space* space, Solution* sln)");
    PrecalcShapeset *pss = new PrecalcShapeset(shapeset);
    if (pss == NULL) error("PrecalcShapeset could not be allocated in PatchSolution::init(Space* patch_space, Space* space, Solution* sln).");

	
  free(); 
  space_type = space->get_type();
  num_components = pss->get_num_components();
  sln_type = HERMES_SLN;
		 int o;

  // copy the mesh  
	 Mesh* big_mesh = space->get_mesh();
 if(mesh==NULL) mesh = patch_space->get_mesh();

  // allocate the coefficient arrays
  num_elems = mesh->get_max_element_id();
	if(elems!=num_elems) 
		error("num_elems != elems in PatchSolution::init(Space* patch_space, Space* space, Solution* sln,scalar* coeffs,  int* list, int elems)");
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

 for(int j = 0; j <elems; j++){  //Alle Elemente durchlaufen
			e = space->get_mesh()->get_element(list[j]);			
    		mode = e->get_mode();
   		   o = space->get_element_order(e->id);
		   patch_space->set_element_order_internal(j, o);	//patch_space_init auf richtige Ordnung setzen.	
           o = std::max(H2D_GET_H_ORDER(o), H2D_GET_V_ORDER(o));
		 for (unsigned int k = 0; k < e->nvert; k++) {
			int eo = space->get_edge_order(e, k);
			 if (eo > o) o = eo;
		 }
           num_coefs += mode ? sqr(o+1) : (o+1)*(o+2)/2;
   		   elem_orders[j] = o;
  }
	patch_space->assign_dofs();


  num_coefs *= num_components;
  if(mono_coefs != NULL)
    delete [] mono_coefs;
  mono_coefs = new scalar[num_coefs];

  // express the solution on elements as a linear combination of monomials
  Quad2D* quad = &g_quad_2d_cheb;
  pss->set_quad_2d(quad);
  scalar* mono = mono_coefs;
 for(int j = 0; j <elems; j++){  //Alle Nachbarelemente durchlaufen
	e = space->get_mesh()->get_element(list[j]);
    mode = e->get_mode();
    quad->set_mode(mode);
    o = elem_orders[j];
    int np = quad->get_num_points(o);

    AsmList al;
    space->get_element_assembly_list(e, &al);
    pss->set_active_element(e);

    for (int l = 0; l < num_components; l++)
    {
      // obtain solution values for the current element
      scalar* val = mono;
      elem_coefs[l][j] = (int) (mono - mono_coefs);
      memset(val, 0, sizeof(scalar)*np);
      for (unsigned int k = 0; k < al.cnt; k++)
      {
        pss->set_active_shape(al.idx[k]);
        pss->set_quad_order(o, H2D_FN_VAL);
        int dof = al.dof[k];        
        scalar coef = al.coef[k] *  coeffs[dof] ;
        double* shape = pss->get_fn_values(l);
        for (int i = 0; i < np; i++)
          val[i] += shape[i] * coef;
      }
      mono += np;

      // solve for the monomial coefficients
      if (mono_lu.mat[mode][o] == NULL)
        mono_lu.mat[mode][o] = calc_mono_matrix(o, mono_lu.perm[mode][o]);
      lubksb(mono_lu.mat[mode][o], np, mono_lu.perm[mode][o], val);
    }
  }

  if(mesh == NULL) error("mesh == NULL.\n");
  init_dxdy_buffer();
  element = NULL;

 delete pss;

}
