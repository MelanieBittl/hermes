#include "lumped_projection.h"



void Lumped_Projection::project_internal(Hermes::vector<Space *> spaces, WeakForm *wf, scalar* target_vec,
                               MatrixSolverType matrix_solver,UMFPackMatrix*  mat ){

  unsigned int n = spaces.size();

  // sanity checks
  if (n <= 0 || n > 10) error("Wrong number of projected functions in project_internal().");
  for (unsigned int i = 0; i < n; i++) if(spaces[i] == NULL) error("this->spaces[%d] == NULL in project_internal().", i);
  if (spaces.size() != n) error("Number of spaces must match number of projected functions in project_internal().");

  // this is needed since spaces may have their DOFs enumerated only locally.
  int ndof = Space::assign_dofs(spaces);

  // Initialize DiscreteProblem.
  DiscreteProblem* dp = new DiscreteProblem(wf, spaces);
    	UMFPackMatrix* matrix = new UMFPackMatrix;	
  	UMFPackVector* rhs = new UMFPackVector(ndof);
	scalar* coeff_vec =NULL; 
	if(mat==NULL) { 
		UMFPackMatrix* lumped_matrix = new UMFPackMatrix;   //M_L 
		dp->assemble(matrix, rhs);  
			//Masslumping		 
		 int size = matrix->get_size();
		 scalar diag[size];
		 int nnz = matrix->get_nnz();
		 int row[size]; 
		int col[size+1];
		 for(int i = 0; i<size; i++){    
		    diag[i] = 0;
		    row[i]= i;
		    col[i]=i;
		 }
		col[size]=size;// letzter Eintrag bezieht sich auf nichts(Ende der Matrix)!
	
		 for(int i = 0; i<nnz; i++){    
		    diag[matrix->get_Ai()[i]] += matrix->get_Ax()[i]; 
		 }

		 lumped_matrix->create(size, size, col, row, diag);  //lumped Matrix aufstellen
		UMFPackLinearSolver* solver = new UMFPackLinearSolver(lumped_matrix,rhs);		
		if(solver->solve()){ 
			coeff_vec = solver->get_solution();		
		}
	  	else error ("Matrix solver failed.\n");
		 if (target_vec != NULL)
    		for (int i=0; i < ndof; i++) target_vec[i] = coeff_vec[i];
		delete solver;
		delete lumped_matrix;
	}else{ 
		dp->assemble(NULL, rhs);
		UMFPackLinearSolver* solver = new UMFPackLinearSolver(mat,rhs);		
		if(solver->solve()) 
			coeff_vec = solver->get_solution();			
	 	 else error ("Matrix solver failed.\n");
		 if (target_vec != NULL)
    		for (int i=0; i < ndof; i++) target_vec[i] = coeff_vec[i];
		delete solver;
	}
  
  
  delete matrix;
  delete rhs;
  delete dp;
  
  
}




void Lumped_Projection::project_lumped( Hermes::vector<Space *> spaces, Hermes::vector<MeshFunction *> source_meshfns,
                             scalar* target_vec, MatrixSolverType matrix_solver, UMFPackMatrix*  mat)
{/*
  _F_
   int n = spaces.size();
  // define temporary projection weak form
  WeakForm* proj_wf = new WeakForm(n);
   if(n>1){ 
	for (int i = 0; i < n; i++) {    
    		proj_wf->add_matrix_form(new ProjectionLumpedMatrixFormVol(i, i));    
    		proj_wf->add_vector_form(new ProjectionLumpedVectorFormVol(i, source_meshfns[i]));
  	}
	 project_internal( spaces, proj_wf, target_vec, matrix_solver,mat);
  }else{
 	int i = 0;
  	ProjectionLumpedMatrixFormVol* form = new ProjectionLumpedMatrixFormVol(i, i);
  	ProjectionLumpedVectorFormVol* vector = new ProjectionLumpedVectorFormVol(i, source_meshfns[i]);
	  proj_wf->add_matrix_form(form);    
    	proj_wf->add_vector_form(vector);
	 project_internal(spaces, proj_wf, target_vec, matrix_solver, mat);
	  delete form;
         delete vector;
  }
 
 
  delete proj_wf;*/
  _F_
  int n = spaces.size();

  // define temporary projection weak form
  WeakForm* proj_wf = new WeakForm(n);
  ProjectionLumpedMatrixFormVol** matrix_form = new ProjectionLumpedMatrixFormVol*[n];
  ProjectionLumpedVectorFormVol** vector_form = new ProjectionLumpedVectorFormVol*[n];
  int found[100];
  for (int i = 0; i < 100; i++) found[i] = 0;
  for (int i = 0; i < n; i++)
  {


    found[i] = 1;
    // Jacobian.
	matrix_form[i] = new ProjectionLumpedMatrixFormVol(i, i);
	proj_wf->add_matrix_form(matrix_form[i]);

    // Residual.
    vector_form[i] = new ProjectionLumpedVectorFormVol(i, source_meshfns[i]);
	proj_wf->add_vector_form(vector_form[i]);

  }
  for (int i=0; i < n; i++)
  {
    if (found[i] == 0)
    {
      warn("Index of component: %d\n", i);
      error("Wrong projection norm in project_global().");
    }
  }

  project_internal(spaces, proj_wf, target_vec, matrix_solver, mat);

  for (int i = 0; i < n; i++)
  {
	delete vector_form[i];
	delete matrix_form[i];
  }
	delete [] vector_form;
	delete [] matrix_form;
  
  delete proj_wf;

}






void Lumped_Projection::project_lumped(Hermes::vector<Space *> spaces, Hermes::vector<Solution *> source_sols,
                             scalar* target_vec, MatrixSolverType matrix_solver, UMFPackMatrix*  mat )
{
  Hermes::vector<MeshFunction *> mesh_fns;
  for(unsigned int i = 0; i < source_sols.size(); i++)
    mesh_fns.push_back(source_sols[i]);
  project_lumped(spaces, mesh_fns, target_vec, matrix_solver,mat);
}


void Lumped_Projection::project_lumped(Space* space, MeshFunction* source_meshfn,
                             scalar* target_vec, MatrixSolverType matrix_solver, UMFPackMatrix*  mat )
{
  Hermes::vector<Space *> spaces;
  spaces.push_back(space);
  Hermes::vector<MeshFunction *> source_meshfns;
  source_meshfns.push_back(source_meshfn);
  project_lumped(spaces, source_meshfns, target_vec, matrix_solver,  mat);
}


void Lumped_Projection::project_lumped(Hermes::vector<Space *> spaces,
                             Hermes::vector<Solution*> sols_src, Hermes::vector<Solution*> sols_dest,
                             MatrixSolverType matrix_solver, UMFPackMatrix*  mat )
{
  _F_

  scalar* target_vec = new scalar[Space::get_num_dofs(spaces)];
  Hermes::vector<MeshFunction *> ref_slns_mf;
  for (unsigned int i = 0; i < sols_src.size(); i++)
    ref_slns_mf.push_back(static_cast<MeshFunction*>(sols_src[i]));

  project_lumped(spaces, ref_slns_mf, target_vec, matrix_solver);
  Solution::vector_to_solutions(target_vec, spaces, sols_dest);

  delete [] target_vec;
}

void Lumped_Projection::project_lumped(Space * space,
                             Solution* sol_src, Solution* sol_dest,
                             MatrixSolverType matrix_solver, UMFPackMatrix*  mat )
{
  Hermes::vector<Space *> spaces;
  spaces.push_back(space);
  Hermes::vector<Solution *> sols_src;
  sols_src.push_back(sol_src);
  Hermes::vector<Solution *> sols_dest;
  sols_dest.push_back(sol_dest);  
  project_lumped(spaces, sols_src, sols_dest, matrix_solver,  mat);
}

void Lumped_Projection::project_lumped(Hermes::vector<Space *> spaces,
                             Hermes::vector<WeakForm::MatrixFormVol *> mfvol,
                             Hermes::vector<WeakForm::VectorFormVol *> vfvol,
                             Hermes::vector<MeshFunction*> source_meshfns,
                             scalar* target_vec, MatrixSolverType matrix_solver, UMFPackMatrix*  mat )
{
  _F_
  unsigned int n = spaces.size();
  unsigned int n_biforms = mfvol.size();
  if (n_biforms == 0)
    error("Please use the simpler version of project_global with the argument Hermes::vector<ProjNormType> proj_norms if you do not provide your own projection norm.");
  if (n_biforms != vfvol.size())
    error("Mismatched numbers of projection forms in project_global().");
  if (n != n_biforms)
    error("Mismatched numbers of projected functions and projection forms in project_global().");

  // This is needed since spaces may have their DOFs enumerated only locally
  // when they come here.
  int ndof = Space::assign_dofs(spaces);

  // Define projection weak form.
  WeakForm* proj_wf = new WeakForm(n);
  for (unsigned int i = 0; i < n; i++) {
    proj_wf->add_matrix_form(mfvol[i]);
    // FIXME
    // proj_wf->add_vector_form(i, proj_liforms[i].first, proj_liforms[i].second, HERMES_ANY, source_meshfns[i]);
  }

  project_internal(spaces, proj_wf, target_vec, matrix_solver,  mat);
	 delete proj_wf;
}



void Lumped_Projection::project_lumped_rhs(Hermes::vector<Space *> spaces, Hermes::vector<Solution *> source_sols,
                             scalar* target_vec, MatrixSolverType matrix_solver, UMFPackMatrix*  mat ){
  Hermes::vector<MeshFunction *> mesh_fns;
  for(unsigned int i = 0; i < source_sols.size(); i++)
    mesh_fns.push_back(source_sols[i]);

	   int n = spaces.size();
  // define temporary projection weak form
  WeakForm* wf = new WeakForm(n);
 	int i = 0;
  	ProjectionLumpedMatrixFormVol* form = new ProjectionLumpedMatrixFormVol(i, i);
  	ProjectionLumpedVectorFormVol* vector = new ProjectionLumpedVectorFormVol(i, mesh_fns[i]);
	  wf->add_matrix_form(form);    
    	wf->add_vector_form(vector);

  // this is needed since spaces may have their DOFs enumerated only locally.
  int ndof = Space::assign_dofs(spaces);

  // Initialize DiscreteProblem.
  DiscreteProblem* dp = new DiscreteProblem(wf, spaces);	
  	Vector* rhs = create_vector(matrix_solver);
	 dp->assemble(NULL, rhs);

    for (int i=0; i < ndof; i++){
		 if(mat->get(i,i)!= 0.0) target_vec[i] = rhs->get(i)/mat->get(i,i);
		else target_vec[i] = 0.0;
	}
  

  delete rhs;
  delete dp;	
	delete form;
 delete vector; 
  delete wf;



}

