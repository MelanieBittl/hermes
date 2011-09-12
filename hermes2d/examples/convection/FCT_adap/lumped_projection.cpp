#include "../hermes_common/common.h"
#include "function/solution.h"
#include "function/forms.h"
#include "hermes2d.h"

class Lumped_Projection
{
public:
  static void project_lumped(Hermes::vector<Space *> spaces, Hermes::vector<MeshFunction *> source_meshfns,
                             scalar* target_vec, MatrixSolverType matrix_solver = SOLVER_UMFPACK);

  static void project_lumped(Hermes::vector<Space *> spaces, Hermes::vector<Solution *> source_sols,
                             scalar* target_vec, MatrixSolverType matrix_solver = SOLVER_UMFPACK);

  static void project_lumped(Space* space, MeshFunction* source_meshfn,
                             scalar* target_vec, MatrixSolverType matrix_solver = SOLVER_UMFPACK);

  static void project_lumped(Hermes::vector<Space *> spaces,
                             Hermes::vector<Solution*> sols_src, Hermes::vector<Solution*> sols_dest,
                             MatrixSolverType matrix_solver = SOLVER_UMFPACK);

  static void project_lumped(Space * space,
                             Solution* sol_src, Solution* sol_dest,
                             MatrixSolverType matrix_solver = SOLVER_UMFPACK);

  static void project_lumped(Hermes::vector<Space *> spaces,
                             Hermes::vector<WeakForm::MatrixFormVol *> mfvol,
                             Hermes::vector<WeakForm::VectorFormVol *> vfvol,
                             Hermes::vector<MeshFunction*> source_meshfns,
                             scalar* target_vec, MatrixSolverType matrix_solver = SOLVER_UMFPACK);

 
protected:
  static void project_internal(Hermes::vector<Space *> spaces, WeakForm *wf, scalar* target_vec,
                               MatrixSolverType matrix_solver);
  
  class ProjectionLumpedMatrixFormVol : public WeakForm::MatrixFormVol
  {
  public:
    ProjectionLumpedMatrixFormVol(int i, int j) : WeakForm::MatrixFormVol(i, j)
    {
      this->adapt_eval = false;      
    }

    scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v,
                 Geom<double> *e, ExtData<scalar> *ext) const
    {      
      
        return lumped_projection_biform<double, scalar>(n, wt, u_ext, u, v, e, ext);
     
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
            Geom<Ord> *e, ExtData<Ord> *ext) const
    {
          
        return lumped_projection_biform<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
     
    }

  private:

    template<typename Real, typename Scalar>
    static Scalar lumped_projection_biform(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u,
                                       Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
    {
      _F_
      Scalar result = 0;
      for (int i = 0; i < n; i++)
        result += wt[i] * (u->val[i] * v->val[i]);	
      return result;
    }

  };

  // Residual.
  class ProjectionLumpedVectorFormVol : public WeakForm::VectorFormVol
  {
  public:
    ProjectionLumpedVectorFormVol(int i, MeshFunction* ext) : WeakForm::VectorFormVol(i)
    {
      this->adapt_eval = false;     
      this->ext = Hermes::vector<MeshFunction *>();
      this->ext.push_back(ext);
    }

    scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v,
                 Geom<double> *e, ExtData<scalar> *ext) const
    {
     
      
        return lumped_projection_residual<double, scalar>(n, wt, u_ext, v, e, ext);
     
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
            Geom<Ord> *e, ExtData<Ord> *ext) const
    {
      
        return lumped_projection_residual<Ord, Ord>(n, wt, u_ext, v, e, ext);
    
    }

  private:

    template<typename Real, typename Scalar>
    Scalar lumped_projection_residual(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v,
                                         Geom<Real> *e, ExtData<Scalar> *ext) const
     {
       _F_
       Scalar result = 0;
       for (int i = 0; i < n; i++)
         //result += wt[i] * (u_ext[this->i]->val[i] - ext->fn[0]->val[i]) * v->val[i];
	result += wt[i] * (ext->fn[0]->val[i]) * v->val[i];
       return result;
    }

   
  };
};



void Lumped_Projection::project_internal(Hermes::vector<Space *> spaces, WeakForm *wf, scalar* target_vec,
                               MatrixSolverType matrix_solver ){

  unsigned int n = spaces.size();

  // sanity checks
  if (n <= 0 || n > 10) error("Wrong number of projected functions in project_internal().");
  for (unsigned int i = 0; i < n; i++) if(spaces[i] == NULL) error("this->spaces[%d] == NULL in project_internal().", i);
  if (spaces.size() != n) error("Number of spaces must match number of projected functions in project_internal().");

  // this is needed since spaces may have their DOFs enumerated only locally.
  int ndof = Space::assign_dofs(spaces);

  // Initialize DiscreteProblem.
  DiscreteProblem* dp = new DiscreteProblem(wf, spaces);

    	UMFPackMatrix* matrix = new UMFPackMatrix(); 
  	Vector* rhs = create_vector(matrix_solver);
	dp->assemble(matrix, rhs);  
//Masslumping
	 UMFPackMatrix* lumped_matrix = new UMFPackMatrix();   //M_L
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
	
	Solver* solver = create_linear_solver(matrix_solver, lumped_matrix, rhs);
  
  scalar* coeff_vec =NULL; 

	if(solver->solve()){ 
		coeff_vec = solver->get_solution();		
	}
	  else error ("Matrix solver failed.\n");

  if (target_vec != NULL)
    for (int i=0; i < ndof; i++) target_vec[i] = coeff_vec[i];

  
  delete solver;
  delete matrix;
  delete rhs;
  delete dp;
  delete lumped_matrix;
  
}




void Lumped_Projection::project_lumped(Hermes::vector<Space *> spaces, Hermes::vector<MeshFunction *> source_meshfns,
                             scalar* target_vec, MatrixSolverType matrix_solver)
{
  _F_
   int n = spaces.size();
  // define temporary projection weak form
  WeakForm* proj_wf = new WeakForm(n);
   if(n>1){ 
	for (int i = 0; i < n; i++) {    
    		proj_wf->add_matrix_form(new ProjectionLumpedMatrixFormVol(i, i));    
    		proj_wf->add_vector_form(new ProjectionLumpedVectorFormVol(i, source_meshfns[i]));
  	}
	 project_internal(spaces, proj_wf, target_vec, matrix_solver);
  }else{
 	int i = 0;
  	ProjectionLumpedMatrixFormVol* form = new ProjectionLumpedMatrixFormVol(i, i);
  	ProjectionLumpedVectorFormVol* vector = new ProjectionLumpedVectorFormVol(i, source_meshfns[i]);
	  proj_wf->add_matrix_form(form);    
    	proj_wf->add_vector_form(vector);
	 project_internal(spaces, proj_wf, target_vec, matrix_solver);
	  delete form;
         delete vector;
  }
 
 
  delete proj_wf;


}






void Lumped_Projection::project_lumped(Hermes::vector<Space *> spaces, Hermes::vector<Solution *> source_sols,
                             scalar* target_vec, MatrixSolverType matrix_solver )
{
  Hermes::vector<MeshFunction *> mesh_fns;
  for(unsigned int i = 0; i < source_sols.size(); i++)
    mesh_fns.push_back(source_sols[i]);
  project_lumped(spaces, mesh_fns, target_vec, matrix_solver);
}


void Lumped_Projection::project_lumped(Space* space, MeshFunction* source_meshfn,
                             scalar* target_vec, MatrixSolverType matrix_solver )
{
  Hermes::vector<Space *> spaces;
  spaces.push_back(space);
  Hermes::vector<MeshFunction *> source_meshfns;
  source_meshfns.push_back(source_meshfn);
  project_lumped(spaces, source_meshfns, target_vec, matrix_solver);
}


void Lumped_Projection::project_lumped(Hermes::vector<Space *> spaces,
                             Hermes::vector<Solution*> sols_src, Hermes::vector<Solution*> sols_dest,
                             MatrixSolverType matrix_solver )
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
                             MatrixSolverType matrix_solver )
{
  Hermes::vector<Space *> spaces;
  spaces.push_back(space);
  Hermes::vector<Solution *> sols_src;
  sols_src.push_back(sol_src);
  Hermes::vector<Solution *> sols_dest;
  sols_dest.push_back(sol_dest);  
  project_lumped(spaces, sols_src, sols_dest, matrix_solver);
}

void Lumped_Projection::project_lumped(Hermes::vector<Space *> spaces,
                             Hermes::vector<WeakForm::MatrixFormVol *> mfvol,
                             Hermes::vector<WeakForm::VectorFormVol *> vfvol,
                             Hermes::vector<MeshFunction*> source_meshfns,
                             scalar* target_vec, MatrixSolverType matrix_solver )
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

  project_internal(spaces, proj_wf, target_vec, matrix_solver);
	 delete proj_wf;
}

