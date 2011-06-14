
#include "../hermes_common/common.h"
#include "function/solution.h"
#include "function/forms.h"
#include "hermes2d.h"

class Lumped_Projection
{
public:
  static void project_lumped(Hermes::vector<Space *> spaces, Hermes::vector<Solution*> source_sols,
                   scalar* target_vec, MatrixSolverType matrix_solver);


 
protected:
  static void project_internal(Hermes::vector<Space *> spaces, WeakForm *wf, scalar* target_vec,
                               MatrixSolverType matrix_solver);



  // Jacobian matrix (same as stiffness matrix since projections are linear).
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
                               MatrixSolverType matrix_solver = SOLVER_UMFPACK){

  unsigned int n = spaces.size();

  // sanity checks
  if (n <= 0 || n > 10) error("Wrong number of projected functions in project_internal().");
  for (unsigned int i = 0; i < n; i++) if(spaces[i] == NULL) error("this->spaces[%d] == NULL in project_internal().", i);
  if (spaces.size() != n) error("Number of spaces must match number of projected functions in project_internal().");

  // this is needed since spaces may have their DOFs enumerated only locally.
  int ndof = Space::assign_dofs(spaces);

  // Initialize DiscreteProblem.
  DiscreteProblem* dp = new DiscreteProblem(wf, spaces);

  //SparseMatrix* matrix = create_matrix(matrix_solver);
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
	}else error ("Matrix solver failed.\n");

  if (target_vec != NULL)
    for (int i=0; i < ndof; i++) target_vec[i] = coeff_vec[i];

  
  delete solver;
  delete matrix;
  delete rhs;
  delete dp;
  delete lumped_matrix;
  
}



void Lumped_Projection::project_lumped(Hermes::vector<Space *> spaces, Hermes::vector<Solution*> source_sols,
                   scalar* target_vec, MatrixSolverType matrix_solver){

  Hermes::vector<MeshFunction *> mesh_fns;
  for(unsigned int i = 0; i < source_sols.size(); i++)
    mesh_fns.push_back(source_sols[i]);

   int n = spaces.size();
  // define temporary projection weak form
  WeakForm* proj_wf = new WeakForm(n);

  for (int i = 0; i < n; i++) {  
    
    proj_wf->add_matrix_form(new ProjectionLumpedMatrixFormVol(i, i));    
    proj_wf->add_vector_form(new ProjectionLumpedVectorFormVol(i, mesh_fns[i]));
  }


  project_internal(spaces, proj_wf, target_vec, matrix_solver);
  delete proj_wf;
}

