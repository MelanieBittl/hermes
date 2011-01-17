#define HERMES_REPORT_INFO

// This file is part of Hermes2D.
//
// Hermes2D is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// Hermes2D is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Hermes2D.  If not, see <http://www.gnu.org/licenses/>.

#include "hermes2d.h"

// TODO: We do not take advantage of the fact that all blocks in the
// Jacobian matrix have the same structure. If the problem does not
// depend explicitly on time, then the blocks in the Jacobian matrix
// are the same up to a multiplicative constant.
void create_stage_wf(double current_time, double time_step, ButcherTable* bt,
                     DiscreteProblem* dp, WeakForm* stage_wf_right,
                     WeakForm* stage_wf_left)
{
  // Number of stages.
  int num_stages = bt->get_size();

  // Original weak formulation.
  WeakForm* wf = dp->get_weak_formulation();

  // Extract mesh from (the first space of) the discrete problem.
  Mesh* mesh = dp->get_space(0)->get_mesh();

  // Create a constant Solution to represent the stage time
  // stage_time = current_time + c_i*time_step.
  // (Temporary workaround. these should be passed as numbers.)
  Solution** stage_time_sol = new Solution*[num_stages];
  for (int i = 0; i < num_stages; i++) {
    stage_time_sol[i] = new Solution(mesh);
    stage_time_sol[i]->set_const(mesh, current_time + bt->get_C(i)*time_step);
  }

  // Extracting volume and surface matrix and vector forms from the
  // original weak formulation.
  if (wf->get_neq() != 1) error("wf->neq != 1 not implemented yet.");
  Hermes::vector<WeakForm::MatrixFormVol> mfvol_base = wf->get_mfvol();
  Hermes::vector<WeakForm::MatrixFormSurf> mfsurf_base = wf->get_mfsurf();
  Hermes::vector<WeakForm::VectorFormVol> vfvol_base = wf->get_vfvol();
  Hermes::vector<WeakForm::VectorFormSurf> vfsurf_base = wf->get_vfsurf();

  // Duplicate matrix volume forms, scale them according
  // to the Butcher's table, enhance them with additional
  // external solutions, and anter them as blocks to the
  // new stage Jacobian.
  for (unsigned int m = 0; m < mfvol_base.size(); m++) {
    WeakForm::MatrixFormVol mfv_base = mfvol_base[m];
    for (int i = 0; i < num_stages; i++) {
      for (int j = 0; j < num_stages; j++) {
        WeakForm::MatrixFormVol mfv_ij;
        mfv_ij.i = i;
        mfv_ij.j = j;
        mfv_ij.sym = mfv_base.sym;
        mfv_ij.area = mfv_base.area;
        mfv_ij.fn = mfv_base.fn;
        mfv_ij.ord = mfv_base.ord;
        std::copy(mfv_base.ext.begin(), mfv_base.ext.end(), mfv_ij.ext.begin());
        mfv_ij.scaling_factor = -time_step * bt->get_A(i, j);

        // Add stage_time_sol[i] as an external function to the form.
        mfv_ij.ext.push_back(stage_time_sol[i]);

        // Add the matrix form to the corresponding block of the
        // stage Jacobian matrix.
        stage_wf_right->add_matrix_form(&mfv_ij);
      }
    }
  }

  // Add mass volumetric form.
  WeakForm::MatrixFormVol mfv_00;
  mfv_00.i = 0;
  mfv_00.j = 0;
  mfv_00.sym = HERMES_SYM;
  mfv_00.area = HERMES_ANY;
  mfv_00.fn = l2_form<double, scalar>;
  mfv_00.ord = l2_form<Ord, Ord>;
  mfv_00.ext = Hermes::vector<MeshFunction*> ();
  mfv_00.scaling_factor = 1.0;
  stage_wf_left->add_matrix_form(&mfv_00);

  // Duplicate matrix surface forms, enhance them with
  // additional external solutions, and anter them as
  // blocks of the stage Jacobian.
  for (unsigned int m = 0; m < mfsurf_base.size(); m++) {
    WeakForm::MatrixFormSurf mfs_base = mfsurf_base[m];
    for (int i = 0; i < num_stages; i++) {
      for (int j = 0; j < num_stages; j++) {
        WeakForm::MatrixFormSurf mfs_ij;
        mfs_ij.i = i;
        mfs_ij.j = j;
        mfs_ij.area = mfs_base.area;
        mfs_ij.fn = mfs_base.fn;
        mfs_ij.ord = mfs_base.ord;
        std::copy(mfs_base.ext.begin(), mfs_base.ext.end(), mfs_ij.ext.begin());
        mfs_ij.scaling_factor = -time_step * bt->get_A(i, j);

        // Add stage_time_sol[i] as an external function to the form.
        mfs_ij.ext.push_back(stage_time_sol[i]);

        // Add the matrix form to the corresponding block of the
        // stage Jacobian matrix.
        stage_wf_right->add_matrix_form_surf(&mfs_ij);
      }
    }
  }

  // Duplicate vector volume forms, enhance them with
  // additional external solutions, and anter them as
  // blocks of the stage residual.
  for (unsigned int m = 0; m < vfvol_base.size(); m++) {
    WeakForm::VectorFormVol vfv_base = vfvol_base[m];
    for (int i = 0; i < num_stages; i++) {
      WeakForm::VectorFormVol vfv_i;
      vfv_i.i = i;
      vfv_i.area = vfv_base.area;
      vfv_i.fn = vfv_base.fn;
      vfv_i.ord = vfv_base.ord;
      std::copy(vfv_base.ext.begin(), vfv_base.ext.end(), vfv_i.ext.begin());
      vfv_i.scaling_factor = -1.0;

      // Add stage_time_sol[i] as an external function to the form.
      vfv_i.ext.push_back(stage_time_sol[i]);

      // Add the matrix form to the corresponding block of the
      // stage Jacobian matrix.
      stage_wf_right->add_vector_form(&vfv_i);
    }
  }

  // Duplicate vector surface forms, enhance them with
  // additional external solutions, and anter them as
  // blocks of the stage residual.
  for (unsigned int m = 0; m < vfsurf_base.size(); m++) {
    WeakForm::VectorFormSurf vfs_base = vfsurf_base[m];
    for (int i = 0; i < num_stages; i++) {
      WeakForm::VectorFormSurf vfs_i;
      vfs_i.i = i;
      vfs_i.area = vfs_base.area;
      vfs_i.fn = vfs_base.fn;
      vfs_i.ord = vfs_base.ord;
      std::copy(vfs_base.ext.begin(), vfs_base.ext.end(), vfs_i.ext.begin());
      vfs_i.scaling_factor = -1.0;

      // Add stage_time_sol[i] as an external function to the form.
      vfs_i.ext.push_back(stage_time_sol[i]);

      // Add the matrix form to the corresponding block of the
      // stage Jacobian matrix.
      stage_wf_right->add_vector_form_surf(&vfs_i);
    }
  }
}

// This takes a matrix, and uses it to formally construct a block-diagonal
// matrix. There are num_blocks blocks on the diagonal. The block diagonal
// matrix is then multiplied with the vector source_vec.
void multiply_as_diagonal_block_matrix(UMFPackMatrix* matrix, int num_blocks,
                                       scalar* source_vec, scalar* target_vec)
{
  int size = matrix->get_size();
  for (int i = 0; i < num_blocks; i++) {
    matrix->multiply(source_vec + i*size, target_vec + i*size);
  }
}

bool HERMES_RESIDUAL_AS_VECTOR_RK = true;
bool rk_time_step(double current_time, double time_step, ButcherTable* const bt,
                  scalar* coeff_vec, DiscreteProblem* dp, MatrixSolverType matrix_solver,
                  bool verbose, double newton_tol, int newton_max_iter,
                  double newton_damping_coeff, double newton_max_allowed_residual_norm)
{
  if (matrix_solver != SOLVER_UMFPACK)
    error("Sorry, rk_time_step() still only works with UMFpack.");
  if (dp->get_weak_formulation()->get_neq() > 1)
    error("Sorry, rk_time_step() does not work with systems yet.");

  // Matrix for the time derivative part of the equation (left-hand side).
  SparseMatrix* matrix_left = create_matrix(matrix_solver);
  //Vector* vector_left = create_vector(matrix_solver);

  // Matrix and vector for the rest (right-hand side).
  SparseMatrix* matrix_right = create_matrix(matrix_solver);
  Vector* vector_right = create_vector(matrix_solver);

  // Create matrix solver.
  Solver* solver = create_linear_solver(matrix_solver, matrix_right, vector_right);

  // Get number of stages from the Butcher's table.
  int num_stages = bt->get_size();

  // Get original space, mesh, and ndof.
  dp->get_space(0);
  Mesh* mesh = dp->get_space(0)->get_mesh();
  int ndof = dp->get_space(0)->get_num_dofs();

  // Create spaces for stage solutions. This is necessary
  // to define a num_stages x num_stages block weak formulation.
  Hermes::vector<Space*> stage_spaces;
  stage_spaces.push_back(dp->get_space(0));
  for (int i = 1; i < num_stages; i++)
    stage_spaces.push_back(dp->get_space(0)->dup(mesh));

  // Create a multistage weak formulation.
  WeakForm stage_wf_left;                   // For the matrix M (size ndof times ndof).
  WeakForm stage_wf_right(num_stages);      // For the rest of equation (written on the right),
                                            // size num_stages*ndof times num_stages*ndof.
  create_stage_wf(current_time, time_step, bt, dp, &stage_wf_right, &stage_wf_left);

  // Initialize discrete problems for the assembling of the
  // matrix M and the stage Jacobian matrix and residual.
  DiscreteProblem stage_dp_left(&stage_wf_left, dp->get_space(0));
  DiscreteProblem stage_dp_right(&stage_wf_right, stage_spaces);

  // Vector K_vector of length num_stages * ndof. will represent
  // the 'k_i' vectors in the usual R-K notation.
  scalar* K_vector = new scalar[num_stages*ndof];
  memset(K_vector, 0, num_stages * ndof * sizeof(scalar));

  // Vector u_prev_vec will represent y_n + h \sum_{j=1}^s a_{ij}k_i
  // in the usual R-K notation.
  scalar* u_prev_vec = new scalar[num_stages*ndof];

  // Vector for the left part of the residual.
  scalar* vector_left = new scalar[num_stages*ndof];

  // Prepare residuals of stage solutions.
  Hermes::vector<Solution*> residuals;
  Hermes::vector<bool> add_dir_lift;
  for (int i = 0; i < num_stages; i++) {
    residuals.push_back(new Solution(mesh));
    add_dir_lift.push_back(false);
  }

  // The Newton's loop.
  double residual_norm;
  int it = 1;
  while (true)
  {
    // Prepare vector Y_n + h\sum_{j=1}^s a_{ij} K_j.
    for (int idx = 0; idx < ndof; idx++) {
      scalar increment_i = 0;
      for (int i = 0; i < num_stages; i++) {
        for (int j = 0; j < num_stages; j++) {
          increment_i += bt->get_A(i, j) * K_vector[j*ndof + idx];
        }
        u_prev_vec[i*ndof + idx] = coeff_vec[idx] + time_step * increment_i;
      }
    }

    // Assemble the block-diagonal mass matrix M of size ndof times ndof.
    // The corresponding part of the global residual vector is obtained
    // just by multiplication.
    stage_dp_left.assemble(matrix_left);

    /*
      FILE* f = fopen("debug-left.txt", "w");
      matrix_left->dump(f, "tmp", DF_MATLAB_SPARSE);
      info("Matrix left dumped.");
      fclose(f);
      exit(0);
    */

    /*
    if (it == 2) {
      printf("K_vector before multiplication with M = ");
      //for (int i=0; i<ndof*num_stages; i++) printf("%g ", vector_left->get(i));
      for (int i=0; i<ndof*num_stages; i++) printf("%g ", K_vector[i]);
      printf("\n");
    }
    */

    multiply_as_diagonal_block_matrix((UMFPackMatrix*)matrix_left, num_stages,
                                      K_vector, vector_left);

    /*
    if (it == 2) {
      printf("vector_left (M times K_vector) = ");
      //for (int i=0; i<ndof*num_stages; i++) printf("%g ", vector_left->get(i));
      for (int i=0; i<ndof*num_stages; i++) printf("%g ", vector_left[i]);
      printf("\n");
    }
    */

    // Assemble the block Jacobian matrix of the stationary residual F
    // Diagonal blocks are created even if empty, so that matrix_left
    // can be added later.
    bool rhs_only = false;
    bool force_diagonal_blocks = true;
    stage_dp_right.assemble(u_prev_vec, matrix_right, vector_right,
                            rhs_only, force_diagonal_blocks);

    /*
      FILE* f = fopen("debug-right.txt", "w");
      matrix_right->dump(f, "tmp", DF_MATLAB_SPARSE);
      info("Matrix right dumped.");
      fclose(f);
    */

    /*
    // debug
    if (it == 2) {
      printf("vector_left = ");
      //for (int i=0; i<ndof*num_stages; i++) printf("%g ", vector_left->get(i));
      for (int i=0; i<ndof*num_stages; i++) printf("%g ", vector_left[i]);
      printf("\n");
      printf("vector_right = ");
      for (int i=0; i<ndof*num_stages; i++) printf("%g ", vector_right->get(i));
      printf("\n");
      //exit(0);

      // Debug.
      FILE* f = fopen("debug-left.txt", "w");
      matrix_left->dump(f, "tmp", DF_MATLAB_SPARSE);
      info("Matrix left dumped.");
      fclose(f);
      f = fopen("debug-right.txt", "w");
      matrix_right->dump(f, "tmp", DF_MATLAB_SPARSE);
      info("Matrix right dumped.");
      fclose(f);
    }
    */

    // Putting the two parts together into matrix_right and rhs_right.
    //((UMFPackMatrix*)matrix_right)->add_matrix((UMFPackMatrix*)matrix_left);
    ((UMFPackMatrix*)matrix_right)->add_to_diagonal_blocks(num_stages,
                                                           (UMFPackMatrix*)matrix_left);
    /*
      f = fopen("debug-composite.txt", "w");
      matrix_right->dump(f, "tmp", DF_MATLAB_SPARSE);
      info("Matrix composite dumped.");
      fclose(f);
      exit(0);
    */

    //((UMFPackMatrix*)matrix_right)->add_to_diagonal(1.0);
    vector_right->add_vector(vector_left);

    /*
    if (it == 2) {
      printf("vector_left = ");
      //for (int i=0; i<ndof*num_stages; i++) printf("%g ", vector_left->get(i));
      for (int i=0; i<ndof*num_stages; i++) printf("%g ", vector_left[i]);
      printf("\n");
      printf("vector_right = ");
      for (int i=0; i<ndof*num_stages; i++) printf("%g ", vector_right->get(i));
      printf("\n");
      exit(0);
    }
    */

    /*
    // Debug.
    if (it == -1) {
      FILE* f = fopen("debug-merged.txt", "w");
      matrix_right->dump(f, "tmp", DF_MATLAB_SPARSE);
      info("Merged matrix dumped.");
      fclose(f);
      exit(0);
    }
    */

    // Multiply the residual vector with -1 since the matrix
    // equation reads J(Y^n) \deltaY^{n+1} = -F(Y^n).
    vector_right->change_sign();

    // Measure the residual norm.
    if (HERMES_RESIDUAL_AS_VECTOR_RK) {
      // Calculate the l2-norm of residual vector.
      residual_norm = get_l2_norm(vector_right);
    }
    else {
      // Translate residual vector into residual functions.
      Solution::vector_to_solutions(vector_right, stage_dp_right.get_spaces(),
                                    residuals, add_dir_lift);
      residual_norm = calc_norms(residuals);
    }

    // Info for the user.
    if (verbose) info("---- Newton iter %d, ndof %d, residual norm %g",
                      it, ndof, residual_norm);

    // If maximum allowed residual norm is exceeded, fail.
    if (residual_norm > newton_max_allowed_residual_norm) {
      if (verbose) {
        info("Current residual norm: %g", residual_norm);
        info("Maximum allowed residual norm: %g", newton_max_allowed_residual_norm);
        info("Newton solve not successful, returning false.");
      }
      return false;
    }

    // If residual norm is within tolerance, or the maximum number
    // of iteration has been reached, then quit.
    if ((residual_norm < newton_tol || it > newton_max_iter) && it > 1) break;

    // Solve the linear system.
    if(!solver->solve()) error ("Matrix solver failed.\n");

    // Add \deltaY^{n+1} to Y^n.
    for (int i = 0; i < num_stages*ndof; i++) {
      K_vector[i] += newton_damping_coeff * solver->get_solution()[i];
      //printf("K_vector[%d] = %g\n", i, K_vector[i]);
    }
    //exit(0);

    // Increase iteration counter.
    it++;
  }

  // If max number of iterations was exceeded, fail.
  if (it >= newton_max_iter) {
    if (verbose) info("Maximum allowed number of Newton iterations exceeded, returning false.");
    return false;
  }

  /* THIS DID NOT WORK
  // Create a metrix solver for the equation M(Y_{n+1} - Y_n) = increment_vector.
  Vector* rk_increment_vector = create_vector(matrix_solver);
  rk_increment_vector->alloc(ndof);
  Solver* rk_increment_solver = create_linear_solver(matrix_solver,
                                matrix_left, rk_increment_vector);
  */

  // Calculate the vector \sum_{j=1}^s b_j k_j.
  Vector* rk_increment_vector = create_vector(matrix_solver);
  rk_increment_vector->alloc(ndof);
  for (int i = 0; i < ndof; i++) {
    rk_increment_vector->set(i, 0);
    for (int j = 0; j < num_stages; j++) {
      rk_increment_vector->add(i, bt->get_B(j) * K_vector[j*ndof + i]);
    }
  }

  /* THIS DID NOT WORK
  // Solve the linear system.
  if(!rk_increment_solver->solve()) error ("Matrix solver failed.\n");

  // Update coeff_vec to new time level.
  for (int i=0; i < ndof; i++) coeff_vec[i] += rk_increment_solver->get_solution()[i];
  */

  // Calculate Y^{n+1} = Y^n + h \sum_{j=1}^s b_j k_j.
  for (int i = 0; i < ndof; i++) coeff_vec[i] += time_step * rk_increment_vector->get(i);

  // Clean up.
  delete matrix_left;
  delete matrix_right;
  delete vector_right;
  delete solver;
  //delete rk_increment_solver;
  delete rk_increment_vector;

  // Delete stage spaces, but not the first (original) one.
  for (int i = 1; i < num_stages; i++) delete stage_spaces[i];

  // Delete all residuals.
  for (int i = 0; i < num_stages; i++) delete residuals[i];

  // TODO: Delete stage_wf, in particular its external solutions
  // stage_time_sol[i], i = 0, 1, ..., num_stages-1.

  // Delete stage_vec and u_prev_vec.
  delete [] K_vector;
  delete [] u_prev_vec;

  // debug
  delete [] vector_left;

  return true;
}