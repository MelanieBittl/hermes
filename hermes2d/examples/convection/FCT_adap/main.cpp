#define HERMES_REPORT_ALL
#include "definitions.h"
#include "lumped_projection.h"

using namespace RefinementSelectors;

const int INIT_REF_NUM =3;                   // Number of initial refinements.
const int P_INIT = 1;                             // Initial polynomial degree.
double time_step = 1e-3;                           // Time step.
const double T_FINAL = 2*PI;                       // Time interval length.
//const double T_FINAL = 0.002;  

const double theta = 0.5;    // time-discretization (theta =0 -> explizit, theta=1 -> implizit)

MatrixSolverType matrix_solver = SOLVER_UMFPACK; 
MatrixSolverType matrix_solver_super = SOLVER_UMFPACK; // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

// Adaptivity
const int UNREF_FREQ = 2;                         // Every UNREF_FREQth time step the mesh is derefined.
const int UNREF_METHOD = 1;                       // 1... mesh reset to basemesh and poly degrees to P_INIT.   
                                                  // 2... one ref. layer shaved off, poly degrees reset to P_INIT.
                                                  // 3... one ref. layer shaved off, poly degrees decreased by one. 
const double THRESHOLD = 0.3;                      // This is a quantitative parameter of the adapt(...) function and
                                                  // it has different meanings for various adaptive strategies (see below).
const int STRATEGY = 1;                           // Adaptive strategy:
                                                  // STRATEGY = 0 ... refine elements until sqrt(THRESHOLD) times total
                                                  //   error is processed. If more elements have similar errors, refine
                                                  //   all to keep the mesh symmetric.
                                                  // STRATEGY = 1 ... refine all elements whose error is larger
                                                  //   than THRESHOLD times maximum element error.
                                                  // STRATEGY = 2 ... refine all elements whose error is larger
                                                  //   than THRESHOLD.
                                                  // More adaptive strategies can be created in adapt_ortho_h1.cpp.
const CandList CAND_LIST = H2D_H_ISO;          // Predefined list of element refinement candidates. Possible values are
                                                  // H2D_P_ISO, H2D_P_ANISO, H2D_H_ISO, H2D_H_ANISO, H2D_HP_ISO,
                                                  // H2D_HP_ANISO_H, H2D_HP_ANISO_P, H2D_HP_ANISO.
                                                  // See the User Documentation for details.
const int MESH_REGULARITY = 1;                   // Maximum allowed level of hanging nodes:
                                                  // MESH_REGULARITY = -1 ... arbitrary level hanging nodes (default),
                                                  // MESH_REGULARITY = 1 ... at most one-level hanging nodes,
                                                  // MESH_REGULARITY = 2 ... at most two-level hanging nodes, etc.
                                                  // Note that regular meshes are not supported, this is due to
                                                  // their notoriously bad performance.
const double CONV_EXP = 1.0;                      // Default value is 1.0. This parameter influences the selection of
                                                  // cancidates in hp-adaptivity. See get_optimal_refinement() for details.
const double ERR_STOP = 5.0;                      // Stopping criterion for adaptivity (rel. error tolerance between the
                                                  // fine mesh and coarse mesh solution in percent).
const int NDOF_STOP = 17000;                      // Adaptivity process stops when the number of degrees of freedom grows
                                                  // over this limit. This is to prevent h-adaptivity to go on forever.

const int ADAPSTEP_MAX = 100;			// maximale Anzahl an Adaptivitaetsschritten


// Boundary markers.
const std::string BDY_IN = "inlet";
const std::string BDY_OUT = "outlet";

#include "fct.cpp"


int main(int argc, char* argv[])
{
  // Load the mesh.
  Mesh mesh, basemesh;
  H2DReader mloader;
  mloader.load("domain.mesh", &basemesh);

  // Perform initial mesh refinements (optional).
  for (int i=0; i < INIT_REF_NUM; i++) basemesh.refine_all_elements();
	//basemesh.refine_towards_vertex(4, 5); //(int vertex_id, int depth)
  mesh.copy(&basemesh);


  // Initialize boundary conditions.
  DefaultEssentialBCConst bc_essential(BDY_IN, 0.0);
  EssentialBCs bcs(&bc_essential);

  // Create an H1 space with default shapeset.
  H1Space space(&mesh, &bcs, P_INIT);
  int ndof = space.get_num_dofs();
  info("ndof = %d", ndof);

 // Initialize coarse and reference mesh solution.
  Solution sln, ref_sln, low_sln;

  // Previous time level solution (initialized by the initial condition).
  CustomInitialCondition sln_prev_time(&mesh);
  ScalarView Aview("Initial condition", new WinGeom(0, 0, 500, 400));
  //Aview.show(&sln_prev_time);

  // Initialize the weak formulation.
CustomWeakFormMassmatrix massmatrix(time_step, &sln_prev_time);
CustomWeakFormConvection convection(&sln_prev_time);

  // Create a refinement selector.
  H1ProjBasedSelector selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);
	//selector.set_option(H2D_PREFER_SYMMETRIC_MESH, true);

  // Output solution in VTK format.
Linearizer lin;
bool mode_3D = true;
//lin.save_solution_vtk(&sln_prev_time, "/space/melli/pics_hermes/init_hadap.vtk", "u", mode_3D);

  // Visualize initial condition.
  char title[100];
  ScalarView view("Initial condition", new WinGeom(0, 500, 500, 400));
	//view.show(&sln_prev_time, HERMES_EPS_HIGH);
  OrderView ordview("Initial mesh", new WinGeom(550, 0, 500, 400)); 
 // ordview.show(&space);
  OrderView ord2view("ref mesh", new WinGeom(550, 500, 500, 400));


// Time stepping loop:
double current_time = 0.0; int ts = 1;
int ref_dof;
scalar* u_L = NULL;
  bool done ; int as;
    double err_est; double err_est_rel_total;

		 UMFPackMatrix* mass_matrix = new UMFPackMatrix;   //M_c/tau
		UMFPackMatrix* conv_matrix = new UMFPackMatrix;   //K
		UMFPackMatrix* low_matrix = new UMFPackMatrix;  
 		UMFPackMatrix* lowmat_rhs = new UMFPackMatrix; 
		 Adapt* adaptivity = new Adapt(&space, HERMES_L2_NORM);
int max_dof = ndof;
do
{
	  if (Space::get_num_dofs(&space) >= NDOF_STOP) mesh.copy(&basemesh);          

    // Periodic global derefinement. 
    if (ts > 1 && ts % UNREF_FREQ == 0) 
    {
      info("Global mesh derefinement.");
      switch (UNREF_METHOD) {
        case 1: mesh.copy(&basemesh);
                space.set_uniform_order(P_INIT);
                break;
        case 2: mesh.unrefine_all_elements();
                space.set_uniform_order(P_INIT);
                break;
        case 3: mesh.unrefine_all_elements();
                //space.adjust_element_order(-1, P_INIT);
                space.adjust_element_order(-1, -1, P_INIT, P_INIT);
                break;
        default: error("Wrong global derefinement method.");
      }

      ndof = Space::get_num_dofs(&space);
    }

   // Spatial adaptivity loop. Note: sln_prev_time must not be changed during spatial adaptivity. 
   done = false; as = 1;
    
 info("Current Time: %f, Time step %d, maxdof %d",current_time, ts, max_dof);
    do {    

	      // Construct globally refined reference mesh and setup reference space.
	      Space* ref_space = Space::construct_refined_space(&space,0);  // 0 erzwingt, dass kein Polynomgrad erhoeht wird
		ref_dof = Space::get_num_dofs(ref_space);
		if(ref_dof>max_dof) max_dof=ref_dof;
		  info("refdof = %d", ref_dof);
//ord2view.show(ref_space);
//info("adaptivity step %d: ndof: %d, ref_ndof: %d", as, Space::get_num_dofs(&space), Space::get_num_dofs(ref_space));
 //-----------------------------------------------------------------------------------		


	      // Initialize discrete problem on reference mesh.
		DiscreteProblem* dp_mass = new DiscreteProblem(&massmatrix, ref_space);
		DiscreteProblem* dp_convection = new DiscreteProblem(&convection, ref_space);  	
	//----------------------MassLumping M_L--------------------------------------------------------------------			  

		dp_mass->assemble(mass_matrix, NULL); //rechte Seite braucht nicht berechnet zu werden!
		UMFPackMatrix* lumped_matrix = massLumping(mass_matrix);
	//------------------------artificial DIFFUSION D---------------------------------------			

		dp_convection->assemble(conv_matrix, NULL,true);  //true => erzwingt Diagonleinträge ggf. = 0!		
		 UMFPackMatrix* diffusion = artificialDiffusion(conv_matrix);		
		//--------------------------------------------------------------------------------------------

		lowmat_rhs->create(conv_matrix->get_size(),conv_matrix->get_nnz(), conv_matrix->get_Ap(), conv_matrix->get_Ai(),conv_matrix->get_Ax());
		lowmat_rhs->add_matrix(diffusion); 	
		low_matrix->create(lowmat_rhs->get_size(),lowmat_rhs->get_nnz(), lowmat_rhs->get_Ap(), lowmat_rhs->get_Ai(),lowmat_rhs->get_Ax());
	//(-theta)(K+D)
		if(theta==0) low_matrix->zero();
		else	low_matrix->multiply_with_scalar(-theta);
	//(1-theta)(K+D)
		if(theta ==1) lowmat_rhs->zero();
		else lowmat_rhs->multiply_with_scalar((1.0-theta));		
	//M_L/tau - theta(D+K)
		low_matrix->add_matrix(lumped_matrix); 		
	//M_L/tau+(1-theta)(K+D)
		lowmat_rhs->add_matrix(lumped_matrix);
		
		mass_matrix->multiply_with_scalar(time_step);  // massmatrix = M_C
		scalar* coeff_vec = new scalar[ref_dof];
		memset(coeff_vec, 0, sizeof(scalar)*ref_dof);
		scalar* coeff_vec_2 = new scalar[ref_dof];				
		UMFPackVector* rhs = new UMFPackVector(ref_dof);
		
		scalar* P_plus = new scalar[ref_dof]; scalar* P_minus = new scalar[ref_dof];
		scalar* Q_plus = new scalar[ref_dof]; scalar* Q_minus = new scalar[ref_dof];	
		scalar* R_plus = new scalar[ref_dof]; scalar* R_minus = new scalar[ref_dof];	


//-------------rhs lower Order M_L/tau+ (1-theta)(K+D) u^n------------	 
		
		// Project the previous solution on the FE space		
		if(ts==1) Lumped_Projection::project_lumped(ref_space, &sln_prev_time, coeff_vec, matrix_solver);
		else{ 
			//OGProjection::project_global(ref_space,&sln_prev_time, coeff_vec, matrix_solver, HERMES_L2_NORM);
			Lumped_Projection::project_lumped(ref_space, &sln_prev_time, coeff_vec, matrix_solver);
			OGProjection::project_global(ref_space,&sln_prev_time, coeff_vec_2, matrix_solver_super, HERMES_L2_NORM);
			lumped_flux_limiter(mass_matrix, lumped_matrix, coeff_vec, coeff_vec_2,	P_plus, P_minus, Q_plus, Q_minus, R_plus, R_minus);
		}
		lowmat_rhs->multiply_with_vector(coeff_vec, coeff_vec_2);  
		rhs->zero(); rhs->add_vector(coeff_vec_2); 


//-------------------------Loesung niedriger Ordnung------------		
		// Solve the linear system and if successful, obtain the solution.
		//Solver* lowOrd = create_linear_solver(matrix_solver,low_matrix,rhs);	
		UMFPackLinearSolver* lowOrd = new UMFPackLinearSolver(low_matrix,rhs);
		if(lowOrd->solve()){ 
			u_L = lowOrd->get_solution(); 
		Solution::vector_to_solution(lowOrd->get_solution(), ref_space, &low_sln);	
		Solution::vector_to_solution(lowOrd->get_solution(), ref_space, &ref_sln);		
	  	}else error ("Matrix solver failed.\n");
//---------------------------------------antidiffusive fluxes-----------------------------------		
	antidiffusiveFlux(mass_matrix,lumped_matrix,conv_matrix,diffusion,rhs, u_L,coeff_vec,P_plus, P_minus, Q_plus, Q_minus, R_plus, R_minus); //Fluesse in coeff_vec	

//-----------------------korrigierte Loesung-----------------------------                   
	for(int i= 0; i<ref_dof; i++) coeff_vec_2[i]= u_L[i]+ (coeff_vec[i]/lumped_matrix->get(i,i));
	Solution::vector_to_solution(coeff_vec_2, ref_space, &ref_sln);
//------------------------------------------------------------	 
	      // Project the fine mesh solution onto the coarse mesh.
	      //info("Projecting fine mesh solution on coarse mesh for error estimation.");
	     OGProjection::project_global(&space, &ref_sln, &sln, matrix_solver_super, HERMES_L2_NORM); 

	      // Calculate element errors and total error estimate.
	      //info("Calculating error estimate.");
	    err_est_rel_total = adaptivity->calc_err_est(&sln, &ref_sln) * 100;

	      // Report results.
	      //info("ndof: %d, ref_ndof: %d, err_est_rel: %g%%", 
		   //Space::get_num_dofs(&space), Space::get_num_dofs(ref_space), err_est_rel_total);
		//info("err_est_rel: %g%%",err_est_rel_total);

	      // If err_est too large, adapt the mesh.
	      if (err_est_rel_total < ERR_STOP){
			//info("Adaptivity-Stop: err_est_rel_total");
		 	done = true;
			
	      }else if(as>=ADAPSTEP_MAX){
			//info("Adaptivity-Stop: ADAPSTEP_MAX");
		 	done = true;		
	     }else{ 	      
		//info("Adapting the coarse mesh.");
		done = adaptivity->adapt(&selector, THRESHOLD, STRATEGY, MESH_REGULARITY);		
		if (Space::get_num_dofs(&space) >= NDOF_STOP){ 
			//info("Adaptivity-Stop:NDOF_STOP,ndof: %d,refdof: %d ",Space::get_num_dofs(&space),Space::get_num_dofs(ref_space));
		  done = true;
		}else
		  // Increase the counter of performed adaptivity steps.
		  as++;
	      }
	
	   //view.show(&ref_sln);	
	

	 
	      // Clean up.
		//matrices			
		delete lumped_matrix; 		
		delete diffusion;
		//solver
		delete lowOrd;		
		//rest
      	delete ref_space;
      	delete dp_mass;
      	delete dp_convection;
		delete rhs;
		delete [] coeff_vec;
		delete [] coeff_vec_2;
		delete [] P_plus;
		delete [] P_minus;
		delete [] Q_plus;
		delete [] Q_minus;
		delete [] R_plus;
		delete [] R_minus;
		
      if(!done)
        delete ref_sln.get_mesh(); 

	  low_matrix->free();
	  lowmat_rhs->free();
    }
    while (done == false);


    // Copy last reference solution into sln_prev_time.
    sln_prev_time.copy(&ref_sln);
		
    // Visualize the solution and mesh.	      
		sprintf(title, "LowOrd-Solution, time %g", current_time);
	      Aview.set_title(title);
		 Aview.show_mesh(false);
	   Aview.show(&low_sln);		
	      sprintf(title, "Ref-Solution, time %g", current_time);
	      view.set_title(title);
		 view.show_mesh(false);	      
	      view.show(&ref_sln);		
	     sprintf(title, "coarse-mesh, time %g", current_time);
	     ordview.set_title(title);
	      ordview.show(&space);
  View::wait();
  // Update global time.
  current_time += time_step;

  // Increase time step counter
  ts++;
if (current_time < T_FINAL) delete ref_sln.get_mesh(); 
	
}
while (current_time < T_FINAL); 

// Visualize the solution and mesh.	      
	/*	sprintf(title, "LowOrd-Solution, time %g", current_time);
	      Aview.set_title(title);
		 Aview.show_mesh(false);
	   Aview.show(&low_sln);		
	      sprintf(title, "Ref-Solution, time %g", current_time);
	      view.set_title(title);
		 view.show_mesh(false);	      
	      view.show(&ref_sln);		
	     sprintf(title, "coarse-mesh, time %g", current_time);
	     ordview.set_title(title);
	      ordview.show(&space);*/
info("End-dofs, max: %d", max_dof);

lin.save_solution_vtk(&sln_prev_time, "/space/melli/pics_hermes/end_hadap.vtk", "u", mode_3D);
	    
		delete ref_sln.get_mesh(); 
		delete mass_matrix;
		delete conv_matrix;	
		delete low_matrix;
		delete lowmat_rhs;	
		delete adaptivity;

  // Wait for the view to be closed.
  View::wait();
  return 0;
}

