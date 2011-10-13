#define HERMES_REPORT_ALL
#include "definitions.h"
#include "lumped_projection.h"
#include "p_only_adapt.h"


// 1. Step: (M_L/tau -theta(K+D)) u^L =   (M_L/tau + (1-theta)(K+D)) u^n
// 2. Step : f_ij = (M_c)_ij (dt_u_L(i)- dt_u_L(j)) + D_ij (u_L(i)- u_L(j)); f_i = sum_(j!=i) alpha_ij f_ij
// 3. Step:  M_L u^(n+1) = M_L u^L + tau * f 

const int INIT_REF_NUM =5;                   // Number of initial refinements.
const int P_INIT = 1;       
const int P_MAX = 4;                        // Initial polynomial degree.
const double time_step = 1e-3;                           // Time step.

const double T_FINAL = 2*PI;                       // Time interval length.
//const double T_FINAL = PI*0.5; 

const double NEWTON_TOL = 1e-5;                   // Stopping criterion for the Newton's method.
const int NEWTON_MAX_ITER = 20;                  // Maximum allowed number of Newton iterations.

const double P_ADAP_TOL_LS = 0.9;
const double P_ADAP_TOL_EX = 0.8;     //*err_max
const double P_ADAP_TOL_EX_2 = 0.5;   //*Anzahl Elementen (moeglichkeit zum verfeinern)
const int P_ADAP_MAX_ITER = 5;

const double EPS = 1e-10;

const int UNREF_FREQ = 5000;                         // Every UNREF_FREQth time step the mesh is derefined.


const double theta = 0.5;    // theta-Schema fuer Zeitdiskretisierung (theta =0 -> explizit, theta=1 -> implizit)

MatrixSolverType matrix_solver = SOLVER_UMFPACK; 

// Boundary markers.
const std::string BDY_IN = "inlet";
const std::string BDY_OUT = "outlet";

// Weak forms & Projection with masslumping
#include "extrapolation.cpp"
#include "patch_solution.cpp"
#include "least_square.cpp"


//FCT & p-Adaptivity
#include "fct.cpp"
#include "p_adapt.cpp"




int main(int argc, char* argv[])
{
  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("domain.mesh", &mesh);

  // Perform initial mesh refinements (optional).
  for (int i=0; i < INIT_REF_NUM; i++) mesh.refine_all_elements();

  // Initialize boundary conditions.
  DefaultEssentialBCConst bc_essential(BDY_IN, 0.0);
  EssentialBCs bcs(&bc_essential);

  
  // Create an H1 space with default shapeset.
  H1Space space(&mesh, &bcs, P_INIT);

  int ndof = space.get_num_dofs();

  info("ndof = %d", ndof);

 // BaseView mview("Hello world!", new WinGeom(700, 700, 500, 500));
  //mview.show(&space);

 // Initialize solution of lower & higher order
  Solution low_sln, ref_sln;

  // Previous time level solution (initialized by the initial condition).
  CustomInitialCondition u_prev_time(&mesh);  

  // Initialize the weak formulation.
CustomWeakFormMassmatrix massmatrix(time_step, &u_prev_time);
CustomWeakFormConvection convection(&u_prev_time);

  // Output solution in VTK format.
Linearizer lin;
bool mode_3D = true;
lin.save_solution_vtk(&u_prev_time, "/space/melli/pics_hermes/pics_padapt/init_padap_neu.vtk", "u", mode_3D);

  // Initialize views.
	ScalarView Lowview("niedriger Ordnung", new WinGeom(500, 500, 500, 400));
	//Lowview.show(&u_prev_time, HERMES_EPS_HIGH);
	ScalarView sview("Solution", new WinGeom(0, 500, 500, 400));
	//sview.show(&u_prev_time, HERMES_EPS_HIGH); 
	ScalarView pview("projezierter Anfangswert", new WinGeom(500, 0, 500, 400));
	OrderView mview("mesh", new WinGeom(0, 0, 500, 400));
	mview.show(&space);

	




		  // Initialize the FE problem.

	UMFPackMatrix* mass_matrix = new UMFPackMatrix;   //M_c/tau
	UMFPackMatrix* conv_matrix = new UMFPackMatrix;   //K
	UMFPackMatrix* low_matrix = new UMFPackMatrix;  
	UMFPackMatrix* lowmat_rhs = new UMFPackMatrix; 
	scalar* u_L = NULL; 

	std::list<int> p1_elements;
	int* elements_to_refine = new int[space.get_mesh()->get_max_element_id()+1];   // 2 = refine, 1 = nothing, 0 = coarse
	double* elements_error_ls_p2 = new double[space.get_mesh()->get_max_element_id()+1];
	double* elements_error_ex = new double[space.get_mesh()->get_max_element_id()+1];	
	double* elements_error_ls = new double[space.get_mesh()->get_max_element_id()+1];

	int* elements_to_refine_old = new int[space.get_mesh()->get_max_element_id()+1];







// Time stepping loop:
	double current_time = 0.0; 
	int ts = 1;
	char title[100];
	bool changed = true;
	int ps = 1;  //p-adapt-schritte
	int max_p = P_MAX;

	


		//Durchlaufen und p bestimmen: in dof_list dof von vertex-basisfunction speichern 
	AsmList al;	
	AsmList dof_list;
	int ref_ndof;


//Timestep loop
do
{	 info(" Time step %d, time %3.5f", ts, current_time); 

	H1Space* ref_space = new H1Space(&mesh, &bcs, P_INIT);
	PonlyAdapt* adapting = new PonlyAdapt(ref_space, HERMES_L2_NORM);
	DiscreteProblem* dp_mass = new DiscreteProblem(&massmatrix, ref_space);
	DiscreteProblem* dp_convection = new DiscreteProblem(&convection, ref_space);

	ps=1; 
	changed = true;
	for(int i = 0; i <= space.get_mesh()->get_max_element_id(); i++) elements_to_refine[i]= 2;
		while((changed ==true)&&(ps<=P_ADAP_MAX_ITER)){
			info("init- adap- step %d, timestep %d", ps, ts);
			ref_ndof = ref_space->get_num_dofs();
			scalar* coeff_vec = new scalar[ref_ndof];
			for(int i=0; i<ref_ndof;i++) coeff_vec[i]=0.0;		
				OGProjection::project_global(ref_space,&u_prev_time, coeff_vec, matrix_solver, HERMES_L2_NORM); 

			 changed = p_adap(ref_space, &u_prev_time, coeff_vec, P_ADAP_TOL_LS,P_ADAP_TOL_EX, P_ADAP_TOL_EX_2,adapting,elements_to_refine, NULL,
										elements_error_ex,elements_error_ls, elements_error_ls_p2, &p1_elements, max_p);
			mview.show(ref_space);
			delete [] coeff_vec;
			ps++;
		}
	changed =true;
	ref_ndof = ref_space->get_num_dofs();

	
		for(int i = 0; i <= space.get_mesh()->get_max_element_id(); i++){
						elements_to_refine_old[i]= elements_to_refine[i];
						 elements_to_refine[i]= 2;
		}

		scalar* coeff_vec = new scalar[ref_ndof];
		for(int i=0; i<ref_ndof;i++) coeff_vec[i]=0.0;	
		p1_list(ref_space, &dof_list, &al,&p1_elements );
		ps=1; 

//Adaptivity loop
	do
	{		info(" adap- step %d, timestep %d", ps, ts); 
			scalar* lumped_scalar = new scalar[ref_ndof];
			UMFPackVector* vec_rhs = new UMFPackVector(ref_ndof);
			scalar* coeff_vec_2 = new scalar[ref_ndof];
			for(int i=0; i<ref_ndof;i++) coeff_vec_2[i]=0.0;
			scalar* flux_scalar = new scalar[ref_ndof]; 
		scalar* P_plus = new scalar[ref_ndof]; scalar* P_minus = new scalar[ref_ndof];
		scalar* Q_plus = new scalar[ref_ndof]; scalar* Q_minus = new scalar[ref_ndof];	
		scalar* R_plus = new scalar[ref_ndof]; scalar* R_minus = new scalar[ref_ndof];	

				//----------------------MassLumping M_L/tau--------------------------------------------------------------------
			  // Set up the matrix, and rhs according to the solver selection.=>For Masslumping
		//	info("Masslumping");

			dp_mass->assemble(mass_matrix, NULL); 	
			UMFPackMatrix* lumped_matrix = massLumping(&dof_list,mass_matrix);

			//------------------------artificial DIFFUSION D---------------------------------------
			  // Set up the solver, matrix, and rhs according to the solver selection.=>artificial Diffusion
			//info("artificialDiffusion");

			dp_convection->assemble(conv_matrix, NULL,true);
			UMFPackMatrix* diffusion = artificialDiffusion(&dof_list,conv_matrix);

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
			low_matrix->add_matrix(lumped_matrix);  //kann nun fuer alle Zeitschritte verwendet werden (ausser bei Adaptivitaet)
			//M_L/tau+(1-theta)(K+D)
			lowmat_rhs->add_matrix(lumped_matrix);	

			lumped_matrix->multiply_with_scalar(time_step);  // M_L
			mass_matrix->multiply_with_scalar(time_step);  // massmatrix = M_C
		
			// Project the initial condition on the FE space->coeff_vec	
		
				/*UMFPackMatrix* proj_matrix = pure_p1_massLumping(&dof_list,mass_matrix);
				Lumped_Projection::project_lumped_rhs(&space, &u_prev_time, coeff_vec, matrix_solver, proj_matrix);
				OGProjection::project_global(&space,&u_prev_time, coeff_vec_2, matrix_solver, HERMES_L2_NORM);		
				lumped_flux_limiter(&dof_list,mass_matrix, proj_matrix, coeff_vec, coeff_vec_2);
				*/		
	/*	if(ts==1){
			Lumped_Projection::project_lumped(ref_space, &u_prev_time, coeff_vec, matrix_solver, lumped_matrix);
					OGProjection::project_global(ref_space,&u_prev_time, coeff_vec_2, matrix_solver, HERMES_L2_NORM);
				lumped_flux_limiter(&dof_list,mass_matrix, lumped_matrix, coeff_vec, coeff_vec_2,
									P_plus, P_minus, Q_plus, Q_minus, R_plus, R_minus );
				Solution::vector_to_solution(coeff_vec, ref_space, &ref_sln);
			pview.show(&ref_sln);
			}else if(changed==true){*/
	if((changed==true)||(ts==1)) {	
		Lumped_Projection::project_lumped(ref_space, &u_prev_time, coeff_vec, matrix_solver, lumped_matrix);
				OGProjection::project_global(ref_space,&u_prev_time, coeff_vec_2, matrix_solver, HERMES_L2_NORM);
			lumped_flux_limiter(&dof_list,mass_matrix, lumped_matrix, coeff_vec, coeff_vec_2,
									P_plus, P_minus, Q_plus, Q_minus, R_plus, R_minus );
				//OGProjection::project_global(ref_space,&u_prev_time, coeff_vec, matrix_solver, HERMES_L2_NORM); 
				Solution::vector_to_solution(coeff_vec, ref_space, &ref_sln);
			pview.show(&ref_sln);
			}
				
	

	//-------------rhs lower Order M_L/tau+ (1-theta)(K+D) u^n------------		
			lowmat_rhs->multiply_with_vector(coeff_vec, lumped_scalar); 
			vec_rhs->zero(); vec_rhs->add_vector(lumped_scalar);
	//-------------------------solution of lower order------------	
				//info("solving low order solultion");
			  // Solve the linear system and if successful, obtain the solution. M_L/tau-theta(D+K) u^n+1=  M_L/tau+ (1-theta)(K+D) u^n
			//Solver* lowOrd = create_linear_solver(matrix_solver,low_matrix,vec_rhs);
			UMFPackLinearSolver* lowOrd = new UMFPackLinearSolver(low_matrix,vec_rhs);	
			if(lowOrd->solve()){ 
				u_L = lowOrd->get_solution();  
				Solution::vector_to_solution(u_L, ref_space, &low_sln);	
			  }else error ("Matrix solver failed.\n");
		//---------------------------------------antidiffusive fluxes-----------------------------------	
				//info("assemble fluxes");	
			 antidiffusiveFlux(&dof_list,mass_matrix,lumped_matrix,conv_matrix,diffusion,vec_rhs, u_L, flux_scalar, 
									P_plus, P_minus, Q_plus, Q_minus, R_plus, R_minus);		
	
			for(int i= 0; i<ref_ndof; i++) coeff_vec[i]= u_L[i]+ (flux_scalar[i]*time_step/lumped_matrix->get(i,i));
			 Solution::vector_to_solution(coeff_vec, ref_space, &ref_sln);

			ps++;

	//	if(ps<=P_ADAP_MAX_ITER)
			changed = p_adap(ref_space, &ref_sln, coeff_vec, P_ADAP_TOL_LS,P_ADAP_TOL_EX,P_ADAP_TOL_EX_2, adapting,elements_to_refine,elements_to_refine_old, elements_error_ex,elements_error_ls, elements_error_ls_p2, &p1_elements, max_p);
		//	else changed = false;

			if(changed==true){ 
					if(coeff_vec!=NULL){ delete [] coeff_vec; 	coeff_vec = NULL;}
					ref_ndof = ref_space->get_num_dofs();
					coeff_vec = new scalar[ref_ndof];
					for(int i=0; i<ref_ndof;i++) coeff_vec[i]=0.0;
					//Durchlaufen und p bestimmen: in dof_list dof von vertex-basisfunction speichern 
					p1_list(ref_space, &dof_list, &al,&p1_elements );
				}


			 // Visualize the solution.
	/* sprintf(title, "low_Ord Time %3.2f", current_time);
			  Lowview.set_title(title);
			 Lowview.show(&low_sln);	 
			  sprintf(title, "korrigierte Loesung: Time %3.2f", current_time);
			  sview.set_title(title);
			  sview.show(&ref_sln);*/
				mview.show(ref_space);


		//View::wait(HERMES_WAIT_KEYPRESS);

			  // Clean up.
				delete lowOrd; 
			delete lumped_matrix; 
			delete diffusion;
			delete[] lumped_scalar;  
			delete[] coeff_vec_2;
			delete[] flux_scalar; 
			delete vec_rhs;
		delete [] P_plus;
		delete [] P_minus;
		delete [] Q_plus;
		delete [] Q_minus;
		delete [] R_plus;
		delete [] R_minus;

		  low_matrix->free();
	 	lowmat_rhs->free();
	

	}
	while((changed ==true)&&(ps<=P_ADAP_MAX_ITER));



			 // Visualize the solution.
sprintf(title, "low_Ord Time %3.2f", current_time);
			  Lowview.set_title(title);
			 Lowview.show(&low_sln);	 
			  sprintf(title, "korrigierte Loesung: Time %3.2f", current_time);
			  sview.set_title(title);
			  sview.show(&ref_sln);
			

	  // Update global time.
  current_time += time_step;

  // Increase time step counter
  ts++;
    // Copy last reference solution into sln_prev_time.
    u_prev_time.copy(&ref_sln);
		if(coeff_vec!=NULL){ delete [] coeff_vec; 	coeff_vec = NULL;}
		delete ref_space;
		delete dp_convection;
		delete dp_mass; 
		delete adapting; 
}
while (current_time < T_FINAL);

lin.save_solution_vtk(&u_prev_time, "/space/melli/pics_hermes/pics_padapt/end_padap_neu.vtk", "solution", mode_3D);
/*sprintf(title, "low_Ord Time %3.2f", current_time);
			  Lowview.set_title(title);
			 Lowview.show(&low_sln);	 
			  sprintf(title, "korrigierte Loesung: Time %3.2f", current_time);
			  sview.set_title(title);
			  sview.show(&ref_sln);*/


	delete mass_matrix;  
	delete conv_matrix;
	delete low_matrix;
	delete lowmat_rhs;

	delete [] elements_to_refine;
	delete [] elements_error_ex;
	delete [] elements_error_ls;
	delete [] elements_error_ls_p2;



  // Wait for the view to be closed.
  View::wait();
  return 0;
}

