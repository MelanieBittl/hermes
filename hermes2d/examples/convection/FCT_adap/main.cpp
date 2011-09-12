#define HERMES_REPORT_ALL
#include "hermes2d.h"
 #define PI (3.141592653589793)     


using namespace RefinementSelectors;

const int INIT_REF_NUM =3;                   // Number of initial refinements.
const int P_INIT = 1;                             // Initial polynomial degree.
double time_step = 1e-3;                           // Time step.
const double T_FINAL = 2*PI;                       // Time interval length.
//const double T_FINAL = 0.1;  

const double theta = 0.5;    // time-discretization (theta =0 -> explizit, theta=1 -> implizit)

MatrixSolverType matrix_solver = SOLVER_UMFPACK; 

// Adaptivity
const int UNREF_FREQ = 100;                         // Every UNREF_FREQth time step the mesh is derefined.
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
const int NDOF_STOP = 6000;                      // Adaptivity process stops when the number of degrees of freedom grows
                                                  // over this limit. This is to prevent h-adaptivity to go on forever.

const int ADAPSTEP_MAX = 100;			// maximale Anzahl an Adaptivitaetsschritten


// Boundary markers.
const std::string BDY_IN = "inlet";
const std::string BDY_OUT = "outlet";

// Weak forms & Projection with masslumping
#include "definitions.cpp"
#include "lumped_projection.cpp"


//Mass lumping
UMFPackMatrix* massLumping(UMFPackMatrix* mass_matrix){
	 UMFPackMatrix* lumped_matrix = new UMFPackMatrix();   //M_L
	 int size = mass_matrix->get_size();
	 scalar diag[size];
	 int nnz = mass_matrix->get_nnz();
	 int row[size]; 
	int col[size+1];
	 for(int i = 0; i<size; i++){    
	    diag[i] = 0;
	    row[i]= i;
	    col[i]=i;
	 }
	col[size]=size;  //= Anzahl der Nichtnulleintraege
	
	 for(int i = 0; i<nnz; i++){    
	    diag[mass_matrix->get_Ai()[i]] += mass_matrix->get_Ax()[i]; 
	 }
	//for(int i = 0; i<size; i++) if(diag[i]<0.0)printf("kleiner 0");	
	 lumped_matrix->create(size, size, col, row, diag);  //lumped Matrix aufstellen
	return lumped_matrix;
}


//artificial Diffusion
UMFPackMatrix* artificialDiffusion(UMFPackMatrix* conv_matrix){
	 int size = conv_matrix->get_size();
	 int nnz = conv_matrix->get_nnz();
	scalar a,b;
	 UMFPackMatrix* diffusion = new UMFPackMatrix();  
	diffusion->create(size, nnz, conv_matrix->get_Ap(), conv_matrix->get_Ai(),conv_matrix->get_Ax());
	diffusion->zero();  //matrix = 0
	for(int i= 0; i<size; i++){
	  for(int j=(i+1);j<size;j++){
	       a= -conv_matrix->get(i,j);
	       b= -conv_matrix->get(j,i);     
	     if((a>=b)&&(a>0.0)){
		 diffusion->add(j,i,a);
		diffusion->add(i,j,a);	
		 diffusion->add(j,j,-a);
		diffusion->add(i,i,-a);	
	     }else if((b>a)&&(b>0.0)){
		 diffusion->add(i,j,b);
		diffusion->add(j,i,b);
	 	 diffusion->add(j,j,-b);
		diffusion->add(i,i,-b);
	     }
	  }
	}
	return diffusion;

}

void lumped_flux_limiter(UMFPackMatrix* mass_matrix,UMFPackMatrix* lumped_matrix, scalar* u_L, scalar* u_H){
	int ndof = mass_matrix->get_size();
	scalar P_plus[ndof], P_minus[ndof];	
	scalar Q_plus[ndof], Q_minus[ndof];	
	scalar R_plus[ndof], R_minus[ndof];	
	scalar alpha,f, plus, minus;
	//Berechnung von P&Q
	for(int i=0; i<ndof;i++){ P_plus[i]=0.0;P_minus[i]=0.0;Q_plus[i]=0.0;Q_minus[i]=0.0;}
	for(int i =0;i<ndof;i++){
		for(int j =(i+1); j<ndof;j++){		
			//f = flux[i][j];		
			f = mass_matrix->get(i,j)*(u_H[i]- u_H[j]);	
			if( (f*(u_H[j]- u_H[i])) > 0.0) f = 0.0; //prelimiting step	
			if(f>0.0)	{ 
				P_plus[i]+=f;
				P_minus[j]-=f;
			}else if (f<0.0){
			 	P_minus[i]+=f;
				P_plus[j]-=f;
			}
			f = lumped_matrix->get(i,i)*(u_L[j]-u_L[i])*time_step; //M_L/tau = lumped_matrix!
			if(f>Q_plus[i]) Q_plus[i] = f;				
			if(f<Q_minus[i]) Q_minus[i] = f;			
			f= lumped_matrix->get(j,j)*(u_L[i]-u_L[j])*time_step; //M_L/tau = lumped_matrix!
			if(f>Q_plus[j]) Q_plus[j] = f;	
			if(f<Q_minus[j]) Q_minus[j] = f;			
		}
	}
	//Berechnung von R
	for(int i=0; i<ndof;i++){
		plus = 1.0; minus = 1.0;		
		if(P_plus[i]!=0.0)  plus = Q_plus[i]/P_plus[i];		
		if(P_minus[i]!=0.0) minus = Q_minus[i]/P_minus[i];			
		if(plus>=1.0) R_plus[i]= 1.0;
		else 	     R_plus[i]= plus;
		if(minus>=1.0) R_minus[i]= 1.0;
		else 	     R_minus[i]= minus;		
	}	
	//Berechnung von alpha & f_i
	
	for(int i =0;i<ndof;i++){
		for(int j =(i+1); j<ndof;j++){				
			f = mass_matrix->get(i,j)*(u_H[i]- u_H[j]);			
			alpha = 1.0;
			if( (f*(u_H[j]- u_H[i])) > 0.0) f = 0.0; //prelimiting step
			if(f>=0.0){					
				if(R_plus[i]>R_minus[j]) alpha = R_minus[j];
				else 	alpha = R_plus[i];
			}else {
				if(R_minus[i]>R_plus[j]) alpha = R_plus[j];
				else 	alpha = R_minus[i]; 
			}
			
			u_L[i] += alpha*f/(lumped_matrix->get(i,i)*time_step);
			u_L[j] -= alpha*f/(lumped_matrix->get(j,j)*time_step);				
		}
	}
}



//Berechnung der antiduffisiven Fluesse
//f_ij und alpha_ij werden nicht explizit berechnet!! da scalar** flux = new_matrix<scalar>(ndof,ndof); zuviel Speicher braucht
void antidiffusiveFlux(UMFPackMatrix* mass_matrix,UMFPackMatrix* lumped_matrix,UMFPackMatrix* conv_matrix,UMFPackMatrix* diffusion,UMFPackVector* flux_dt_rhs, scalar* u_L, scalar* flux_scalar ){
	int ndof = conv_matrix->get_size();		
	Solver* flux_dt;
	scalar* dt_u_L = NULL;
	scalar P_plus[ndof], P_minus[ndof];	
	scalar Q_plus[ndof], Q_minus[ndof];	
	scalar R_plus[ndof], R_minus[ndof];	
	scalar alpha,f, plus, minus;

	conv_matrix->multiply_with_vector(u_L, flux_scalar);
	flux_dt_rhs->zero(); flux_dt_rhs->add_vector(flux_scalar);  //K u^L	
	flux_dt = create_linear_solver(matrix_solver,mass_matrix,flux_dt_rhs); //M_c u_t = K u^L
	if(flux_dt->solve())	dt_u_L = flux_dt->get_solution();	
	  else error ("Matrix solver failed.\n");

	//Berechnung von P&Q
	for(int i=0; i<ndof;i++){ P_plus[i]=0.0;P_minus[i]=0.0;Q_plus[i]=0.0;Q_minus[i]=0.0;}
	for(int i =0;i<ndof;i++){
		for(int j =(i+1); j<ndof;j++){		
			//f = flux[i][j];		
			f = mass_matrix->get(i,j)*(dt_u_L[i]- dt_u_L[j]) + diffusion->get(i,j)*(u_L[i]- u_L[j]);
			if( (f*(u_L[j]- u_L[i])) > 0.0) f = 0.0;  //prelimiting step
			if(f>0.0)	{ 
				P_plus[i]+=f;
				P_minus[j]-=f;
			}else if (f<0.0){
			 	P_minus[i]+=f;
				P_plus[j]-=f;
			}
			f = lumped_matrix->get(i,i)*(u_L[j]-u_L[i]); //M_L/tau = lumped_matrix!
			if(f>Q_plus[i]) Q_plus[i] = f;				
			if(f<Q_minus[i]) Q_minus[i] = f;			
			f= lumped_matrix->get(j,j)*(u_L[i]-u_L[j]); //M_L/tau = lumped_matrix!
			if(f>Q_plus[j]) Q_plus[j] = f;	
			if(f<Q_minus[j]) Q_minus[j] = f;			
		}
	}
	//Berechnung von R
	for(int i=0; i<ndof;i++){
		plus = 1.0; minus = 1.0;		
		if(P_plus[i]!=0.0)  plus = Q_plus[i]/P_plus[i];		
		if(P_minus[i]!=0.0) minus = Q_minus[i]/P_minus[i];			
		if(plus>=1.0) R_plus[i]= 1.0;
		else 	     R_plus[i]= plus;
		if(minus>=1.0) R_minus[i]= 1.0;
		else 	     R_minus[i]= minus;		
	}	
	//Berechnung von alpha & f_i
	for(int i=0; i<ndof;i++) flux_scalar[i]=0.0;
	for(int i =0;i<ndof;i++){
		for(int j =(i+1); j<ndof;j++){				
			f= mass_matrix->get(i,j)*(dt_u_L[i]- dt_u_L[j]) + diffusion->get(i,j)*(u_L[i]- u_L[j]);	
			if( (f*(u_L[j]- u_L[i])) > 0.0) f = 0.0;  //prelimiting step
			alpha = 1.0;
			if(f>0.0){					
				if(R_plus[i]>R_minus[j]) alpha = R_minus[j];
				else 	alpha = R_plus[i];
			}else if (f<0.0){
				if(R_minus[i]>R_plus[j]) alpha = R_plus[j];
				else 	alpha = R_minus[i]; 
			}
			
			flux_scalar[i] += alpha*f;
			flux_scalar[j] -= alpha*f;				
		}
	}

	//Cleanup	
	delete flux_dt;
}



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


  // Visualize initial condition.
  char title[100];
  ScalarView view("Initial condition", new WinGeom(0, 500, 500, 400));
  OrderView ordview("Initial mesh", new WinGeom(550, 0, 500, 400)); 
  //ordview.show(&space);
  OrderView ord2view("ref mesh", new WinGeom(550, 500, 500, 400));
 


// Time stepping loop:
double current_time = 0.0; int ts = 1;
int ref_dof;
scalar* u_L = NULL;
  bool done ; int as;
    double err_est; double err_est_rel_total;

		  UMFPackMatrix* mass_matrix = new UMFPackMatrix();   //M_c/tau
		UMFPackMatrix* conv_matrix = new UMFPackMatrix();   //K
		 UMFPackMatrix* low_matrix = new UMFPackMatrix();  
 		UMFPackMatrix* lowmat_rhs = new UMFPackMatrix(); 
		 Adapt* adaptivity = new Adapt(&space, HERMES_L2_NORM);

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
    
 info("Current Time: %f, Time step %d",current_time, ts);
    do {    

	      // Construct globally refined reference mesh and setup reference space.
	      Space* ref_space = Space::construct_refined_space(&space,0);  // 0 erzwingt, dass kein Polynomgrad erhoeht wird
		ref_dof = Space::get_num_dofs(ref_space);
//info("adaptivity step %d: ndof: %d, ref_ndof: %d", as, Space::get_num_dofs(&space), Space::get_num_dofs(ref_space));
 
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
		scalar* coeff_vec_2 = new scalar[ref_dof];				
		UMFPackVector* rhs = new UMFPackVector(ref_dof);

//-------------rhs lower Order M_L/tau+ (1-theta)(K+D) u^n------------	 
		
		// Project the previous solution on the FE space
		for(int i=0; i<ref_dof;i++) coeff_vec[i]=0.0; 
		if(ts==1) Lumped_Projection::project_lumped(ref_space, &sln_prev_time, coeff_vec, matrix_solver);
		else{ 
			//OGProjection::project_global(ref_space,&sln_prev_time, coeff_vec, matrix_solver, HERMES_L2_NORM);
			Lumped_Projection::project_lumped(ref_space, &sln_prev_time, coeff_vec, matrix_solver);
			 OGProjection::project_global(ref_space,&sln_prev_time, coeff_vec_2, matrix_solver, HERMES_L2_NORM);
			lumped_flux_limiter(mass_matrix, lumped_matrix, coeff_vec, coeff_vec_2);
		}
		lowmat_rhs->multiply_with_vector(coeff_vec, coeff_vec_2);  
		rhs->zero(); rhs->add_vector(coeff_vec_2); 
//-----------------------------------------------------------------------------------		
   // Now we can deallocate the previous fine mesh.
      		if(as > 1) delete ref_sln.get_mesh();  

//-------------------------Loesung niedriger Ordnung------------		
		// Solve the linear system and if successful, obtain the solution.
		Solver* lowOrd = create_linear_solver(matrix_solver,low_matrix,rhs);	
		if(lowOrd->solve()){ 
			u_L = lowOrd->get_solution(); 
		Solution::vector_to_solution(lowOrd->get_solution(), ref_space, &low_sln);	
		//Solution::vector_to_solution(lowOrd->get_solution(), ref_space, &ref_sln);		
	  	}else error ("Matrix solver failed.\n");
//---------------------------------------antidiffusive fluxes-----------------------------------		
	antidiffusiveFlux(mass_matrix,lumped_matrix,conv_matrix,diffusion,rhs, u_L,coeff_vec); //Fluesse in coeff_vec	

//-----------------------korrigierte Loesung-----------------------------                   
	for(int i= 0; i<ref_dof; i++) coeff_vec_2[i]= u_L[i]+ (coeff_vec[i]/lumped_matrix->get(i,i));
	Solution::vector_to_solution(coeff_vec_2, ref_space, &ref_sln);
//------------------------------------------------------------	 
	      // Project the fine mesh solution onto the coarse mesh.
	      //info("Projecting fine mesh solution on coarse mesh for error estimation.");
	     OGProjection::project_global(&space, &ref_sln, &sln, matrix_solver, HERMES_L2_NORM); 

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
		
	  
    }
    while (done == false);


    // Copy last reference solution into sln_prev_time.
    sln_prev_time.copy(&ref_sln);

   /*   // Visualize the solution and mesh.	      
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
	      ordview.show(&space);*/

  // Update global time.
  current_time += time_step;

  // Increase time step counter
  ts++;
}
while (current_time < T_FINAL); 

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

	    delete mass_matrix;	
				delete conv_matrix;	
		delete low_matrix;
		delete lowmat_rhs;
	     	delete adaptivity;

  

  // Wait for the view to be closed.
  View::wait();
  return 0;
}

