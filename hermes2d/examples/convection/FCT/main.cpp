#define HERMES_REPORT_ALL
#include "hermes2d.h"
 #define PI (3.141592653589793) 

// 1. Step: (M_L/tau -theta(K+D)) u^L =   (M_L/tau + (1-theta)(K+D)) u^n
// 2. Step : f_ij = (M_c)_ij (dt_u_L(i)- dt_u_L(j)) + D_ij (u_L(i)- u_L(j)); f_i = sum_(j!=i) alpha_ij f_ij
// 3. Step:  M_L u^(n+1) = M_L u^L + tau * f 




const int INIT_REF_NUM =6;                   // Number of initial refinements.
const int P_INIT = 1;                             // Initial polynomial degree.
const double time_step = 1e-3;                           // Time step.
const double T_FINAL = 2*PI;                       // Time interval length.


const double NEWTON_TOL = 1e-5;                   // Stopping criterion for the Newton's method.
const int NEWTON_MAX_ITER = 20;                  // Maximum allowed number of Newton iterations.


const double theta = 0.5;    // theta-Schema fuer Zeitdiskretisierung (theta =0 -> explizit, theta=1 -> implizit)

MatrixSolverType matrix_solver = SOLVER_UMFPACK; 

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
	col[size]=size;// = Anzahl der Nichtnulleintraege (bezieht sich auf theoretisch als naechstes kommenden Eintrag)
	
	 for(int i = 0; i<nnz; i++){    
	    diag[mass_matrix->get_Ai()[i]] += mass_matrix->get_Ax()[i]; 
	 }

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
//Assemble antidiffusive fluxes & Limiter
//f_ij und alpha_ij werden nicht explizit berechnet!! da scalar** flux = new_matrix<scalar>(ndof,ndof); zuviel Speicher braucht
void antidiffusiveFlux(UMFPackMatrix* mass_matrix,UMFPackMatrix* lumped_matrix,UMFPackMatrix* conv_matrix,UMFPackMatrix* diffusion,UMFPackVector* flux_dt_rhs, scalar* u_L, scalar* flux_scalar ){
	int ndof = conv_matrix->get_size();
	scalar flux_dt_scalar[ndof];	
	Solver* flux_dt;
	scalar* dt_u_L = NULL;
	scalar P_plus[ndof], P_minus[ndof];	
	scalar Q_plus[ndof], Q_minus[ndof];	
	scalar R_plus[ndof], R_minus[ndof];	
	scalar alpha,f, plus, minus;

	conv_matrix->multiply_with_vector(u_L, flux_dt_scalar);
	flux_dt_rhs->zero(); flux_dt_rhs->add_vector(flux_dt_scalar);  //K u^L	
	flux_dt = create_linear_solver(matrix_solver,mass_matrix,flux_dt_rhs); //M_c u_t = K u^L
	if(flux_dt->solve())	dt_u_L = flux_dt->get_solution();	
	  else error ("Matrix solver failed.\n");

	//Berechnung von P&Q
	for(int i=0; i<ndof;i++){ P_plus[i]=0.0;P_minus[i]=0.0;Q_plus[i]=0.0;Q_minus[i]=0.0;}
	for(int i =0;i<ndof;i++){
		for(int j =(i+1); j<ndof;j++){		
			//f = flux[i][j];		
			f = mass_matrix->get(i,j)*(dt_u_L[i]- dt_u_L[j]) + diffusion->get(i,j)*(u_L[i]- u_L[j]);
			if( (f*(u_L[j]- u_L[i])) > 0.0) f = 0.0; //prelimiting step
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
	alpha = 0.0;
	for(int i =0;i<ndof;i++){
		for(int j =(i+1); j<ndof;j++){				
			f= mass_matrix->get(i,j)*(dt_u_L[i]- dt_u_L[j]) + diffusion->get(i,j)*(u_L[i]- u_L[j]);	
			if( (f*(u_L[j]- u_L[i])) > 0.0) f = 0.0; //prelimiting step
			if(f>0){					
				if(R_plus[i]>R_minus[j]) alpha = R_minus[j];
				else 	alpha = R_plus[i];
			}else{
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


bool newton(scalar* coeff_vec, DiscreteProblem* dp, Solver* solver,  UMFPackMatrix* jac, UMFPackVector* rhs,
                            UMFPackVector* rhs_const,  double newton_tol, int newton_max_iter, 
                            bool verbose) 
{
//rhs_const = (M/tau + (1-theta)K)  *u^n  
//rhs = rhs_const + jac*Y^n
//jac  = (M/tau - theta*K) 
//Jac(Y^n) \deltaY^{n+1} = -F(Y^n) (=rhs).
double max_allowed_residual_norm = 1e6;
double damping_coeff = 1.0; 

Hermes2D hermes2d;

// Obtain the number of degrees of freedom.
    int ndof = dp->get_num_dofs();

scalar* rhs_var = new scalar[ndof];
for(int i=0; i<ndof;i++) rhs_var[i]=0.0;

  // The Newton's loop.
  double residual_norm;
  int it = 1;
  while (1)
  {	
	rhs->zero(); rhs->add_vector(rhs_const);
    	jac->multiply_with_vector(coeff_vec, rhs_var);
	for(int i = 0; i<ndof;i++) rhs_var[i] *= (-1.0);	
	rhs->add_vector(rhs_var); // rhs= -F(Y^n);
   
      // Calculate the l2-norm of residual vector
      //residual_norm = Hermes2D::get_l2_norm(rhs);
   residual_norm = hermes2d.get_l2_norm(rhs);

    // Info for the user.
    if (it == 1) {
      if (verbose) info("---- Newton initial residual norm: %g", residual_norm);
    }
    else if (verbose) info("---- Newton iter %d, residual norm: %g", it-1, residual_norm);

 // If maximum allowed residual norm is exceeded, fail.
    if (residual_norm > max_allowed_residual_norm) {
      if (verbose) {
        info("Current residual norm: %g", residual_norm);
        info("Maximum allowed residual norm: %g", max_allowed_residual_norm);
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
        for (int i = 0; i < ndof; i++) coeff_vec[i] += damping_coeff * solver->get_solution()[i];

    it++;
  }
  if (it >= newton_max_iter) {
    if (verbose) info("Maximum allowed number of Newton iterations exceeded, returning false.");
    return false;
  }
delete [] rhs_var;
  return true;
}


int main(int argc, char* argv[])
{
  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("domain.mesh", &mesh);

  // Perform initial mesh refinements (optional).
  for (int i=0; i < INIT_REF_NUM; i++) mesh.refine_all_elements();
//   MeshView mview("Hello world!", new WinGeom(0, 0, 350, 350));
 // mview.show(&mesh);

  // Initialize boundary conditions.
  DefaultEssentialBCConst bc_essential(BDY_IN, 0.0);
  EssentialBCs bcs(&bc_essential);

  // Create an H1 space with default shapeset.
  H1Space space(&mesh, &bcs, P_INIT);
  int ndof = space.get_num_dofs();
  info("ndof = %d", ndof);

 // Initialize solution of lower & higher order
  Solution low_sln;


  // Previous time level solution (initialized by the initial condition).
  CustomInitialCondition u_prev_time(&mesh);
  

  // Initialize the weak formulation.
CustomWeakFormMassmatrix massmatrix(time_step, &u_prev_time);
CustomWeakFormConvection convection(&u_prev_time);

  // Initialize the FE problem.
  DiscreteProblem dp_mass(&massmatrix, &space);
  DiscreteProblem dp_convection(&convection, &space);



//----------------------MassLumping M_L/tau--------------------------------------------------------------------

  // Set up the matrix, and rhs according to the solver selection.=>For Masslumping
  UMFPackMatrix* mass_matrix = new UMFPackMatrix();   //M_c/tau
  //Vector* mass_rhs = create_vector(matrix_solver);  
   //dp_mass.assemble(mass_matrix, mass_rhs); 
	dp_mass.assemble(mass_matrix, NULL); 	
UMFPackMatrix* lumped_matrix = massLumping(mass_matrix);



//------------------------artificial DIFFUSION D---------------------------------------

  // Set up the solver, matrix, and rhs according to the solver selection.=>artificial Diffusion
  UMFPackMatrix* conv_matrix = new UMFPackMatrix();   //K
  //Vector* conv_rhs = create_vector(matrix_solver);  
   //dp_convection.assemble(conv_matrix, conv_rhs,true);  //true => erzwingt Diagonleinträge ggf. = 0!
	dp_convection.assemble(conv_matrix, NULL,true);
 UMFPackMatrix* diffusion = artificialDiffusion(conv_matrix);

/*
//Test for Diffusion : dij=dji
for(int i=0; i<ndof;i++){
	for(int j=(i+1);j<ndof;j++){
		if(diffusion->get(i,j)!=diffusion->get(j,i) ){
			printf("(i,j)=(%i,%i)\n", i,j);
			printf("diff(i,j)= %f; diff(j,i)=%f", diffusion->get(i,j),diffusion->get(j,i) );
			
		}
	}
}
//Test for Diffusion : sum = 0 
scalar row_sum, col_sum;
for(int i=0; i<ndof;i++){
	row_sum = 0.0; col_sum = 0.0;
	for(int j=0;j<ndof;j++){
		row_sum += diffusion->get(i,j);
		col_sum += diffusion->get(j,i);	
	}
	if((row_sum < -1e-15)||(row_sum >1e-15)) printf("row= %i, sum =%f\n",i, row_sum);
	if((col_sum < -1e-15)||(col_sum >1e-15)) printf("col= %i, sum =%f\n",i, col_sum);	
}
*/
//--------------------------------------------------------------------------------------------
 UMFPackMatrix* low_matrix = new UMFPackMatrix();  
 UMFPackMatrix* K_D = new UMFPackMatrix();  
K_D->create(conv_matrix->get_size(),conv_matrix->get_nnz(), conv_matrix->get_Ap(), conv_matrix->get_Ai(),conv_matrix->get_Ax());
K_D->add_matrix(diffusion); 

/*  // Test fuer matrix_vector_multiplication
UMFPackMatrix* test_matrix = new UMFPackMatrix(); 
scalar* entry = new scalar[4];
int rows[4] = {0,1,0,1};
int cols[3] = {0,2,4};
for (int i=0; i<4; i++) entry[i] = (scalar) i;
test_matrix->create(2,4,cols,rows,entry);

// Print matrix 
for (int i=0;i<2;i++) {
  for (int j=0;j<2;j++)
    printf("%f ", test_matrix->get(i,j));
  printf("\n");
}
   
scalar* vector_in = new scalar[2];
for (int i=0;i<2;i++) vector_in[i] = 1.;

scalar* vector_out = new scalar[2];
test_matrix->multiply_with_vector(vector_in, vector_out);

for (int i=0;i<2;i++) printf("%f\n", vector_out[i]);
*/

/*
 

//Test for negative off-diagonal in K+D
for(int i=0; i<ndof;i++){
	for(int j=(i+1);j<ndof;j++){
		if(K_D->get(i,j)<0 ){
			printf("(i,j)=(%i,%i):::::", i,j);
			printf("K_D(i,j)= %f\n", K_D->get(i,j) );
			
		}
		if(K_D->get(j,i)<0 ){
			printf("(j,i)=(%i,%i):::::", j,i);
			printf("K_D(j,i)= %f\n", K_D->get(j,i) );
			
		}
	}
}*/

low_matrix->create(K_D->get_size(),K_D->get_nnz(), K_D->get_Ap(), K_D->get_Ai(),K_D->get_Ax());
//(-theta)(K+D)
if(theta==0) low_matrix->zero();
else	low_matrix->multiply_with_scalar(-theta);
//(1-theta)(K+D)
if(theta ==1) K_D->zero();
else K_D->multiply_with_scalar((1.0-theta));

/*
//Test fuer theta = 0.5
for(int i=0; i<ndof;i++){
	for(int j=0;j<ndof;j++){
		if(K_D->get(i,j)!=(-low_matrix->get(i,j)) ){
			printf("(i,j)=(%i,%i)\n", i,j);
			printf("KD= %f; low_matr=%f", K_D->get(i,j),low_matrix->get(i,j) );
			
		}
	}
}*/

//M_L/tau - theta(D+K)
low_matrix->add_matrix(lumped_matrix);  //kann nun fuer alle Zeitschritte verwendet werden (ausser bei Adaptivitaet)

//M_L/tau+(1-theta)(K+D)
K_D->add_matrix(lumped_matrix);


  // Initialize views.
   ScalarView Lowview("niedriger Ordnung", new WinGeom(500, 500, 500, 400));
  Lowview.show(&u_prev_time);
  ScalarView sview("Solution", new WinGeom(0, 500, 500, 400));
  sview.show(&u_prev_time); 

	
 	 // Previous time level solution (initialized by the initial condition).
 		 CustomInitialCondition u_prev_time_high(&mesh);
		ConvectionForm high_convection(time_step, &u_prev_time_high);
		  // Instantiate a class with global functions.
		  Hermes2D hermes2d;
		  // Project the initial condition on the FE space to obtain initial
		  // coefficient vector for the Newton's method.		
		  scalar* coeff_vec_newton = new scalar[ndof];
		  OGProjection::project_global(&space, &u_prev_time_high, coeff_vec_newton, matrix_solver);
		  // Initialize the FE problem.
		  DiscreteProblem dp(&high_convection, &space);
		  // Set up the solver, matrix, and rhs according to the solver selection.
		  SparseMatrix* matrix = create_matrix(matrix_solver);
		  Vector* rhs = create_vector(matrix_solver);
		  Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);
		   ScalarView Hermview("hohe Ordnung", new WinGeom(0, 0, 500, 400));
		  Hermview.show(&u_prev_time_high);
		  bool jacobian_changed = true;
		




// Time stepping loop:
double current_time = 0.0; 
int ts = 1;

scalar* lumped_scalar = new scalar[ndof];
UMFPackVector*  vec_rhs = new UMFPackVector(ndof);
scalar* coeff_vec = new scalar[ndof];
for(int i=0; i<ndof;i++) coeff_vec[i]=0.0;

scalar* u_L = NULL;
//scalar* fct_scalar = new scalar[ndof];
scalar* flux_scalar = new scalar[ndof];; 
Solver* lowOrd;
//Solver* fct;
char title[100];


UMFPackVector* newton_rhs = new UMFPackVector(ndof);
//newton_rhs->zero(); lumped_rhs->zero();
lowOrd = create_linear_solver(matrix_solver,low_matrix,vec_rhs);
//fct  = create_linear_solver(matrix_solver,lumped_matrix,vec_rhs); 
bool verbose = true;

// Project the initial condition on the FE space->coeff_vec		
	Lumped_Projection::project_lumped(&space, &u_prev_time, coeff_vec, matrix_solver);
	 mass_matrix->multiply_with_scalar(time_step);  // massmatrix = M_C

	
do
{
	  info(" Time step %d, time %3.5f", ts, current_time);  
//-----------------solution of higher order with Hermes  	
//info("------HERMES: ");
	// Perform Newton's iteration.		
	    if (!hermes2d.solve_newton(coeff_vec_newton, &dp, solver, matrix, rhs, jacobian_changed,
				       NEWTON_TOL, NEWTON_MAX_ITER, verbose)) error("Newton's iteration failed.");
	    // Update previous time level solution.
	    Solution::vector_to_solution(coeff_vec_newton, &space, &u_prev_time_high);
	

//-------------rhs lower Order M_L/tau+ (1-theta)(K+D) u^n------------		
	K_D->multiply_with_vector(coeff_vec, lumped_scalar); 
	vec_rhs->zero(); vec_rhs->add_vector(lumped_scalar);	

//-------------------------solution of lower order------------
	
	  // Solve the linear system and if successful, obtain the solution. M_L/tau-theta(D+K) u^n+1=  M_L/tau+ (1-theta)(K+D) u^n	
	if(lowOrd->solve()){ 
		u_L = lowOrd->get_solution();  
		if (ts > 1 && ts % 10 == 0) Solution::vector_to_solution(u_L, &space, &low_sln);	
	  }else error ("Matrix solver failed.\n");

	
  // if (!newton(u_L, &dp_mass, lowOrd, low_matrix,newton_rhs,lumped_rhs,NEWTON_TOL, NEWTON_MAX_ITER, verbose)) error("Newton's iteration failed.");	
	// 	Solution::vector_to_solution(u_L, &space, &low_sln);	



	//---------------------------------------antidiffusive fluxes-----------------------------------
		
	 antidiffusiveFlux(mass_matrix,lumped_matrix,conv_matrix,diffusion,vec_rhs, u_L, flux_scalar);		
	

	//-----------------corrected solution-----------------------------------

	/*vec_rhs->zero();
	lumped_matrix->multiply_with_vector(u_L,fct_scalar);
	vec_rhs->add_vector(fct_scalar); vec_rhs->add_vector(flux_scalar);	 
	for(int i=0; i<ndof;i++) coeff_vec[i]=0.0;		
	if(fct->solve()){ 
	coeff_vec =fct->get_solution();	
	 Solution::vector_to_solution(coeff_vec, &space, &u_prev_time);	
	  }else error ("Matrix solver failed.\n");	
	*/
    	
	//vielleicht geht das so schneller??????
	for(int i= 0; i<ndof; i++) coeff_vec[i]= u_L[i]+ (flux_scalar[i]/lumped_matrix->get(i,i));
	if (ts > 1 && ts % 10 == 0) Solution::vector_to_solution(coeff_vec, &space, &u_prev_time);
	


if (ts > 1 && ts % 10 == 0) {
	  	// Visualize the solution.	
		  sprintf(title, "high_Ord_Hermes Time %3.2f", current_time);
		  Hermview.set_title(title);
		  Hermview.show(&u_prev_time_high);
  // Visualize the solution.
	  sprintf(title, "low_Ord Time %3.2f", current_time);
	  Lowview.set_title(title);
	  Lowview.show(&low_sln);
	  // Visualize the solution.	 
	  sprintf(title, "korrigierte Loesung: Time %3.2f", current_time);
	  sview.set_title(title);
	  sview.show(&u_prev_time);
}

	  // Update global time.
	  current_time += time_step;

	  // Increase time step counter
	  ts++;

    info("Ende Zeitschritt: %i", ts-1);


}
while (current_time < T_FINAL);

  // Clean up.
  delete mass_matrix;  
  delete lumped_matrix; 
  delete conv_matrix;
  delete diffusion;
  delete low_matrix;
  delete K_D;
	delete matrix;
	delete solver;
	delete rhs;
	delete[] coeff_vec_newton;  

  delete[] lumped_scalar;  
  delete[] coeff_vec;  
 delete[] flux_scalar; 
//delete[] fct_scalar;

  delete vec_rhs;

	 delete lowOrd;
	//delete fct;

delete newton_rhs;

 

  // Wait for the view to be closed.
  View::wait();
  return 0;
}

