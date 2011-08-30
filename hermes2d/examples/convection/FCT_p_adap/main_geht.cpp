#define HERMES_REPORT_ALL
#include "hermes2d.h"
 #define PI (3.141592653589793) 

// 1. Step: (M_L/tau -theta(K+D)) u^L =   (M_L/tau + (1-theta)(K+D)) u^n
// 2. Step : f_ij = (M_c)_ij (dt_u_L(i)- dt_u_L(j)) + D_ij (u_L(i)- u_L(j)); f_i = sum_(j!=i) alpha_ij f_ij
// 3. Step:  M_L u^(n+1) = M_L u^L + tau * f 

const int INIT_REF_NUM =2;                   // Number of initial refinements.
const int P_INIT = 2;                             // Initial polynomial degree.
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
#include "extrapolation.cpp"

//Mass lumping an den Stellen von AsmList, sonst standard Massmatrix
UMFPackMatrix* massLumping(AsmList* al,UMFPackMatrix* mass_matrix){  //al=NULL=>lumped=mass
	 UMFPackMatrix* lumped_matrix = new UMFPackMatrix();   //M_L
	int size = mass_matrix->get_size();
	int nnz = mass_matrix->get_nnz();
	lumped_matrix->create(size, nnz, mass_matrix->get_Ap(), mass_matrix->get_Ai(),mass_matrix->get_Ax());
	scalar a =0.0;
	if(al!=NULL){  //nur wenn Liste fuer Vertex-DOFS vorhanden
	for(unsigned int i = 0; i < al->cnt; i ++){
		for(unsigned int j = (i+1); j < al->cnt; j ++){	//Mc_ij = Mc_ji				
			if(al->dof[i]!=al->dof[j]){  //kein Diagonaleintrag!
				if((al->dof[i]>=size)||(al->dof[j]>=size)) error("massLumping:wrong DOF");				
				a = lumped_matrix->get(al->dof[i],al->dof[j]);
				if(a!=0.0){
					lumped_matrix->add(al->dof[i],al->dof[i],a);    //zur Diagonale hinzufuegen
					lumped_matrix->add(al->dof[j],al->dof[j],a);    //zur Diagonale hinzufuegen
					lumped_matrix->add(al->dof[i],al->dof[j],-a);	//i,j Eintrag auf 0 setzen
					lumped_matrix->add(al->dof[j],al->dof[i],-a);	//i,j Eintrag auf 0 setzen
				}				
			}
		}			
	}
	}
	
	return lumped_matrix;
}

//Mass lumping an den Stellen von AsmList, sonst 0
UMFPackMatrix* pure_p1_massLumping(AsmList* al,UMFPackMatrix* mass_matrix){  //al=NULL=>lumped=mass
	 UMFPackMatrix* lumped_matrix = new UMFPackMatrix();   //M_L
	int size = mass_matrix->get_size();
	int nnz = mass_matrix->get_nnz();
	lumped_matrix->create(size, nnz, mass_matrix->get_Ap(), mass_matrix->get_Ai(),mass_matrix->get_Ax());
	lumped_matrix->zero();
	scalar a =0.0;
	if(al!=NULL){  //nur wenn Liste fuer Vertex-DOFS vorhanden
	for(unsigned int i = 0; i < al->cnt; i ++){
		for(unsigned int j = i; j < al->cnt; j ++){	//Mc_ij = Mc_ji				
				if((al->dof[i]>=size)||(al->dof[j]>=size)) error("massLumping:wrong DOF");				
				a = mass_matrix->get(al->dof[i],al->dof[j]);
				if(a!=0.0){
					lumped_matrix->add(al->dof[i],al->dof[i],a);    //zur Diagonale hinzufuegen
					if(al->dof[i]!=al->dof[j]) lumped_matrix->add(al->dof[j],al->dof[j],a);    //zur Diagonale hinzufuegen
				}				
			
		}			
	}
	}
	
	return lumped_matrix;
}


//artificial Diffusion
UMFPackMatrix* artificialDiffusion(AsmList* al,UMFPackMatrix* conv_matrix){ //al=NULL => diffusion=0
	 int size = conv_matrix->get_size();
	 int nnz = conv_matrix->get_nnz();
	scalar a,b;
	 UMFPackMatrix* diffusion = new UMFPackMatrix();  
	diffusion->create(size, nnz, conv_matrix->get_Ap(), conv_matrix->get_Ai(),conv_matrix->get_Ax());
	diffusion->zero();  //matrix = 0
	if(al!=NULL){ //nur wenn Liste fuer Vertex-DOFS vorhanden
	for(unsigned int i = 0; i < al->cnt; i ++){
	 for(unsigned int j = (i+1); j < al->cnt; j ++){
		if(al->dof[i]!=al->dof[j]){		
				       a= -conv_matrix->get(al->dof[i],al->dof[j]);
				      	b= -conv_matrix->get(al->dof[j],al->dof[i]);     
				     if((a>=b)&&(a>0.0)){
					 diffusion->add(al->dof[j],al->dof[i],a);
					diffusion->add(al->dof[i],al->dof[j],a);	
					 diffusion->add(al->dof[j],al->dof[j],-a);
					diffusion->add(al->dof[i],al->dof[i],-a);	
				     }else if((b>a)&&(b>0.0)){
					 diffusion->add(al->dof[i],al->dof[j],b);
					diffusion->add(al->dof[j],al->dof[i],b);
				 	 diffusion->add(al->dof[j],al->dof[j],-b);
					diffusion->add(al->dof[i],al->dof[i],-b);
	     				}
			}
		}
	}
	}

	return diffusion;

}
//Assemble antidiffusive fluxes & Limiter
//f_ij und alpha_ij werden nicht explizit berechnet!! da scalar** flux = new_matrix<scalar>(ndof,ndof); zuviel Speicher braucht
void antidiffusiveFlux(AsmList* al,UMFPackMatrix* mass_matrix,UMFPackMatrix* lumped_matrix,UMFPackMatrix* conv_matrix,UMFPackMatrix* diffusion,UMFPackVector* flux_dt_rhs, scalar* u_L, scalar* flux_scalar ){ //al==NULL =>flux=0
	int ndof = conv_matrix->get_size();
	scalar flux_dt_scalar[ndof];	
	Solver* flux_dt;
	scalar* dt_u_L = NULL;
	scalar P_plus[ndof], P_minus[ndof];	
	scalar Q_plus[ndof], Q_minus[ndof];	
	scalar R_plus[ndof], R_minus[ndof];	
	scalar alpha,f, plus, minus;
	for(int i=0; i<ndof;i++) flux_scalar[i]=0.0;
if(al!=NULL){
	conv_matrix->multiply_with_vector(u_L, flux_dt_scalar);
	flux_dt_rhs->zero(); flux_dt_rhs->add_vector(flux_dt_scalar);  //K u^L	
	flux_dt = create_linear_solver(matrix_solver,mass_matrix,flux_dt_rhs); //M_c u_t = K u^L
	if(flux_dt->solve())	dt_u_L = flux_dt->get_solution();	
	  else error ("Matrix solver failed.\n");

	//Berechnung von P&Q
	for(int i=0; i<ndof;i++){ P_plus[i]=0.0;P_minus[i]=0.0;Q_plus[i]=0.0;Q_minus[i]=0.0;}
	for(unsigned int i = 0; i < al->cnt; i ++){
	 	for(unsigned int j = (i+1); j < al->cnt; j ++){
			if(al->dof[i]!=al->dof[j]){		
				f = mass_matrix->get(al->dof[i],al->dof[j])*(dt_u_L[al->dof[i]]- dt_u_L[al->dof[j]]) 
					+ diffusion->get(al->dof[i],al->dof[j])*(u_L[al->dof[i]]- u_L[al->dof[j]]);
				if( (f*(u_L[al->dof[j]]- u_L[al->dof[i]])) > 0.0) f = 0.0; //prelimiting step
				if(f>0.0)	{ 
					P_plus[al->dof[i]]+=f;
					P_minus[al->dof[j]]-=f;
				}else if (f<0.0){
				 	P_minus[al->dof[i]]+=f;
					P_plus[al->dof[j]]-=f;
				}
				f = lumped_matrix->get(al->dof[i],al->dof[i])*(u_L[al->dof[j]]-u_L[al->dof[i]])/time_step; 
				if(f>Q_plus[al->dof[i]]) Q_plus[al->dof[i]] = f;				
				if(f<Q_minus[al->dof[i]]) Q_minus[al->dof[i]] = f;			
				f= lumped_matrix->get(al->dof[j],al->dof[j])*(u_L[al->dof[i]]-u_L[al->dof[j]])/time_step; 
				if(f>Q_plus[al->dof[j]]) Q_plus[al->dof[j]] = f;	
				if(f<Q_minus[al->dof[j]]) Q_minus[al->dof[j]] = f;
			}
		}
	}

	//Berechnung von R	
	for(unsigned int i = 0; i < al->cnt; i ++){
		plus = 1.0; minus = 1.0;		
		if(P_plus[al->dof[i]]!=0.0)  plus = Q_plus[al->dof[i]]/P_plus[al->dof[i]];		
		if(P_minus[al->dof[i]]!=0.0) minus = Q_minus[al->dof[i]]/P_minus[al->dof[i]];			
		if(plus>=1.0) R_plus[al->dof[i]]= 1.0;
		else 	     R_plus[al->dof[i]]= plus;
		if(minus>=1.0) R_minus[al->dof[i]]= 1.0;
		else 	     R_minus[al->dof[i]]= minus;	
	
	}
	
	//Berechnung von alpha & f_i
	alpha = 1.0;
	for(unsigned int i = 0; i < al->cnt; i ++){
	 	for(unsigned int j = (i+1); j < al->cnt; j ++){	
		 if(al->dof[i]!=al->dof[j]){			
			f= mass_matrix->get(al->dof[i],al->dof[j])*(dt_u_L[al->dof[i]]- dt_u_L[al->dof[j]]) 
				+ diffusion->get(al->dof[i],al->dof[j])*(u_L[al->dof[i]]- u_L[al->dof[j]]);	
			if( (f*(u_L[al->dof[j]]- u_L[al->dof[i]])) > 0.0) f = 0.0; //prelimiting step
			if(f>0.0){					
				if(R_plus[al->dof[i]]>R_minus[al->dof[j]]) alpha = R_minus[al->dof[j]];
				else 	alpha = R_plus[al->dof[i]];
			}else if (f<0.0){
				if(R_minus[al->dof[i]]>R_plus[al->dof[j]]) alpha = R_plus[al->dof[j]];
				else 	alpha = R_minus[al->dof[i]]; 
			}
			flux_scalar[al->dof[i]] += alpha*f;
			flux_scalar[al->dof[j]] -= alpha*f;
		  }				
		}
	}

	//Cleanup	
	delete flux_dt;
}	

}




void lumped_flux_limiter(AsmList* al,UMFPackMatrix* mass_matrix,UMFPackMatrix* lumped_matrix, scalar* u_L, scalar* u_H){
	int ndof = mass_matrix->get_size();
	scalar P_plus[ndof], P_minus[ndof];	
	scalar Q_plus[ndof], Q_minus[ndof];	
	scalar R_plus[ndof], R_minus[ndof];	
	scalar alpha,f, plus, minus;
	if(al!=NULL){
		//Berechnung von P&Q
		for(int i=0; i<ndof;i++){ P_plus[i]=0.0;P_minus[i]=0.0;Q_plus[i]=0.0;Q_minus[i]=0.0;}
		for(unsigned int i = 0; i < al->cnt; i ++){
		 	for(unsigned int j = (i+1); j < al->cnt; j ++){
				if(al->dof[i]!=al->dof[j]){		
					f = mass_matrix->get(al->dof[i],al->dof[j])*(u_H[al->dof[i]]- u_H[al->dof[j]]);
					if( (f*(u_H[al->dof[j]]- u_H[al->dof[i]])) > 0.0) f = 0.0; //prelimiting step
					if(f>0.0)	{ 
						P_plus[al->dof[i]]+=f;
						P_minus[al->dof[j]]-=f;
					}else if (f<0.0){
					 	P_minus[al->dof[i]]+=f;
						P_plus[al->dof[j]]-=f;
					}				
					f = (u_L[al->dof[j]]-u_L[al->dof[i]]); 
					if(f>Q_plus[al->dof[i]]) Q_plus[al->dof[i]] = f;				
					if(f<Q_minus[al->dof[i]]) Q_minus[al->dof[i]] = f;
					f= (u_L[al->dof[i]]-u_L[al->dof[j]]);
					if(f>Q_plus[al->dof[j]]) Q_plus[al->dof[j]] = f;	
					if(f<Q_minus[al->dof[j]]) Q_minus[al->dof[j]] = f;
				}
			}
		}
	
		//Berechnung von R
	for(unsigned int i = 0; i < al->cnt; i ++){
		plus = 1.0; minus = 1.0;		
		if(P_plus[al->dof[i]]!=0.0)  plus = lumped_matrix->get(al->dof[i],al->dof[i])*Q_plus[al->dof[i]]/P_plus[al->dof[i]];		
		if(P_minus[al->dof[i]]!=0.0) minus = lumped_matrix->get(al->dof[i],al->dof[i])*Q_minus[al->dof[i]]/P_minus[al->dof[i]];			
		if(plus>=1.0) R_plus[al->dof[i]]= 1.0;
		else 	     R_plus[al->dof[i]]= plus;
		if(minus>=1.0) R_minus[al->dof[i]]= 1.0;
		else 	     R_minus[al->dof[i]]= minus;	
	
	}	
		//Berechnung von alpha & f_i
		alpha = 1.0;
		for(unsigned int i = 0; i < al->cnt; i ++){
		 	for(unsigned int j = (i+1); j < al->cnt; j ++){	
			 if(al->dof[i]!=al->dof[j]){			
				f= mass_matrix->get(al->dof[i],al->dof[j])*(u_H[al->dof[i]]- u_H[al->dof[j]]) ;	
				if( (f*(u_H[al->dof[j]]- u_H[al->dof[i]])) > 0.0) f = 0.0; //prelimiting step				
				if(f>0.0){					
					if(R_plus[al->dof[i]]>R_minus[al->dof[j]]) alpha = R_minus[al->dof[j]];
					else 	alpha = R_plus[al->dof[i]];
				}else{
					if(R_minus[al->dof[i]]>R_plus[al->dof[j]]) alpha = R_plus[al->dof[j]];
					else 	alpha = R_minus[al->dof[i]]; 
				}
				u_L[al->dof[i]] += alpha*f/(lumped_matrix->get(al->dof[i],al->dof[i]));
				u_L[al->dof[j]] -= alpha*f/(lumped_matrix->get(al->dof[j],al->dof[j]));
			  }				
			}
		}
	
		
		
	}



}

//Durchlaufen und p bestimmen: in dof_list dof von vertex-basisfunction speichern 
//Quad hat p-Ordnung fuer horizontal und vertical zusammen space->get_element_order(e->id) = H2D_MAKE_QUAD_ORDER(h-ord, v-ord)
//H2D_GET_H_ORDER(id),  H2D_GET_V_ORDER(id)
void p1_list(H1Space* space, AsmList* dof_list, Element* e,AsmList* al ){
	bool more = false;
	for_all_active_elements(e, space->get_mesh()){
	 //if( space->get_element_order(e->id)== H2D_MAKE_QUAD_ORDER(1, 1)){	//Ordnung soll 1 
	space->get_element_assembly_list(e, al);
	  	for (unsigned int iv = 0; iv < e->nvert; iv++){   		
		  int index =  space->get_shapeset()->get_vertex_index(iv);
			Node* vn = e->vn[iv];
		  if (space->get_element_order(e->id) == 0) break;
		  if (!vn->is_constrained_vertex()){  //unconstrained ->kein haengender Knoten!!!
			for(unsigned int j = 0; j < al->cnt; j ++){			 
		    		if((al->idx[j]==index)&&(al->dof[j]!=-1.0)){ 
					for(unsigned int i = 0; i < dof_list->cnt; i ++){ 
						if(dof_list->dof[i]==al->dof[j]){ more =true; break;}	//ueberpruefen ob dof schon in liste enhalten
					}
				if(more==false) dof_list->add_triplet(index, al->dof[j], 1.0);  //dof=-1 =>dirichlet
				more = false;
				}
			}
		   }
		}
	//}	

	 //printf("OrdV=%d  : ", H2D_GET_V_ORDER(space.get_element_order(e->id)));
	//printf("OrdH=%d \n", H2D_GET_H_ORDER(space.get_element_order(e->id)));
	}

}



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
  Solution low_sln;

  // Previous time level solution (initialized by the initial condition).
  CustomInitialCondition u_prev_time(&mesh);  

  // Initialize the weak formulation.
CustomWeakFormMassmatrix massmatrix(time_step, &u_prev_time);
CustomWeakFormConvection convection(&u_prev_time);

  // Initialize the FE problem.
  DiscreteProblem dp_mass(&massmatrix, &space);
  DiscreteProblem dp_convection(&convection, &space);

	//Durchlaufen und p bestimmen: in dof_list dof von vertex-basisfunction speichern 
	Element* e = NULL;  
	AsmList* al = new AsmList();	
	AsmList* dof_list = new AsmList();
	p1_list(&space, dof_list, e,al );
	//for(unsigned int j = 0; j < dof_list->cnt; j ++)printf("DOFS: %d \n",dof_list->dof[j]);			
		printf("p1-DOFS: %d \n",dof_list->cnt);	
	//dof_list = NULL;
//--------------------------------------------------------------------------------------------------------------------------------------------

Node* vn; Node* edge;
scalar* coeff_vec = new scalar[ndof];
for(int i=0; i<ndof;i++) coeff_vec[i]=0.0;
//OGProjection::project_global(&space,&u_prev_time, coeff_vec, matrix_solver, HERMES_L2_NORM);
//Solution::vector_to_solution(coeff_vec, &space, &u_prev_time);
//initview.show(&u_prev_time);
MeshView mview("Hello world!", new WinGeom(700, 700, 350, 350));
ScalarView exview("extrapolierte Loesung", new WinGeom(0, 0, 500, 400));
ScalarView initview("Solution", new WinGeom(700, 0, 500, 400));

int counter = 0;


for_all_active_elements(e, space.get_mesh()){// e ist Ausgangselement
	
	int v_ord = 0.0; int h_ord =0.0; int index =0.0; 
	int max_ord = 0.0;   //Polynomgrad in Element e, bei Quad: max(v_ord,h_ord)	
	int no_of_tris = 0;
	int no_of_quads = 0;


//----------------------
//Polynomgrad bestimmen
	if (space.get_element_order(e->id) == 0.0) break;    
	if(e->is_triangle()==true){ //Triangle
		max_ord = space.get_element_order(e->id);			
		no_of_tris++;	
	}else if(e->is_quad()==true){ //Quad
		v_ord = H2D_GET_V_ORDER(space.get_element_order(e->id));
		h_ord = H2D_GET_H_ORDER(space.get_element_order(e->id));
		if(v_ord>= h_ord) max_ord = v_ord;
		else max_ord = h_ord;  //max Polynomgrad im Quad		
		no_of_quads++;		
	}else break;

//Liste zusammenzustellen fuer vertex-dofs von e	
	AsmList* e_list = new AsmList();  //Liste von vertex-dofs von e 	
	double* x_coord = new double[e->nvert];
	double* y_coord= new double[e->nvert]; //Koordinaten des Elements (jedes Knotenpunkts)
	int* number_of_elem= new int[e->nvert];  //Anzahl der Elemente, die Knoten benutzen
	space.get_element_assembly_list(e, al);
	for (unsigned int iv = 0; iv < e->nvert; iv++){  
		index =  space.get_shapeset()->get_vertex_index(iv);
		vn = e->vn[iv];	
		if(vn->ref>=TOP_LEVEL_REF) number_of_elem[iv] = vn->ref- TOP_LEVEL_REF; 	//keine Ahnung, ob das immer funzt.!!..			
		else 	number_of_elem[iv] = vn->ref;	
		x_coord[iv] = vn->x;
		y_coord[iv] = vn->y;
		//printf("(x,y):(%f,%f), %i \n",x_coord[iv],y_coord[iv],number_of_elem[iv] );
		for(unsigned int j = 0; j < al->cnt; j ++)			 
	    		if(al->idx[j]==index) e_list->add_triplet(index, al->dof[j], 1.0);  //dof=-1 =>dirichlet
			
	} 

 	//printf("nur e : size: e_list %i, quads: %i, tris: %i\n", e_list->cnt,no_of_quads,no_of_tris );

	
// nun Patch (=Vereinigung aller Nachbar-elemente) bestimmen
	AsmList* ne_list = new AsmList();  //speichert in index = Elementid
	bool in_list = false;
	Element* ne = NULL;
	int no_of_nb= 0;
	for(int i = 0; i<e->nvert; i++){			
		 if(number_of_elem[i]>1){ //error("Element vergessen");  
			vn = e->vn[i];						
			for_all_active_elements(ne, space.get_mesh()){
				if(ne!=e){
					for(int c=0; c < ne->nvert; c++){
						in_list = false;
						if(ne->vn[c]==vn){
							for(unsigned int l = 0;l<ne_list->cnt; l++) //ueberpruefen ob schon in Liste
								if(ne->id== ne_list->idx[l]){in_list = true; break;}//schon in liste
							if(in_list==true) break;
							else{//in Liste noch aufnehmen							
								for(int k=0; k< ne->nvert; k++)//ueberpruefen ob wirklich Nachbar
									for(int j = 0; j<e->nvert; j ++)
										if(ne->vn[k]==e->vn[j]) number_of_elem[j]--;
								if(ne->is_triangle()==true)no_of_tris++;
								else if(ne->is_quad()==true)no_of_quads++; 		
								ne_list->add_triplet(ne->id, no_of_nb, 1.0);
								no_of_nb++;
							}
						  break;	
						}

							
					}
										
				}
				if(number_of_elem[i]==1) break;
				
			}
		}else if(number_of_elem[i]<1) error("zu oft durchlaufen");
	}
			


	//Patch speichern als mesh
	int nv = 4* no_of_quads + 3* no_of_tris;	//erstmal zu gross schaetzen :)
	double2* verts = new double2[nv];
	int4* tris = new int4[no_of_tris];
	int5* quads = new int5[no_of_quads];
	int nm = 0;
	int3* mark = new int3[nv];
	int ver = 0;
	int tri = 0;
	int qua= 0;
	
	
	//erst e hinzufuegen
	for (unsigned int iv = 0; iv < e->nvert; iv++){ 			
		vn = e->vn[iv];	
		edge = e->en[iv];				
		verts[ver][0] = vn->x;
		verts[ver][1] = vn->y;
		if(e->is_triangle()==true) tris[tri][iv]= ver;
		else if(e->is_quad()==true) quads[qua][iv]=ver;	
		if(edge->bnd==1){   //boundary edge
			mark[nm][0]= ver; 
			if(iv== ((e->nvert)-1.0)) mark[nm][1]=0;
			else mark[nm][1]= ver+1; 
			mark[nm][2] = edge->marker;				
			 nm++;
		} 		
		ver++;	
	}
	if(ver < e->nvert) error("Anzahl der Vertices kleiner als e->nvert");

	if(e->is_triangle()==true){ tris[tri][3]= 1;tri++;}  //Element markieren! mark = 1 -> e, sonst mark =HERMES_ANY_INT
	else if(e->is_quad()==true){ quads[qua][4]=1;qua++;}

	bool p2 =true; //von Randkante 2.Element bestimmt
	for(unsigned int j = 0; j <ne_list->cnt; j++){  //Alle Nachbarelemente durchlaufen
			ne = space.get_mesh()->get_element(ne_list->idx[j]);
			for (int iv = 0; iv < ne->nvert; iv++){ 		//Alle Knoten von ne durchlaufen	
				vn = ne->vn[iv];
				in_list = false;
				for(int k =0; k<ver; k ++){  //ueberpruefen ob Knoten schon in Liste enthalten
					if((vn->x==verts[k][0])&&(vn->y==verts[k][1])){  //Knoten schon enthalten
						in_list =true;
						if(ne->is_triangle()==true) tris[tri][iv]= k;   //in Element eintragen
						else if(ne->is_quad()==true) quads[qua][iv]=k;						
						if((k<e->nvert)&&(e->en[k]->bnd!=1)&&(e->en[e->prev_vert(k)]->bnd!=1)){// Knoten von e und kein Rand!							
									p2 = true;	
						} else{
							if(p2==false) { //erst noch p2 setzen
								mark[nm][1] = k;
								nm++;
							}						
							mark[nm][0] = k;   //p1 setzen
							mark[nm][2] = 2; //Marker setzen
							p2 = false;	
							if(iv == (ne->nvert-1)){ //letzter Knoten des Elements erreicht! p2 evtl vn[0] 
								vn = ne->vn[0];
								for(int k =0; k<ver; k ++){  //bestimmen von k
										if((vn->x==verts[k][0])&&(vn->y==verts[k][1])){  //Knoten schon enthalten
											in_list = true;
											if((k<e->nvert)&&(e->en[k]->bnd!=1)&&(e->en[e->prev_vert(k)]->bnd!=1)){// Knoten von e und kein Randknoten!!
												 //do nothing
											} else{	mark[nm][1] = k;
													nm++;
											}
											p2 = true;									
									 		break; //Schleife nicht weiter durchlaufen!
									 	} 
								}
								if(in_list == false) error("p2 konnte nicht gefunden werden!");							
							}
						}				
					 	break;  //Schleife nicht weiter durchlaufen!
					 } 
				}
				if(in_list==false) { //neuer Knoten
					verts[ver][0] = vn->x;  //in Knotenliste eintragen
					verts[ver][1] = vn->y;
					if(p2 == false){ mark[nm][1] = ver; nm++;}
					mark[nm][0] = ver;   //p1 setzen
					mark[nm][2] = 2; //Marker setzen
					p2 = false;
					if(ne->is_triangle()==true) tris[tri][iv]= ver;   //in Element eintragen
					else if(ne->is_quad()==true) quads[qua][iv]=ver;
					if(iv == (ne->nvert-1)){ //letzter Knoten des Elements erreicht! p2 evtl vn[0] 
						vn = ne->vn[0];
						for(int k =0; k<ver; k ++){  //bestimmen von k
								if((vn->x==verts[k][0])&&(vn->y==verts[k][1])){  //Knoten schon enthalten
									in_list = true;
									if((k<e->nvert)&&(e->en[k]->bnd!=1)&&(e->en[e->prev_vert(k)]->bnd!=1)){// Knoten von e und kein Randknoten!!
										 //do nothing
									} else{	mark[nm][1] = k;
											nm++;
									}
									p2 = true;									
							 		break;
							 	} 
						}
						if(in_list == false) error("p2 konnte nicht gefunden werden!");
					}	
					ver++;
				}
			}
		//komplettes Element durchlaufen, Marker setzen
		if(ne->is_triangle()==true){ tris[tri][3]= HERMES_ANY_INT;tri++; } 
		else if(ne->is_quad()==true){ quads[qua][4]=HERMES_ANY_INT;qua++; }
		if(p2==false) error("p2 ist nicht gesetzt, obwohl Element schon durchlaufen");
	}

	nv = ver;	

	Mesh* patch_mesh = new Mesh();
	patch_mesh->create(nv, verts, tri, tris, qua, quads,nm,mark );  //Patch erstellen

	H1Space* patch_space = new H1Space(patch_mesh, P_INIT);
	for_all_active_elements(ne, patch_space->get_mesh()){
		if(ne->marker == 1) e = ne; 		//Element e markiert mit mark = 1
		
	}	
//Patch auf richtige Ordnung setzen (Achtung: Ordnung wird an Ordnung von e angepasst!!!)
	patch_space->set_uniform_order(max_ord);  //setzt alle Elemente auf Ordnung max_ord
	if(e->is_quad()&&(v_ord!=h_ord)){  //braucht man nur wenn e ein Quad und vertikale und horizontale Ordnung sich unterscheiden
		patch_space->set_uniform_order_internal(Ord2(h_ord, v_ord), 1);  //setzt e (marker =1) auf v_ord, h_ord
		patch_space->assign_dofs();
	}	

	CustomInitialCondition* u_init= new CustomInitialCondition(patch_mesh);
	ExtrapolationSolution* u_ex = new ExtrapolationSolution(patch_mesh, e);

	u_ex->set_extra_fn(u_init);
	u_ex->init_extra_fn(patch_space);

	ErrorEstimation* err = new ErrorEstimation(patch_space, HERMES_L2_NORM);
	exview.set_min_max_range(-0.05, 0.1); 
	initview.set_min_max_range(-0.05, 0.1); 
	exview.show(u_ex);
	//exview.show(u_init);
	//mview.show(patch_mesh);

	OGProjection::project_global(patch_space,u_init, coeff_vec, matrix_solver, HERMES_L2_NORM);
	Solution::vector_to_solution(coeff_vec, patch_space, u_init);
	initview.show(u_init);

	printf("error(%i): %f\n",counter, err->calc_err_extrapolation(u_init, u_ex, false, HERMES_TOTAL_ERROR_ABS));
	
	//if (counter==14)
	// View::wait();

	delete  e_list;
	delete ne_list;
	delete [] x_coord ;
	delete [] y_coord;
	
	delete patch_mesh;
	delete patch_space;
	delete u_init;
	delete u_ex;
	delete err;
	delete [] number_of_elem;
	delete [] verts;
	delete [] tris;
	delete [] quads;

 counter++;
}






//-----------------------------------------------------------------------------------------------------------------------------------------
 	 // Previous time level solution (initialized by the initial condition).
 		 CustomInitialCondition u_prev_time_high(&mesh);
		ConvectionForm high_convection(time_step, &u_prev_time_high);
		  // Instantiate a class with global functions.
		  Hermes2D hermes2d;
		  // Project the initial condition on the FE space to obtain initial
		  // coefficient vector for the Newton's method.		
		  scalar* coeff_vec_newton = new scalar[ndof];
		  OGProjection::project_global(&space, &u_prev_time_high, coeff_vec_newton, matrix_solver, HERMES_L2_NORM);

		
//Solution::vector_to_solution(coeff_vec_newton, &space, &u_prev_time_high);

		  // Initialize the FE problem.
		  DiscreteProblem dp(&high_convection, &space);
		  // Set up the solver, matrix, and rhs according to the solver selection.
		  SparseMatrix* matrix = create_matrix(matrix_solver);
		  Vector* rhs = create_vector(matrix_solver);
		  Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);
		   ScalarView Hermview("hohe Ordnung", new WinGeom(0, 0, 500, 400));
		  //Hermview.show(&u_prev_time_high);
		  bool jacobian_changed = true;
		bool verbose = true;
		


//----------------------MassLumping M_L/tau--------------------------------------------------------------------

  // Set up the matrix, and rhs according to the solver selection.=>For Masslumping
  UMFPackMatrix* mass_matrix = new UMFPackMatrix();   //M_c/tau
	dp_mass.assemble(mass_matrix, NULL); 	
UMFPackMatrix* lumped_matrix = massLumping(dof_list,mass_matrix);


//------------------------artificial DIFFUSION D---------------------------------------

  // Set up the solver, matrix, and rhs according to the solver selection.=>artificial Diffusion
  UMFPackMatrix* conv_matrix = new UMFPackMatrix();   //K
	dp_convection.assemble(conv_matrix, NULL,true);
 UMFPackMatrix* diffusion = artificialDiffusion(dof_list,conv_matrix);

//--------------------------------------------------------------------------------------------
 UMFPackMatrix* low_matrix = new UMFPackMatrix();  
 UMFPackMatrix* K_D = new UMFPackMatrix();  
K_D->create(conv_matrix->get_size(),conv_matrix->get_nnz(), conv_matrix->get_Ap(), conv_matrix->get_Ai(),conv_matrix->get_Ax());
K_D->add_matrix(diffusion); 

//K_D->zero();


low_matrix->create(K_D->get_size(),K_D->get_nnz(), K_D->get_Ap(), K_D->get_Ai(),K_D->get_Ax());
//(-theta)(K+D)
if(theta==0) low_matrix->zero();
else	low_matrix->multiply_with_scalar(-theta);
//(1-theta)(K+D)
if(theta ==1) K_D->zero();
else K_D->multiply_with_scalar((1.0-theta));


//M_L/tau - theta(D+K)
low_matrix->add_matrix(lumped_matrix);  //kann nun fuer alle Zeitschritte verwendet werden (ausser bei Adaptivitaet)
//low_matrix->add_matrix(mass_matrix); 
//M_L/tau+(1-theta)(K+D)
K_D->add_matrix(lumped_matrix);
//K_D->add_matrix(mass_matrix);






/*

  // Initialize views.
   ScalarView Lowview("niedriger Ordnung", new WinGeom(500, 500, 500, 400));
 // Lowview.show(&u_prev_time, HERMES_EPS_HIGH);
  ScalarView sview("Solution", new WinGeom(0, 500, 500, 400));
  //sview.show(&u_prev_time, HERMES_EPS_HIGH); 

		




// Time stepping loop:
double current_time = 0.0; 
int ts = 1;

scalar* lumped_scalar = new scalar[ndof];
UMFPackVector*  vec_rhs = new UMFPackVector(ndof);
//scalar* coeff_vec = new scalar[ndof];
//for(int i=0; i<ndof;i++) coeff_vec[i]=0.0;
scalar* coeff_vec_2 = new scalar[ndof];
for(int i=0; i<ndof;i++) coeff_vec_2[i]=0.0;

scalar* u_L = NULL;

scalar* flux_scalar = new scalar[ndof];; 
Solver* lowOrd;

char title[100];


	lumped_matrix->multiply_with_scalar(time_step);  // M_L
	mass_matrix->multiply_with_scalar(time_step);  // massmatrix = M_C

// Project the initial condition on the FE space->coeff_vec		
	//Lumped_Projection::project_lumped(&space, &u_prev_time, coeff_vec, matrix_solver, lumped_matrix);
	//Lumped_Projection::project_lumped(&space, &u_prev_time, coeff_vec, matrix_solver, NULL);

UMFPackMatrix* proj_matrix = pure_p1_massLumping(dof_list,mass_matrix);

			Lumped_Projection::project_lumped_rhs(&space, &u_prev_time, coeff_vec, matrix_solver, proj_matrix);
			//Lumped_Projection::project_lumped(&space, &u_prev_time, coeff_vec, matrix_solver, lumped_matrix);
			//OGProjection::project_global(&space,&u_prev_time, coeff_vec_2, matrix_solver, HERMES_L2_NORM);
			//lumped_flux_limiter(dof_list,mass_matrix, lumped_matrix, coeff_vec, coeff_vec_2);
			//lumped_flux_limiter(dof_list,mass_matrix, proj_matrix, coeff_vec, coeff_vec_2);
		
	//OGProjection::project_global(&space,&u_prev_time, coeff_vec, matrix_solver, HERMES_L2_NORM); 



Solution::vector_to_solution(coeff_vec, &space, &u_prev_time);
ScalarView pview("projezierter Anfangswert", new WinGeom(500, 0, 500, 400));
	 // pview.show(&u_prev_time);



do
{	//if(ts==2) for(int i=0;i<ndof;i++) coeff_vec[i]= coeff_vec_newton[i];
	  info(" Time step %d, time %3.5f", ts, current_time);  
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
	lowOrd = create_linear_solver(matrix_solver,low_matrix,vec_rhs);	
	if(lowOrd->solve()){ 
		u_L = lowOrd->get_solution();  
		//coeff_vec = lowOrd->get_solution();  
		Solution::vector_to_solution(u_L, &space, &low_sln);	
	  }else error ("Matrix solver failed.\n");

	//---------------------------------------antidiffusive fluxes-----------------------------------		
	 antidiffusiveFlux(dof_list,mass_matrix,lumped_matrix,conv_matrix,diffusion,vec_rhs, u_L, flux_scalar);	
  	
	//vielleicht geht das so schneller??????
	for(int i= 0; i<ndof; i++) coeff_vec[i]= u_L[i]+ (flux_scalar[i]*time_step/lumped_matrix->get(i,i));
	 Solution::vector_to_solution(coeff_vec, &space, &u_prev_time);
	

	  	// Visualize the solution.	
		  sprintf(title, "high_Ord_Hermes Time %3.2f", current_time);
		  Hermview.set_title(title);
		  Hermview.show(&u_prev_time_high);

	
  // Visualize the solution.
	  sprintf(title, "low_Ord Time %3.2f", current_time);
	  Lowview.set_title(title);
	 // Lowview.show(&low_sln);
	  // Visualize the solution.	 
	  sprintf(title, "korrigierte Loesung: Time %3.2f", current_time);
	  sview.set_title(title);
	 // sview.show(&u_prev_time);


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
  delete[] lumped_scalar;  
  delete[] coeff_vec;  
	  delete[] coeff_vec_2;
 delete[] flux_scalar; 
  delete vec_rhs;
 delete lowOrd;
 delete proj_matrix;


delete dof_list;
delete al; 
 */

  // Wait for the view to be closed.
  View::wait();
  return 0;
}

