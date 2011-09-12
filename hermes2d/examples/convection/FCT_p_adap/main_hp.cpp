#define HERMES_REPORT_ALL
#include "hermes2d.h"
 #define PI (3.141592653589793)     


using namespace RefinementSelectors;

const int INIT_REF_NUM =4;                   // Number of initial refinements.
const int P_INIT = 3;                             // Initial polynomial degree.
double time_step = 1e-3;                           // Time step.
//const double T_FINAL = 2*PI;                       // Time interval length.
const double T_FINAL = 0.1;  

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
const CandList CAND_LIST = H2D_P_ISO;          // Predefined list of element refinement candidates. Possible values are
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
const double ERR_STOP = 1.0;                      // Stopping criterion for adaptivity (rel. error tolerance between the
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
#include "extrapolation.cpp"
#include "least_square.cpp"
#include "error_est.cpp"
#include "patch_solution.cpp"


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


void p_adap(Space* space, Solution* sln, scalar* sln_coeff){
	info("p-adap");

	if(sln->get_type() == HERMES_EXACT) info("exakte Loesung");

	ScalarView exview("extrapolierte Loesung", new WinGeom(0, 0, 500, 400));
	ScalarView lsview("LS- Loesung", new WinGeom(0, 700, 500, 400));
	ScalarView initview("Solution", new WinGeom(700, 0, 500, 400));
	Node* vn; Node* edge;
	AsmList* al = new AsmList();	
	int counter = 0;
	int v_ord = 0.0; int h_ord =0.0; int index =0.0; 
	int max_ord = 0.0;   //Polynomgrad in Element e, bei Quad: max(v_ord,h_ord)	
	int no_of_tris = 0;
	int no_of_quads = 0;
	Element* e =NULL;

	bool* elements_to_refine = new bool[space->get_mesh()->get_max_element_id()];

	for_all_active_elements(e, space->get_mesh()){// e ist Ausgangselement
	
	 v_ord = 0.0;  h_ord =0.0;  index =0.0; 
	 max_ord = 0.0;   //Polynomgrad in Element e, bei Quad: max(v_ord,h_ord)	
	 no_of_tris = 0;
	 no_of_quads = 0;
	int element_id = e->id;
	

//----------------------
//Polynomgrad bestimmen
	if (space->get_element_order(element_id) == 0.0) break;    
	if(e->is_triangle()==true){ //Triangle
		max_ord = space->get_element_order(element_id);			
		no_of_tris++;	
	}else if(e->is_quad()==true){ //Quad
		v_ord = H2D_GET_V_ORDER(space->get_element_order(element_id));
		h_ord = H2D_GET_H_ORDER(space->get_element_order(element_id));
		if(v_ord>= h_ord) max_ord = v_ord;
		else max_ord = h_ord;  //max Polynomgrad im Quad		
		no_of_quads++;		
	}else break;

//Liste zusammenzustellen fuer vertex-dofs von e	
	AsmList* e_list = new AsmList();  //Liste von vertex-dofs von e 	
	double* x_coord = new double[e->nvert];
	double* y_coord= new double[e->nvert]; //Koordinaten des Elements (jedes Knotenpunkts)
	int* number_of_elem= new int[e->nvert];  //Anzahl der Elemente, die Knoten benutzen
	space->get_element_assembly_list(e, al);
	for (unsigned int iv = 0; iv < e->nvert; iv++){  
		index =  space->get_shapeset()->get_vertex_index(iv);
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
			for_all_active_elements(ne, space->get_mesh()){
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
	
	int* elem_ids = new int[no_of_nb];  //id in space, id in patch_space
	int elems = 0;
	
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
	elem_ids[elems]=  e->id;
	elems++;

	bool p2 =true; //von Randkante 2.Element bestimmt
	for(unsigned int j = 0; j <ne_list->cnt; j++){  //Alle Nachbarelemente durchlaufen
			ne = space->get_mesh()->get_element(ne_list->idx[j]);
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
		if(ne->is_triangle()==true){ tris[tri][3]= HERMES_ANY_INT; tri++; } 
		else if(ne->is_quad()==true){ quads[qua][4]=HERMES_ANY_INT; qua++; }
		elem_ids[elems]=  ne->id;	
		elems++;
		if(p2==false) error("p2 ist nicht gesetzt, obwohl Element schon durchlaufen");
	}

	nv = ver;	

	Mesh* patch_mesh = new Mesh();
	patch_mesh->create(nv, verts, tri, tris, qua, quads,nm,mark );  //Patch erstellen

	H1Space* patch_space = new H1Space(patch_mesh, P_INIT);

	Element* e_patch =NULL;
	for_all_active_elements(ne, patch_space->get_mesh()){
		if(ne->marker == 1) e_patch = ne; 		//Element e markiert mit mark = 1
	}	
		
//Patch auf richtige Ordnung setzen (Achtung: Ordnung wird an Ordnung von e angepasst!!!)
	patch_space->set_uniform_order(max_ord);  //setzt alle Elemente auf Ordnung max_ord
	if(e_patch->is_quad()&&(v_ord!=h_ord)){  //braucht man nur wenn e ein Quad und vertikale und horizontale Ordnung sich unterscheiden
		patch_space->set_uniform_order_internal(Ord2(h_ord, v_ord), 1);  //setzt e (marker =1) auf v_ord, h_ord
		patch_space->assign_dofs();
	}	
	//printf("error in step %i \n",counter);
	
	//Solution* u_init;
	 CustomInitialCondition* u_init= new CustomInitialCondition(patch_mesh);
	if(sln->get_type() != HERMES_EXACT) {	
		delete u_init;
		PatchSolution* u_init = new PatchSolution();
		if(sln_coeff==NULL)	u_init->init(patch_space, space, sln, elem_ids, elems);
		else u_init->init(patch_space, space, sln, sln_coeff, elem_ids, elems);
	}


	ExtrapolationSolution* u_ex = new ExtrapolationSolution(patch_mesh, e_patch);

	u_ex->set_fn(u_init);
	u_ex->init_extrapolation(patch_space);		

	//OGProjection::project_global(patch_space,u_init, coeff_vec, matrix_solver, HERMES_L2_NORM);
	//Solution::vector_to_solution(coeff_vec, patch_space, u_init);	
		
	LSSolution* u_LS = new LSSolution(patch_mesh,e_patch);		
	u_LS->LS_project(u_init, patch_space,matrix_solver);


	//Views:
		exview.set_min_max_range(-0.1, 1); 
	initview.set_min_max_range(-0.1, 1); 
	lsview.set_min_max_range(-0.1, 1); 
	//initview.show(u_init);		
	//lsview.show(u_LS);
	//exview.show(u_ex);

	ErrorEstimation* err_ex = new ErrorEstimation(patch_space, HERMES_L2_NORM, e_patch);
	ErrorEstimation* err_ls = new ErrorEstimation(patch_space, HERMES_L2_NORM, e_patch);
	
	double error_ex = err_ex->calc_err_est(u_ex, u_init,true, HERMES_TOTAL_ERROR_ABS);
	double error_ls = err_ls->calc_err_est(u_LS, u_init,true, HERMES_TOTAL_ERROR_ABS);

	//printf("LS-error: %f ; Ex-error: %f \n ",error_ls,error_ex );	
	//	printf("error(LS/ex)= %f \n",error_ls/error_ex );
	
	double error_limit = 0.1;

	if((error_ls/error_ex < error_limit)&&(error_ex!=0.0)) elements_to_refine[element_id] = true;
	else elements_to_refine[element_id] = false;
	
	//if(err_ex->get_element_error_squared(0,e->id)>1e-5) 
	/*for_all_active_elements(ne, patch_space->get_mesh())
		printf("error_ex(%i): %f\n",ne->id, err_ex->get_element_error_squared(0,ne->id)); //Elementweiser Fehler
	*/

	/*for_all_active_elements(ne, patch_space->get_mesh())
		printf("error_ls(%i): %f\n",ne->id, err_LS->get_element_error_squared(0,ne->id)); 
	*/
	
	//printf("---------------------------- \n");
	//View::wait(HERMES_WAIT_KEYPRESS);

	delete  e_list;
	delete ne_list;
	delete [] x_coord ;
	delete [] y_coord;
	delete [] number_of_elem;
	delete [] verts;
	delete [] tris;
	delete [] quads;
	delete []  mark;
	delete [] elem_ids;
	
	delete patch_mesh;
	delete patch_space;
	delete u_init;
	delete u_LS;
	delete u_ex;
	
	
	delete err_ex;
	delete err_ls;


 counter++;
}

ErrorEstimation* adapt = new ErrorEstimation(space, HERMES_L2_NORM, e);
	adapt->adapt(elements_to_refine);

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
  Aview.show(&sln_prev_time);

  // Initialize the weak formulation.
CustomWeakFormMassmatrix massmatrix(time_step, &sln_prev_time);
CustomWeakFormConvection convection(&sln_prev_time);

  // Create a refinement selector.
  H1ProjBasedSelector selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);



  // Visualize initial condition.
  char title[100];
  ScalarView view("Initial condition", new WinGeom(0, 500, 500, 400));
  OrderView ordview("Initial mesh", new WinGeom(550, 0, 500, 400)); 
  ordview.show(&space);
  OrderView ord2view("ref mesh", new WinGeom(550, 500, 500, 400));
 

	p_adap(&space, &sln_prev_time, NULL);

// Time stepping loop:
double current_time = 0.0; int ts = 1;
int ref_dof;
scalar* u_L = NULL;

do
{
    // Periodic global derefinement. 
    if ((ts > 1 && ts % UNREF_FREQ == 0)||(Space::get_num_dofs(&space) >= NDOF_STOP)) 
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
    bool done = false; int as = 1;
    double err_est;
 info("Current Time: %f, Time step %d",current_time, ts);
    do {
	      //info("Time step %d, adaptivity step %d:", ts, as);

	     // Construct globally refined reference mesh and setup reference space.
		 Space* ref_space = Space::construct_refined_space(&space, 0);//verfeinert 1*in h & 0*in p 

	//Eigentlich will ich nur mit p-Verfeinerung testen, aber das funzt nicht->fehlermeldung:  	
  	/*Mesh* ref_mesh = new Mesh;
  	ref_mesh->copy(&basemesh);
	Space* ref_space = space.dup(ref_mesh, 2); //nur p-Verfeinerung	 */    	
		ref_dof = Space::get_num_dofs(ref_space);

info("adaptivity step %d: ndof: %d, ref_ndof: %d", as, Space::get_num_dofs(&space), Space::get_num_dofs(ref_space));




//Durchlaufen und p bestimmen:
//Quad hat p-Ordnung fuer horizontal und vertical zusammen ref_space->get_element_order(e->id) = H2D_MAKE_QUAD_ORDER(h-ord, v-ord)
//H2D_GET_H_ORDER(id),  H2D_GET_V_ORDER(id)
//
//Problem: Momentan werden Elemente durchlaufen, d.h. Knotenpunkte, die zu mehreren Elementen gehoeren werden mehrmals durchlaufen!
//und in dof_list gespeichert
Element* e;  
AsmList* al = new AsmList();
AsmList* dof_list = new AsmList();
bool more = false;

for_all_active_elements(e, space.get_mesh()){
	 //if( space.get_element_order(e->id)== H2D_MAKE_QUAD_ORDER(1, 1)){	//Ordnung soll 1 
	space.get_element_assembly_list(e, al);
	  	for (unsigned int iv = 0; iv < e->nvert; iv++){   		
		  int index =  space.get_shapeset()->get_vertex_index(iv);
			Node* vn = e->vn[iv];
		  if (space.get_element_order(e->id) == 0) break;
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
//dof_list = NULL;

 
	      // Initialize discrete problem on reference mesh.
		DiscreteProblem* dp_mass = new DiscreteProblem(&massmatrix, ref_space);
		DiscreteProblem* dp_convection = new DiscreteProblem(&convection, ref_space);  	
	//----------------------MassLumping M_L--------------------------------------------------------------------			  
		  UMFPackMatrix* mass_matrix = new UMFPackMatrix();   //M_c/tau
		dp_mass->assemble(mass_matrix, NULL); //rechte Seite braucht nicht berechnet zu werden!
		UMFPackMatrix* lumped_matrix = massLumping(dof_list,mass_matrix);
	//------------------------artificial DIFFUSION D---------------------------------------			
		  UMFPackMatrix* conv_matrix = new UMFPackMatrix();   //K
		dp_convection->assemble(conv_matrix, NULL,true);  //true => erzwingt Diagonleinträge ggf. = 0!		
		 UMFPackMatrix* diffusion = artificialDiffusion(dof_list,conv_matrix);		
	//--------------------------------------------------------------------------------------------
		 UMFPackMatrix* low_matrix = new UMFPackMatrix();  
 		UMFPackMatrix* lowmat_rhs = new UMFPackMatrix(); 
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
		lumped_matrix->multiply_with_scalar(time_step);  // M_L

//-------------rhs lower Order M_L/tau+ (1-theta)(K+D) u^n------------		 
		scalar* coeff_vec = new scalar[ref_dof];
		scalar* coeff_vec_2 = new scalar[ref_dof];				
		UMFPackVector* rhs = new UMFPackVector(ref_dof);
		// Project the previous solution on the FE space
		for(int i=0; i<ref_dof;i++) coeff_vec[i]=0.0; 
		//OGProjection::project_global(ref_space,&sln_prev_time, coeff_vec, matrix_solver, HERMES_L2_NORM);
		if(ts==1)Lumped_Projection::project_lumped(ref_space, &sln_prev_time, coeff_vec, matrix_solver,lumped_matrix);
		else{ 	
			Lumped_Projection::project_lumped(ref_space, &sln_prev_time, coeff_vec, matrix_solver, lumped_matrix);
			 OGProjection::project_global(ref_space,&sln_prev_time, coeff_vec_2, matrix_solver, HERMES_L2_NORM);
			lumped_flux_limiter(dof_list,mass_matrix, lumped_matrix, coeff_vec, coeff_vec_2);
		}
		lowmat_rhs->multiply_with_vector(coeff_vec, coeff_vec_2);  
		rhs->zero(); rhs->add_vector(coeff_vec_2); 
		
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
		
	antidiffusiveFlux(dof_list,mass_matrix,lumped_matrix,conv_matrix,diffusion,rhs, u_L,coeff_vec); //Fluesse in coeff_vec	

//-----------------------korrigierte Loesung-----------------------------	
                   
	for(int i= 0; i<ref_dof; i++) coeff_vec_2[i]= u_L[i]+ (coeff_vec[i]*time_step/lumped_matrix->get(i,i));
	Solution::vector_to_solution(coeff_vec_2, ref_space, &ref_sln);



	      // Project the fine mesh solution onto the coarse mesh.
	      //info("Projecting fine mesh solution on coarse mesh for error estimation.");
	     OGProjection::project_global(&space, &ref_sln, &sln, matrix_solver, HERMES_L2_NORM); 

	      // Calculate element errors and total error estimate.
	      //info("Calculating error estimate.");
	      Adapt* adaptivity = new Adapt(&space, HERMES_L2_NORM);
	    double err_est_rel_total = adaptivity->calc_err_est(&sln, &ref_sln) * 100;

	      // Report results.
	      //info("ndof: %d, ref_ndof: %d, err_est_rel: %g%%", 
		   //Space::get_num_dofs(&space), Space::get_num_dofs(ref_space), err_est_rel_total);
		info("err_est_rel: %g%%",err_est_rel_total);

	      // If err_est too large, adapt the mesh.
	      if (err_est_rel_total < ERR_STOP){
			info("Adaptivity-Stop: err_est_rel_total");
		 	done = true;
			
	      }else if(as>=ADAPSTEP_MAX){
			info("Adaptivity-Stop: ADAPSTEP_MAX");
		 	done = true;		
	     }else{ 	      
		info("Adapting the coarse mesh.");
		done = adaptivity->adapt(&selector, THRESHOLD, STRATEGY, MESH_REGULARITY);	
		info("Adapted.");	
		if (Space::get_num_dofs(&space) >= NDOF_STOP){ 
			info("Adaptivity-Stop:NDOF_STOP,ndof: %d,refdof: %d ",Space::get_num_dofs(&space),Space::get_num_dofs(ref_space));
		  done = true;
		}else
		  // Increase the counter of performed adaptivity steps.
		  as++;
	      }
		
	   
	      // Visualize the solution and mesh.
	      char title[100];
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
	      ord2view.show(ref_space);

	      // Clean up.
		//matrices
		delete mass_matrix;		
		delete lumped_matrix; 
		delete conv_matrix;		
		delete diffusion;
		delete low_matrix;
		delete lowmat_rhs;
		//solver
		delete lowOrd;
		
		//rest
	     	delete adaptivity;
	      	delete ref_space;
		//delete ref_mesh;
	      	delete dp_mass;
	      	delete dp_convection;
		delete rhs;
		delete [] coeff_vec;
		delete [] coeff_vec_2;
		delete al;
		delete dof_list;				
		
	  
    }
    while (done == false);


    // Copy last reference solution into sln_prev_time.
    sln_prev_time.copy(&ref_sln);


  // Update global time.
  current_time += time_step;


  // Increase time step counter
  ts++;
}
while (current_time < T_FINAL); 
  

  // Wait for the view to be closed.
  View::wait();
  return 0;
}

