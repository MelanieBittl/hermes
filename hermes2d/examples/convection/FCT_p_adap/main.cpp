#define HERMES_REPORT_ALL
#include "hermes2d.h"
 #define PI (3.141592653589793) 

// 1. Step: (M_L/tau -theta(K+D)) u^L =   (M_L/tau + (1-theta)(K+D)) u^n
// 2. Step : f_ij = (M_c)_ij (dt_u_L(i)- dt_u_L(j)) + D_ij (u_L(i)- u_L(j)); f_i = sum_(j!=i) alpha_ij f_ij
// 3. Step:  M_L u^(n+1) = M_L u^L + tau * f 

const int INIT_REF_NUM =5;                   // Number of initial refinements.
const int P_INIT = 2;                             // Initial polynomial degree.
const double time_step = 1e-3;                           // Time step.

const double T_FINAL = 2*PI;                       // Time interval length.

const double NEWTON_TOL = 1e-5;                   // Stopping criterion for the Newton's method.
const int NEWTON_MAX_ITER = 20;                  // Maximum allowed number of Newton iterations.

const double P_ADAP_TOL_INIT = 90;
const double P_ADAP_TOL =90;
const int P_ADAP_MAX_ITER = 5;


const double theta = 0.5;    // theta-Schema fuer Zeitdiskretisierung (theta =0 -> explizit, theta=1 -> implizit)

MatrixSolverType matrix_solver = SOLVER_UMFPACK; 

// Boundary markers.
const std::string BDY_IN = "inlet";
const std::string BDY_OUT = "outlet";

// Weak forms & Projection with masslumping
#include "definitions.cpp"
#include "lumped_projection.cpp"
#include "extrapolation.cpp"
#include "least_square.cpp"
#include "p_only_adapt.cpp"
#include "patch_solution.cpp"


//Mass lumping an den Stellen von AsmList, sonst standard Massmatrix
UMFPackMatrix* massLumping(AsmList* al,UMFPackMatrix* mass_matrix)
{  //al=NULL=>lumped=mass
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
UMFPackMatrix* pure_p1_massLumping(AsmList* al,UMFPackMatrix* mass_matrix)
{  //al=NULL=>lumped=mass
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
UMFPackMatrix* artificialDiffusion(AsmList* al,UMFPackMatrix* conv_matrix)
{ //al=NULL => diffusion=0
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
void antidiffusiveFlux(AsmList* al,UMFPackMatrix* mass_matrix,UMFPackMatrix* lumped_matrix,UMFPackMatrix* conv_matrix,UMFPackMatrix* diffusion,UMFPackVector* flux_dt_rhs, scalar* u_L, scalar* flux_scalar )
{ //al==NULL =>flux=0
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




void lumped_flux_limiter(AsmList* al,UMFPackMatrix* mass_matrix,UMFPackMatrix* lumped_matrix, scalar* u_L, scalar* u_H)
{
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
void p1_list(H1Space* space, AsmList* dof_list, AsmList* al )
{
	Element* e =NULL;
	dof_list->clear();
	bool more = false;
	for_all_active_elements(e, space->get_mesh()){
	 	if((space->get_element_order(e->id)== H2D_MAKE_QUAD_ORDER(1, 1))||(space->get_element_order(e->id)==1)){	//Ordnung soll 1 
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
		}

	}

}


bool p_adap(Space* space, Solution* sln, scalar* sln_coeff, double error_limit, PonlyAdapt* adapt)
{

	//if(sln->get_type() == HERMES_EXACT) info("exakte Loesung");
	

	if(sln_coeff!=NULL){	
		bool all_zero = true;
		for(int i= 0; i<space->get_num_dofs(); i++){
			if(sln_coeff[i]!=0.0) all_zero =false;
		} 
		if(all_zero==true) info("p_adap: all zero");
	}


	/*ScalarView exview("extrapolierte Loesung", new WinGeom(0, 0, 500, 400));
	ScalarView lsview("LS- Loesung", new WinGeom(0, 700, 500, 400));	
	OrderView mview("mesh", new WinGeom(0, 0, 350, 350));
	ScalarView initview("Solution", new WinGeom(700, 0, 500, 400));*/

	int counter = 0;
	Element* e =NULL;
	bool coarse = false;  //->vergroebern ?

	double error_ex = 0.0;
	double error_ls = 0.0;
	double err_max =0.0;
	 double l2_norm =1.0;

	bool* elements_to_refine = new bool[space->get_mesh()->get_max_element_id()+1];
	memset(elements_to_refine, false, (space->get_mesh()->get_max_element_id()+1)*sizeof(bool));

	for(int i=0; i<=space->get_mesh()->get_max_element_id(); i++) if(elements_to_refine[i]!=false) info("memset");

	for_all_active_elements(e, space->get_mesh()){	// e ist Ausgangselement
		int v_ord = 0.0; int h_ord =0.0; int index =0.0; 
		int max_ord = 0.0;   //Polynomgrad in Element e, bei Quad: max(v_ord,h_ord)	
		int no_of_tris = 0;
		int no_of_quads = 0;
		int element_id = e->id;
		Node* vn; Node* edge;		

	//----------------------
	//Polynomgrad bestimmen
		/*if(l2_norm->get_element_error_squared(0,e->id)==0.0){
				info("l2_norm ==0->lsg =0");
					elements_to_refine[element_id] = false;
		
				continue; 
		}*/


		if (space->get_element_order(element_id) <=1){
				elements_to_refine[element_id] = false;		
				continue; 
		}   
		if(e->is_triangle()==true){ //Triangle
			max_ord = space->get_element_order(element_id);			
			no_of_tris++;	
		}else if(e->is_quad()==true){ //Quad
			v_ord = H2D_GET_V_ORDER(space->get_element_order(element_id));
			h_ord = H2D_GET_H_ORDER(space->get_element_order(element_id));
			if(v_ord >= h_ord) max_ord = v_ord;
			else max_ord = h_ord;  //max Polynomgrad im Quad		
			no_of_quads++;		
		}else {
			elements_to_refine[element_id] = false;
			continue; 
		}   
		if(max_ord <=1) {
			elements_to_refine[element_id] = false;	
			 continue;
		}   
	

	//bestimme die Anzahl der Elemente die Knoten von e benutzen
		int* number_of_elem= new int[e->nvert];  //Anzahl der Elemente, die Knoten benutzen	
		for (unsigned int iv = 0; iv < e->nvert; iv++){  
			index =  space->get_shapeset()->get_vertex_index(iv);
			vn = e->vn[iv];	
			if(vn->ref>=TOP_LEVEL_REF) number_of_elem[iv] = vn->ref- TOP_LEVEL_REF; 	//keine Ahnung, ob das immer funzt.!!..			
			else 	number_of_elem[iv] = vn->ref;			
		} 

	
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
							if((k<e->nvert)&&(e->en[k]->bnd!=1)&&(e->en[e->prev_vert(k)]->bnd!=1)){
											// Knoten von e und kein Rand!							
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
												if((k<e->nvert)&&(e->en[k]->bnd!=1)&&(e->en[e->prev_vert(k)]->bnd!=1)){
														// Knoten von e und kein Randknoten!!
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
	//mview.show(patch_space);

		ExtrapolationSolution* u_ex = new ExtrapolationSolution(patch_mesh, e_patch);
		LSSolution* u_LS = new LSSolution(patch_mesh,e_patch);		
		//Adapt* err_ex = new Adapt(patch_space, HERMES_L2_NORM);
		//Adapt* err_ls = new Adapt(patch_space, HERMES_L2_NORM);
		
		ErrorEstimation* err_ex = new ErrorEstimation(patch_space, HERMES_L2_NORM,e_patch);
		ErrorEstimation* err_ls = new ErrorEstimation(patch_space, HERMES_L2_NORM, e_patch);



		if(sln->get_type() == HERMES_EXACT) {	
			CustomInitialCondition* u_init= new CustomInitialCondition(patch_mesh);

			  //hier ist es moeglich, dass error_ls>error_ex, da u_init exact, u_ex und u_ls aber ein projeziertes u_init verwenden
		/*	scalar* coeff_vec = new scalar[patch_space->get_num_dofs()];				
			OGProjection::project_global(patch_space,u_init, coeff_vec, matrix_solver, HERMES_L2_NORM); 
			Solution::vector_to_solution(coeff_vec,patch_space, u_init);
*/

			u_ex->set_fn(u_init);
			u_ex->init_extrapolation(patch_space);		
			u_LS->LS_project(u_init, patch_space,matrix_solver);

			error_ex = err_ex->calc_err_est(u_ex, u_init,true, HERMES_TOTAL_ERROR_ABS);
			error_ls = err_ls->calc_err_est(u_LS, u_init,true, HERMES_TOTAL_ERROR_ABS);	
		
			l2_norm = err_ex->calc_l2_norm( u_ex, u_init,patch_space, true);

/*	if(error_ls/error_ex > 1.0)
	{
		for_all_active_elements(ne, patch_space->get_mesh())
			printf("error_ex(%i): %f, error_ls: %f \n",ne->id, err_ex->get_element_error_squared(0,ne->id),err_ls->get_element_error_squared(0,ne->id));
		
		printf("LS-error: %f ; Ex-error: %f \n ",error_ls,error_ex );
		exview.set_min_max_range(-0.1, 1); 
			initview.set_min_max_range(-0.1, 1); 
			lsview.set_min_max_range(-0.1, 1); 
			initview.show(u_init);		
			lsview.show(u_LS);
			exview.show(u_ex);
		View::wait(HERMES_WAIT_KEYPRESS);
	}*/

		//Views:
		/*		exview.set_min_max_range(-0.1, 1); 
			initview.set_min_max_range(-0.1, 1); 
			lsview.set_min_max_range(-0.1, 1); 
			initview.show(u_init);		
			lsview.show(u_LS);
			exview.show(u_ex);*/

			//delete [] coeff_vec;

			if(u_init!=NULL){ delete u_init; u_init = NULL;}

		}else{
				PatchSolution* u_init = new PatchSolution(patch_mesh);
				H1Space* patch_space_init = new H1Space(patch_mesh, P_INIT);


			if(sln_coeff==NULL)	u_init->init(patch_space_init, space, sln, elem_ids, elems);
			else u_init->init(patch_space_init, space, sln, sln_coeff, elem_ids, elems);

			u_ex->set_fn(u_init);
			u_ex->init_extrapolation(patch_space);		
			u_LS->LS_project(u_init, patch_space,matrix_solver);	

			error_ex = err_ex->calc_err_est(u_ex, u_init,true, HERMES_TOTAL_ERROR_ABS);
			error_ls = err_ls->calc_err_est(u_LS, u_init,true, HERMES_TOTAL_ERROR_ABS);	


			//error_ex = err_ex->l2_error(u_ex, u_ex, u_init,u_init,patch_space, true);
			//error_ls = err_ls->l2_error(u_LS, u_LS, u_init,u_init, patch_space, true);

			l2_norm = err_ex->calc_l2_norm( u_ex, u_init,patch_space, true);
	

			if(u_init!=NULL){ delete u_init; u_init = NULL;}
			if(patch_space_init!=NULL){ delete patch_space_init; patch_space_init = NULL;}

			/*if(error_ex > 1e-5) {
			for_all_active_elements(ne, patch_space->get_mesh())
				printf("error_ex(%i): %f\n",ne->id, err_ex->get_element_error_squared(0,ne->id)); //Elementweiser Fehler
			//Views:
					exview.set_min_max_range(-0.1, 1);			
				lsview.set_min_max_range(-0.1, 1); 				
				lsview.show(u_LS);
				exview.show(u_ex);
			initview.set_min_max_range(-0.1, 1);
			initview.show(u_init);	
			View::wait(HERMES_WAIT_KEYPRESS);			


			}*/

			//View::wait(HERMES_WAIT_KEYPRESS);
		}	

		if(error_ls/error_ex > err_max) err_max = error_ls/error_ex;		

		//if((error_ls/(error_ex)< error_limit)&&((error_ex >NEWTON_TOL)&&(error_ls>NEWTON_TOL))){  //error_ex = 0 => u= 0
		if((error_ex/l2_norm > error_limit)&&(l2_norm >0.0)){ 
			coarse = true;
			elements_to_refine[element_id] = true; 
			for(unsigned int j = 0; j <ne_list->cnt; j++)  //Alle Nachbarelemente durchlaufen
				elements_to_refine[ne_list->idx[j]] = true;		
		}
			
		

		//printf("LS-error: %f ; Ex-error: %f \n ",error_ls,error_ex );	
		//if(error_ex!=0.0) { printf("LS-error: %f ; Ex-error: %f  ",error_ls,error_ex );	printf("error(LS/ex)= %f \n",error_ls/error_ex );}
	/*for_all_active_elements(ne, patch_space->get_mesh())
		printf("error_ex(%i): %f, error_ls: %f \n",ne->id, err_ex->get_element_error_squared(0,ne->id),err_ls->get_element_error_squared(0,ne->id));	
		printf("LS-error: %f ; Ex-error: %f  ",error_ls,error_ex );*/
	
		
				

		delete ne_list;
		delete [] number_of_elem;
		delete [] verts;
		delete [] tris;
		delete [] quads;
		delete []  mark;
		delete [] elem_ids;	
		delete patch_mesh;
		delete patch_space;
		delete u_LS;
		delete u_ex;	
		delete err_ex;
		delete err_ls;

	 	counter++;
	}

	info("max(err) = %f", err_max);
	if(coarse==true) coarse = adapt->adapt(elements_to_refine);
	
   

	delete [] elements_to_refine;
	



	return coarse;

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
  Solution low_sln, ref_sln;

  // Previous time level solution (initialized by the initial condition).
  CustomInitialCondition u_prev_time(&mesh);  

  // Initialize the weak formulation.
CustomWeakFormMassmatrix massmatrix(time_step, &u_prev_time);
CustomWeakFormConvection convection(&u_prev_time);


  // Initialize views.
	ScalarView Lowview("niedriger Ordnung", new WinGeom(500, 500, 500, 400));
	//Lowview.show(&u_prev_time, HERMES_EPS_HIGH);
	ScalarView sview("Solution", new WinGeom(0, 500, 500, 400));
	//sview.show(&u_prev_time, HERMES_EPS_HIGH); 
	ScalarView pview("projezierter Anfangswert", new WinGeom(500, 0, 500, 400));
	OrderView mview("mesh", new WinGeom(0, 0, 350, 350));
	mview.show(&space);


//-----------------------------------------------------------------------------------------------------------------------------------------
 	 // Previous time level solution (initialized by the initial condition).
 /*		 CustomInitialCondition u_prev_time_high(&mesh);
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
	*/	

//------------------------------------
	




		  // Initialize the FE problem.
	DiscreteProblem* dp_mass = new DiscreteProblem(&massmatrix, &space);
	DiscreteProblem* dp_convection = new DiscreteProblem(&convection, &space);
	UMFPackMatrix* mass_matrix = new UMFPackMatrix();   //M_c/tau
	UMFPackMatrix* conv_matrix = new UMFPackMatrix();   //K
	UMFPackMatrix* low_matrix = new UMFPackMatrix();  
	UMFPackMatrix* K_D = new UMFPackMatrix(); 
	scalar* u_L = NULL; 

// Time stepping loop:
	double current_time = 0.0; 
	int ts = 1;
	char title[100];
	bool changed = true;
	int ps = 1;

	PonlyAdapt* adapting = new PonlyAdapt(&space, HERMES_L2_NORM);
	if(P_INIT>1){
		while((changed ==true)&&(ps<P_ADAP_MAX_ITER)){
			 changed = p_adap(&space, &u_prev_time, NULL, P_ADAP_TOL_INIT, adapting);
			mview.show(&space);
			ps++;
		}
		if(ndof!=space.get_num_dofs()){ changed =true;ndof = space.get_num_dofs();}
		else changed = false;
	}

	
	//Lowview.show(&u_prev_time, HERMES_EPS_HIGH);

    scalar* coeff_vec = new scalar[ndof];
	for(int i=0; i<ndof;i++) coeff_vec[i]=0.0;	

		//Durchlaufen und p bestimmen: in dof_list dof von vertex-basisfunction speichern 
	AsmList al;	
	AsmList dof_list;
	p1_list(&space, &dof_list, &al );


do
{		 info(" Time step %d, time %3.5f", ts, current_time); 
	ps=1; 	

	do
	{
			scalar* lumped_scalar = new scalar[ndof];
			UMFPackVector* vec_rhs = new UMFPackVector(ndof);
			scalar* coeff_vec_2 = new scalar[ndof];
			for(int i=0; i<ndof;i++) coeff_vec_2[i]=0.0;
			scalar* flux_scalar = new scalar[ndof]; 

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

			K_D->create(conv_matrix->get_size(),conv_matrix->get_nnz(), conv_matrix->get_Ap(), conv_matrix->get_Ai(),conv_matrix->get_Ax());
			K_D->add_matrix(diffusion); 
			low_matrix->create(K_D->get_size(),K_D->get_nnz(), K_D->get_Ap(), K_D->get_Ai(),K_D->get_Ax());
			//(-theta)(K+D)
			if(theta==0) low_matrix->zero();
			else	low_matrix->multiply_with_scalar(-theta);
			//(1-theta)(K+D)
			if(theta ==1) K_D->zero();
			else K_D->multiply_with_scalar((1.0-theta));

			//M_L/tau - theta(D+K)
			low_matrix->add_matrix(lumped_matrix);  //kann nun fuer alle Zeitschritte verwendet werden (ausser bei Adaptivitaet)
			//M_L/tau+(1-theta)(K+D)
			K_D->add_matrix(lumped_matrix);	

			lumped_matrix->multiply_with_scalar(time_step);  // M_L
			mass_matrix->multiply_with_scalar(time_step);  // massmatrix = M_C
		
			// Project the initial condition on the FE space->coeff_vec	
			if((changed==true)||(ts==1)){ 				
				/*UMFPackMatrix* proj_matrix = pure_p1_massLumping(&dof_list,mass_matrix);
				Lumped_Projection::project_lumped_rhs(&space, &u_prev_time, coeff_vec, matrix_solver, proj_matrix);
				OGProjection::project_global(&space,&u_prev_time, coeff_vec_2, matrix_solver, HERMES_L2_NORM);		
				lumped_flux_limiter(&dof_list,mass_matrix, proj_matrix, coeff_vec, coeff_vec_2);
				*/				
				Lumped_Projection::project_lumped(&space, &u_prev_time, coeff_vec, matrix_solver, lumped_matrix);
				OGProjection::project_global(&space,&u_prev_time, coeff_vec_2, matrix_solver, HERMES_L2_NORM);
				lumped_flux_limiter(&dof_list,mass_matrix, lumped_matrix, coeff_vec, coeff_vec_2);	
				//OGProjection::project_global(&space,&u_prev_time, coeff_vec, matrix_solver, HERMES_L2_NORM); 
				Solution::vector_to_solution(coeff_vec, &space, &u_prev_time);
				pview.show(&u_prev_time);
				
			}



		/*
			info("------HERMES: ");
		// Perform Newton's iteration.		
			if (!hermes2d.solve_newton(coeff_vec_newton, &dp, solver, matrix, rhs, jacobian_changed,
						   NEWTON_TOL, NEWTON_MAX_ITER, verbose)) error("Newton's iteration failed.");
			// Update previous time level solution.
			Solution::vector_to_solution(coeff_vec_newton, &space, &u_prev_time_high);
			// Visualize the solution.	
			  sprintf(title, "high_Ord_Hermes Time %3.2f", current_time);
			  Hermview.set_title(title);
			  Hermview.show(&u_prev_time_high);

		*/



	//-------------rhs lower Order M_L/tau+ (1-theta)(K+D) u^n------------		
			K_D->multiply_with_vector(coeff_vec, lumped_scalar); 
			vec_rhs->zero(); vec_rhs->add_vector(lumped_scalar);
	//-------------------------solution of lower order------------	
				//info("solving low order solultion");
			  // Solve the linear system and if successful, obtain the solution. M_L/tau-theta(D+K) u^n+1=  M_L/tau+ (1-theta)(K+D) u^n
			Solver* lowOrd = create_linear_solver(matrix_solver,low_matrix,vec_rhs);	
			if(lowOrd->solve()){ 
				u_L = lowOrd->get_solution();  
				Solution::vector_to_solution(u_L, &space, &low_sln);	
			  }else error ("Matrix solver failed.\n");
		//---------------------------------------antidiffusive fluxes-----------------------------------	
				//info("assemble fluxes");	
			 antidiffusiveFlux(&dof_list,mass_matrix,lumped_matrix,conv_matrix,diffusion,vec_rhs, u_L, flux_scalar);		
	
			for(int i= 0; i<ndof; i++) coeff_vec[i]= u_L[i]+ (flux_scalar[i]*time_step/lumped_matrix->get(i,i));
			 Solution::vector_to_solution(coeff_vec, &space, &ref_sln);

			changed = p_adap(&space, &ref_sln, coeff_vec, P_ADAP_TOL, adapting);

			if(changed==true){ 
					if(coeff_vec!=NULL){ delete [] coeff_vec; 	coeff_vec = NULL;}
					ndof = space.get_num_dofs();
					coeff_vec = new scalar[ndof];
					for(int i=0; i<ndof;i++) coeff_vec[i]=0.0;
					//Durchlaufen und p bestimmen: in dof_list dof von vertex-basisfunction speichern 
					p1_list(&space, &dof_list, &al );
				}


			 // Visualize the solution.
			  sprintf(title, "low_Ord Time %3.2f", current_time);
			  Lowview.set_title(title);
			 Lowview.show(&low_sln);	 
			  sprintf(title, "korrigierte Loesung: Time %3.2f", current_time);
			  sview.set_title(title);
			  sview.show(&ref_sln);
				mview.show(&space);

			  // Update global time.
			  current_time += time_step;

			  // Increase time step counter
			  ts++;
		//View::wait(HERMES_WAIT_KEYPRESS);

			  // Clean up.
			if(lowOrd!=NULL)		delete lowOrd; 
			if(lumped_matrix!=NULL)	delete lumped_matrix; 
			if(diffusion!=NULL)		delete diffusion;
			if(lumped_scalar!=NULL)	delete[] lumped_scalar;  
			if(coeff_vec_2!=NULL)	delete[] coeff_vec_2;
			if(flux_scalar!=NULL)	delete[] flux_scalar; 
			if(vec_rhs!=NULL)		delete vec_rhs;

	
	ps++;

	}
	while((changed ==true)&&(ps<P_ADAP_MAX_ITER));

    // Copy last reference solution into sln_prev_time.
    u_prev_time.copy(&ref_sln);
}
while (current_time < T_FINAL);

	if(coeff_vec!=NULL) delete [] coeff_vec;
	if(mass_matrix!=NULL) 	delete mass_matrix;  
	if(conv_matrix!=NULL)	delete conv_matrix;
	if(low_matrix!=NULL)	delete low_matrix;
	if(K_D!=NULL)			delete K_D;
	if(dp_convection!=NULL)	delete dp_convection;
	if(dp_mass!=NULL)		delete dp_mass; 
	if(adapting !=NULL) delete adapting;




  // Wait for the view to be closed.
  View::wait();
  return 0;
}

