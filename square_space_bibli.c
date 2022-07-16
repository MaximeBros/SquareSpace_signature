#include <time.h>
#include <flint/fq_poly.h>
#include <stdio.h>
#include <bsd/stdlib.h>

/*
#define Q 11
#define M 5
#define R 2
#define T (R * (R+1)) / 2
#define D1 0
#define D2 2
#define D3 5
#define C1 9
#define C2 10
#define C3 1
*/


#define Q 1451
#define M 131
#define R 4
#define T (R * (R+1)) / 2
#define D1 0
#define D2 1
#define D3 131
#define C1 60
#define C2 2
#define C3 1
	
void random_fqm_list(fq_t* tab_poly, fmpz_t p, fq_ctx_t ctx_fqm, fmpz_mod_ctx_t ctx_fq){
	slong *random_p_tab = malloc(sizeof(slong) * R * M);

	for(int i = 0; i < R * M; i++){
		random_p_tab[i] = arc4random_uniform(Q);
	}
		
    fmpz_mod_poly_t poly; 
    fq_t elt_fqm; 
    fq_init(elt_fqm, ctx_fqm); 
    fmpz_mod_poly_init(poly, ctx_fq); 
     
	for(int i = 0; i < R; i++){
		fmpz_mod_poly_zero(poly, ctx_fq);
        for(int j = 0; j < M; j++){
    		fmpz_mod_poly_set_coeff_ui(poly, j, random_p_tab[i * M + j], ctx_fq);
    	}
    	fq_set_fmpz_mod_poly(elt_fqm, poly, ctx_fqm);
		fq_set(tab_poly[i], elt_fqm, ctx_fqm);
    }
    fmpz_mod_poly_clear(poly, ctx_fq);
    fq_clear(elt_fqm, ctx_fqm);
	free(random_p_tab);
}

void square(fq_t* tab_square, const fq_t* tab_elt, slong nb_elts, fq_ctx_t ctx_fqm){
	int k = 0;
	for(int i = 0; i < nb_elts; i++){
		for(int j = i; j < nb_elts; j++){
    		fq_mul(tab_square[k], tab_elt[i], tab_elt[j], ctx_fqm);
    		k++;
		}
    }
}

void multiply(fq_t* tab_res, const fq_t* tab_elt_1, const fq_t* tab_elt_2, slong nb_elts, fq_ctx_t ctx_fqm){
	int k = 0;
	for(int i = 0; i < nb_elts; i++){
		for(int j = 0; j < nb_elts; j++){
    		fq_mul(tab_res[k], tab_elt_1[i], tab_elt_2[j], ctx_fqm);
    		k++;
		}
    }
}

long fqm_list_to_mat(nmod_mat_t mat, fmpz_t p, const fq_t* fqm_list, slong length, fq_ctx_t ctx_fqm){
    fmpz_mod_mat_t mat_buffer;
    fmpz_mod_mat_init(mat_buffer, M, 1, p);

    for(int j = 0; j < length; j++){
        fq_get_fmpz_mod_mat(mat_buffer, fqm_list[j], ctx_fqm);
		for(int i = 0; i < M; i++){
		    nmod_mat_entry(mat,i,j) = *fmpz_mod_mat_entry(mat_buffer,i,0);
		}
    }
    fmpz_mod_mat_clear(mat_buffer);
   
    return nmod_mat_rref(mat);
}

int main(void)
{
	srand(time(NULL));
    clock_t t;
    t = clock();
    
    fmpz_t p;
    fmpz_init_set_ui(p, Q);
    
    // creation du module et du ctx fqm
    fmpz_mod_ctx_t ctx_fq;
    fmpz_mod_ctx_init(ctx_fq, p);

    fmpz_mod_poly_t modulus;   
    fmpz_mod_poly_init(modulus, ctx_fq);
    
    fmpz_mod_poly_set_coeff_si(modulus, D1, C1, ctx_fq);
    fmpz_mod_poly_set_coeff_si(modulus, D2, C2, ctx_fq);
    fmpz_mod_poly_set_coeff_si(modulus, D3, C3, ctx_fq);
    fq_ctx_t ctx_fqm;
    fq_ctx_init_modulus(ctx_fqm, modulus, ctx_fq, "x");
    fq_ctx_print(ctx_fqm);

    //////////    E generation
    fq_t  *tab_E = malloc(sizeof(fq_t) * R);
    for(int i = 0; i < R; i++){
        fq_init(tab_E[i], ctx_fqm);
    } 
    random_fqm_list(tab_E, p, ctx_fqm, ctx_fq);

    //////////    U generation
    fq_t *tab_U = malloc(sizeof(fq_t) * T);
    for(int i = 0; i < T; i++){
        fq_init(tab_U[i], ctx_fqm);
    }
    square(tab_U, tab_E, R, ctx_fqm);
	
	/*
    for(int i = 0; i < T; i++){
    	fq_print_pretty(tab_U[i], ctx_fqm);printf("\n");
    }
    for(int i = 0; i < R; i++){
    	fq_print_pretty(tab_E[i], ctx_fqm);printf("\n");
    }*/

	//////////    E in matrix form
    ulong q = Q;
    nmod_mat_t mat_E;
    nmod_mat_init(mat_E, M, R, q);
    long rank = fqm_list_to_mat(mat_E, p, tab_E, R, ctx_fqm);
    printf("rank = %ld\n", rank);
    
    //////////    U in matrix form
    nmod_mat_t mat_U;
    nmod_mat_init(mat_U, M, T, q);
    rank = fqm_list_to_mat(mat_U, p, tab_U, T, ctx_fqm);
    printf("rank = %ld\n", rank); 
    
    t = clock() - t;
    double time_taken = ((double)t)/CLOCKS_PER_SEC; // calculate the elapsed time
    printf("The program took %f seconds to generate keys\n", time_taken);
    
    t = clock() - t;
    uint8_t *hash_message = malloc(sizeof(uint8_t) * 128);
	for(int i = 0; i < 128; i++){
		hash_message[i] = arc4random() % 2;
	}
	
    for(int z = 0; z < 1; z++){
        //////////    X generation
		fq_t *tab_X = malloc(sizeof(fq_t) * R);
		for(int i = 0; i < R; i++){
		    fq_init(tab_X[i], ctx_fqm);
		} 
		random_fqm_list(tab_X, p, ctx_fqm, ctx_fq);
		
		//////////    K generation
		fq_t *tab_K = malloc(sizeof(fq_t) * T);
		for(int i = 0; i < T; i++){
		    fq_init(tab_K[i], ctx_fqm);
		}
		square(tab_K, tab_X, R, ctx_fqm);
		
		//////////    K in matrix form
		nmod_mat_t mat_K;
		nmod_mat_init(mat_K, M, T, q);
		rank = fqm_list_to_mat(mat_K, p, tab_K, T, ctx_fqm);
		//printf("rank = %ld\n", rank);
		
		if(hash_message[z] == 1){
			//////////    XE generation
			fq_t *tab_XE = malloc(sizeof(fq_t) * R * R);
			for(int i = 0; i < R * R; i++){
				fq_init(tab_XE[i], ctx_fqm);
			}
			multiply(tab_XE, tab_X, tab_E, R, ctx_fqm);
			//////////    XE in matrix form
			nmod_mat_t mat_XE;
			nmod_mat_init(mat_XE, M, R * R, q);
			rank = fqm_list_to_mat(mat_XE, p, tab_XE, R * R, ctx_fqm);
			printf("rank = %ld\n", rank); 
			
			//////////    Y^2 generation
			fq_t *tab_Y2 = malloc(sizeof(fq_t) * (int)(((R * R + 1) * R * R) / 2));
			for(int i = 0; i < (int)(((R * R + 1) * R * R) / 2); i++){
				fq_init(tab_Y2[i], ctx_fqm);
			}
			square(tab_Y2, tab_XE, R * R, ctx_fqm);
			
			//////////    Y^2 in matrix form
			nmod_mat_t mat_Y2;
			nmod_mat_init(mat_Y2, M, (int)(((R * R + 1) * R * R) / 2), q);
			rank = fqm_list_to_mat(mat_Y2, p, tab_Y2, (int)(((R * R + 1) * R * R) / 2), ctx_fqm);
			printf("rank Y^2 = %ld\n", rank);
		}
    }

    t = clock() - t;
    time_taken = ((double)t)/CLOCKS_PER_SEC; // calculate the elapsed time
    printf("The program took %f seconds to generate keys\n", time_taken);
    /*
    for(int i = 0; i < R * R; i++){
    	fq_print_pretty(tab_XE[i], ctx_fqm);printf("\n");
    }*/

    
    // ne pas oublier les clears 
    //fq_ctx_clear(ctx);
    fmpz_clear(p);
	
    return EXIT_SUCCESS;
}

