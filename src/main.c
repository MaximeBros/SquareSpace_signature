#include <time.h>
#include "flint_utils.h"
#include "square_space.h"
#include "parameters.h"

int main(void)
{
	srand(time(NULL));
    clock_t t = clock();
    
    //////////    Initialization
    fmpz_t p;
    fmpz_init_set_ui(p, P);	    
    fmpz_mod_ctx_t ctx_fp;
    fmpz_mod_poly_t modulus;
    fq_ctx_t ctx_fpm;
    init_modulus(&modulus, &ctx_fp, &ctx_fpm, p);

	fq_t *array_E, *array_U, *array_X, *array_K, *array_XE, *array_Y2;
	nmod_mat_t mat_E, mat_U, mat_K, mat_XE, mat_Y2;
	init_array(&array_E, R, ctx_fpm); 
	init_array(&array_U, T, ctx_fpm);
	init_array(&array_X, R, ctx_fpm);
	init_array(&array_K, T, ctx_fpm);
	init_array(&array_XE, R * R, ctx_fpm);
	init_array(&array_Y2, T2, ctx_fpm);
    
    //////////    E generation
    random_fqm_list(array_E, ctx_fpm, ctx_fp);

    //////////    U generation
    square(array_U, array_E, R, ctx_fpm);
    
	//////////    E in matrix form
    long rank = fqm_list_to_mat(mat_E, p, array_E, R, ctx_fpm);
    printf("rank = %ld\n", rank);
     
    //////////    U in matrix form
    rank = fqm_list_to_mat(mat_U, p, array_U, T, ctx_fpm);
    printf("rank = %ld\n", rank);
    
    t = clock() - t;
    printf("The program took %f seconds to generate keys\n", ((double)t)/CLOCKS_PER_SEC);
    t = clock();
    for(int z = 0; z < 128; z++){
        //////////    X generation
		random_fqm_list(array_X, ctx_fpm, ctx_fp);
		
		//////////    K generation
		square(array_K, array_X, R, ctx_fpm);
		
		//////////    K in matrix form
		rank = fqm_list_to_mat(mat_K, p, array_K, T, ctx_fpm);
		

		if(arc4random() % 2){	
			//////////    XE generation
			multiply(array_XE, array_X, array_E, R, ctx_fpm);
			
			//////////    XE in matrix form
			rank = fqm_list_to_mat(mat_XE, p, array_XE, R * R, ctx_fpm);
    		
			//////////    Y^2 generation
			square(array_Y2, array_XE, R * R, ctx_fpm);

			//////////    Y^2 in matrix form
			rank = fqm_list_to_mat(mat_Y2, p, array_Y2, T2, ctx_fpm);
		}
    }
    t = clock() - t;
    printf("The program took %f seconds to sign\n", ((double)t)/CLOCKS_PER_SEC);
    
	clear_vars(&array_E, &array_U, &array_X, &array_K, &array_XE, &array_Y2, &modulus, &ctx_fp, &ctx_fpm, &p);
    return EXIT_SUCCESS;
}

/*
printf("E = \n");
for(int i = 0; i < R; i++){
	fq_print_pretty(array_E[i], ctx_fqm);printf("\n");
}
*/
