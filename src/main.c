#include <time.h>
#include "flint_utils.h"
#include "square_space.h"
#include "parameters.h"

int main(void)
{
	srand(time(NULL));
    clock_t t = clock();
    
    ////////////////////////////    KEY GEN     ///////////////////////
    
    //////////    Initialization
    fmpz_t p;
    fmpz_init_set_ui(p, P);	    
    fmpz_mod_ctx_t ctx_fp;
    fmpz_mod_poly_t modulus;
    fq_ctx_t ctx_fpm;
    init_modulus(&modulus, &ctx_fp, &ctx_fpm, p);

	fq_t *array_E, *array_U;
	nmod_mat_t mat_E, mat_U;
	
	init_array(&array_E, R, ctx_fpm); 
	init_array(&array_U, T, ctx_fpm);
    
    //////////    E generation
    random_fqm_list(array_E, ctx_fpm, ctx_fp);

    //////////    U generation
    square(array_U, array_E, R, ctx_fpm);
    
	//////////    E in matrix form
    fqm_list_to_mat(mat_E, p, array_E, R, ctx_fpm);
     
    //////////    U in matrix form
    fqm_list_to_mat(mat_U, p, array_U, T, ctx_fpm);
    
    t = clock() - t;
    printf("The program took %f seconds to generate keys\n", ((double)t)/CLOCKS_PER_SEC);

    ////////////////////////////    SIGNATURE     ///////////////////////
    
    t = clock();
    fq_t *array_k, *array_xe;
	nmod_mat_t mat_X, mat_K, mat_XE;
	
	fq_t*      *list_x    = malloc(sizeof(fq_t*) * 128);
	nmod_mat_t *commits   = malloc(sizeof(nmod_mat_t) * 128);
	nmod_mat_t *responses = malloc(sizeof(nmod_mat_t) * 128);
	init_array(&array_k, T, ctx_fpm);
	init_array(&array_xe, R * R, ctx_fpm);
	
    for(int z = 0; z < 128; z++){
        //////////    X generation
        init_array(&list_x[z], R, ctx_fpm);
		random_fqm_list(list_x[z], ctx_fpm, ctx_fp);
		
		//////////    K generation
		square(array_k, list_x[z], R, ctx_fpm);
		
		//////////    K in matrix form
		fqm_list_to_mat(mat_K, p, array_k, T, ctx_fpm);
		*commits[z] = *mat_K;
    }    
    
	//////////    Hash from commits, public key and message
	uint8_t *commits_bytes, *mat_U_bytes, *message_bytes, *hash;
    generate_hash(&hash, &commits_bytes, &mat_U_bytes, &message_bytes, commits, mat_U);
  
  	//////////    Responses generation
	int count_chall_1 = 0;
	for(int i = 0; i < 128; i++){
		//////////    If i-th bit is set to 0, compute X matrix, else compute XE matrix
		if((hash[i / 8] >> i % 8) & 1){
			//////////    X in matrix form
			fqm_list_to_mat(mat_X, p, list_x[i], R, ctx_fpm);
			*responses[i] = *mat_X;
		} else {
			//////////    XE generation
			multiply(array_xe, list_x[i], array_E, R, ctx_fpm);
			
			//////////    XE in matrix form
			fqm_list_to_mat(mat_XE, p, array_xe, R * R, ctx_fpm);
    		*responses[i] = *mat_XE;
    		count_chall_1++;
		}
	}

    t = clock() - t;
    printf("The program took %f seconds to sign\n", ((double)t)/CLOCKS_PER_SEC);
    
    //////////    Writing signature into file
    slong size_responses = (11 * M * R2 * count_chall_1) / 8 + (11 * M * R * (128 - count_chall_1)) / 8;
    if(count_chall_1 % 2 == 1) size_responses++; // Last byte truncated
    
    uint8_t *responses_bytes = malloc(sizeof(uint8_t) * size_responses);
	responses_to_bytes(responses_bytes, responses, hash);
    write_signature(commits_bytes, responses_bytes, size_responses);

	clear_vars(&array_E, &array_U, list_x, &array_k, &array_xe, &modulus, &ctx_fp, &ctx_fpm, &p, &hash, &commits_bytes, &mat_U_bytes, &message_bytes);
    return EXIT_SUCCESS;
}
