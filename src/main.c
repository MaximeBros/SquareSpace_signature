#include <time.h>
#include "flint_utils.h"
#include "square_space.h"
#include "parameters.h"
#include "sha2.h"
    
fmpz_t p;    
fmpz_mod_ctx_t ctx_fp;
fmpz_mod_poly_t modulus;
fq_ctx_t ctx_fpm;
    
void key_gen(fq_t *array_E, fq_t *array_U, nmod_mat_t mat_E, nmod_mat_t mat_U){
	//////////    E generation
    random_fqm_list(array_E, ctx_fpm, ctx_fp);

    //////////    U generation
    square(array_U, array_E, R, ctx_fpm);
    
	//////////    E in matrix form
    long rank = fqm_list_to_mat(mat_E, p, array_E, R, ctx_fpm);
    printf("rank = %ld\n", rank);
    
    //////////    U in matrix form
    rank = fqm_list_to_mat(mat_U, p, array_U, T, ctx_fpm);

}

void signature(nmod_mat_t *commits, nmod_mat_t *responses, fq_t *array_E, nmod_mat_t mat_E, nmod_mat_t mat_U){
	fq_t *array_k, *array_xe;
    fq_t* *list_x = malloc(sizeof(fq_t*) * 128);
	nmod_mat_t mat_X, mat_K, mat_XE;
	
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

    uint8_t *hash;
	generate_hash(&hash, commits, mat_U);
	
	int count_chall_1 = 0;
	for(int i = 0; i < 128; i++){
		if((hash[i / 8] >> i % 8) & 1){
			//////////    X in matrix form
			long rank = fqm_list_to_mat(mat_X, p, list_x[i], R, ctx_fpm);
				
			*responses[i] = *mat_X;
		} else {
			//////////    XE generation
			multiply(array_xe, list_x[i], array_E, R, ctx_fpm);
			
			//////////    XE in matrix form
			fqm_list_to_mat(mat_XE, p, array_xe, R2, ctx_fpm);
    		*responses[i] = *mat_XE;

    		count_chall_1++;
		}
		clear_array(&list_x[i], R, ctx_fpm);
	}

	clear_array(&array_k, T, ctx_fpm);
	clear_array(&array_xe, R2, ctx_fpm);
	free(hash);
	free(list_x);
}

void verification(nmod_mat_t *commits, nmod_mat_t *responses, fq_t *array_U, nmod_mat_t mat_U){
    uint8_t *hash;
	generate_hash(&hash, commits, mat_U);
    
	nmod_mat_t mat_response, mat_KU;
	
	//nmod_mat_set_entry(responses[64], 1, 2, 1);
	//nmod_mat_set_entry(commits[64], 1, 2, 1);
	for(int i = 0; i < 128; i++){
		fq_t *array_commit, *array_response, *array_square;
		//////////    Get commitment
		init_array(&array_commit, T, ctx_fpm); 
		mat_to_fqm_list(commits[i], array_commit, T, ctx_fpm, ctx_fp);
		
		//////////    If i-th hash bit is set to 0, compute X matrix squared
		//////////    and compare with commitment K
		if((hash[i / 8] >> i % 8) & 1){
			init_array(&array_response, R, ctx_fpm); 
			init_array(&array_square, T, ctx_fpm);
			
			mat_to_fqm_list(responses[i], array_response, R, ctx_fpm, ctx_fp);
			square(array_square, array_response, R, ctx_fpm);
			long rank = fqm_list_to_mat(mat_response, p, array_square, T, ctx_fpm);
			
			if(nmod_mat_equal(commits[i], mat_response)){
				printf("\nChallenge X  %d : Correct.", i + 1);
				
			} else {
				printf("\nChallenge X  %d : FALSE.", i + 1);	
			}
			clear_array(&array_response, R, ctx_fpm); 
			clear_array(&array_square, T, ctx_fpm);
		//////////   Otherwise, compute XE matrix squared, compute KU, and compare them
		} else {
			init_array(&array_response, R2, ctx_fpm); 
			init_array(&array_square, R2_2, ctx_fpm);
			
			mat_to_fqm_list(responses[i], array_response, R2, ctx_fpm, ctx_fp);
			square(array_square, array_response, R2, ctx_fpm);
			long rank = fqm_list_to_mat_T2(mat_response, p, array_square, R2_2, ctx_fpm);

			multiply(array_square, array_commit, array_U, T, ctx_fpm);
			rank = fqm_list_to_mat(mat_KU, p, array_square, T2, ctx_fpm);
								printf("rank = %ld\n", rank);
			if(nmod_mat_equal(mat_KU, mat_response)){
				printf("\nChallenge XE %d : Correct.", i + 1);
			} else {
				printf("\nChallenge XE %d : FALSE.", i + 1);
			}
			clear_array(&array_response, R2, ctx_fpm); 
			clear_array(&array_square, R2_2, ctx_fpm);
		}
		clear_array(&array_commit, T, ctx_fpm);
		nmod_mat_clear(commits[i]);
		nmod_mat_clear(responses[i]); 
	}
	nmod_mat_clear(mat_response);
	nmod_mat_clear(mat_KU);
	free(hash);
}

int main(void)
{
	fmpz_init_set_ui(p, P);	
    init_modulus(&modulus, &ctx_fp, &ctx_fpm, p);
    fq_t *array_E, *array_U;
	nmod_mat_t mat_E, mat_U;
	nmod_mat_t *commits   = malloc(sizeof(nmod_mat_t) * 128);
	nmod_mat_t *responses = malloc(sizeof(nmod_mat_t) * 128);
	
	for(int i = 0; i < 100; i++){
	    ////////////////////////////     KEY GEN       ///////////////////////
    
		clock_t t = clock();
		init_array(&array_E, R, ctx_fpm); 
		init_array(&array_U, T, ctx_fpm);
		key_gen(array_E, array_U, mat_E, mat_U);
		t = clock() - t;
		printf("The program took %f seconds to generate keys\n", ((double)t)/CLOCKS_PER_SEC);

		////////////////////////////    SIGNATURE     ///////////////////////
		
		t = clock();
		signature(commits, responses, array_E, mat_E, mat_U);
		t = clock() - t;
		printf("The program took %f seconds to sign\n", ((double)t)/CLOCKS_PER_SEC);
	   
		//////////////////////////   VERIFICATION     ///////////////////////
		
	  	t = clock();
		verification(commits, responses, array_U, mat_U);
		t = clock() - t;
		printf("\nThe program took %f seconds to verify\n", ((double)t)/CLOCKS_PER_SEC);
	}

    
    //////////////////////    WRITING INTO FILE     /////////////////////
    
    /*
    slong size_responses = (11 * M * R2 * count_chall_1) / 8 + (11 * M * R * (128 - count_chall_1)) / 8 + 1;
    if(count_chall_1 % 2 == 1) size_responses++; // Last byte truncated
    
    uint8_t *responses_bytes = malloc(sizeof(uint8_t) * size_responses);
	responses_to_bytes(responses_bytes, responses, hash);
    write_signature(commits_bytes, responses_bytes, size_responses);
    */
    
	//////////    Free memory for every tabs, context and modulus
	
    clear_array(&array_E, R, ctx_fpm);
	clear_array(&array_U, T, ctx_fpm);
	nmod_mat_clear(mat_E);
	nmod_mat_clear(mat_U);
	free(commits);
	free(responses);
    fmpz_mod_poly_clear(modulus, ctx_fp);
    fq_ctx_clear(ctx_fpm);
    fmpz_mod_ctx_clear(ctx_fp);
    fmpz_clear(p);
    
    return EXIT_SUCCESS;
}

