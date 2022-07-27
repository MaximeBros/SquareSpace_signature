#include "square_space.h"
#include "arc4random.h"
#include "sha2.h"
#include "flint_utils.h"

/**
 * Generate the private key, which consists of a random list of M elements
 * (polynomials) in Fp^m
 * @param fq_t*           array_poly :  List of polynomials returned
 * @param fmpz_mod_ctx_t  ctx_fp     :  Context over Fp
 * @param fq_ctx_t        ctx_fpm    :  Context over Fp^m (modulus)
 */
void random_fqm_list(fq_t* array_poly, fq_ctx_t ctx_fpm, fmpz_mod_ctx_t ctx_fp){
	// Using the arc4random function from BSD library to generate
	// uniformly R * M random numbers between 0 and M - 1
	uint16_t *random_p_array = malloc(sizeof(uint16_t) * R * M);
	for(int i = 0; i < R * M; i++){
		random_p_array[i] = arc4random_uniform(P);
	}
	
	// Initialization
    fmpz_mod_poly_t poly; 
    fq_t elt_fqm; 
    fq_init(elt_fqm, ctx_fpm); 
    fmpz_mod_poly_init(poly, ctx_fp); 
     
    // Fill all coefficients for all polynomials
	for(int i = 0; i < R; i++){
		fmpz_mod_poly_zero(poly, ctx_fp);
		
        for(int j = 0; j < M; j++)
    		fmpz_mod_poly_set_coeff_ui(poly, j, random_p_array[i * M + j], ctx_fp);
    	
    	fq_set_fmpz_mod_poly(elt_fqm, poly, ctx_fpm);
		fq_set(array_poly[i], elt_fqm, ctx_fpm);
    }
    
    // Free variables
    fmpz_mod_poly_clear(poly, ctx_fp);
    fq_clear(elt_fqm, ctx_fpm);
	free(random_p_array);
}


/**
 * Generate the public key, which consists to make the vector space of
 * all unique couples of elements multiplicated themselves from the 
 * private key
 * @param fq_t*     array_square :  Final list returned
 * @param fq_t*     array_elt    :  List of polynomials to square
 * @param slong     nb_elts      :  # elements in list
 * @param fq_ctx_t  ctx_fpm      :  Context over Fp^m (modulus)
 */
void square(fq_t* array_square, const fq_t* array_elt, slong nb_elts, fq_ctx_t ctx_fpm){
	int k = 0;
	for(int i = 0; i < nb_elts; i++){
		for(int j = i; j < nb_elts; j++){
    		fq_mul(array_square[k], array_elt[i], array_elt[j], ctx_fpm);
    		k++;
		}
    }
}


/**
 * Generate the list of polynomials resulted from the multiplication
 * of all differents couples of elements from two vector spaces
 * (with same number of elements)
 * @param fq_t*     array_res    :  Final list returned
 * @param fq_t*     array_elt_1  :  List of polynomials to square
 * @param fq_t*     array_elt_2  :  List of polynomials to square
 * @param slong     nb_elts      :  # elements in list
 * @param fq_ctx_t  ctx_fpm      :  Context over Fp^m (modulus)
 */
void multiply(fq_t* array_res, const fq_t* array_elt_1, const fq_t* array_elt_2, slong nb_elts, fq_ctx_t ctx_fpm){
	int k = 0;
	for(int i = 0; i < nb_elts; i++){
		for(int j = 0; j < nb_elts; j++){
    		fq_mul(array_res[k], array_elt_1[i], array_elt_2[j], ctx_fpm);
    		k++;
		}
    }
}


/**
 * Generate the SHA256 from commits, public key and message
 * @param uint8_t*     hash           :  Hash returned
 * @param nmod_mat_t*  commits        :  List of matrices K commited
 * @param uint8_t*     commits_bytes  :  List of matrices K commited in 8-byte list
 * @param nmod_mat_t   mat_U          :  Public key U in matrix
 * @param uint8_t*     mat_U_bytes    :  Public key U in 8-byte list
 * @param uint8_t*     message_bytes  :  Message in 8-byte list
 */
void generate_hash(uint8_t* *hash, uint8_t* *commits_bytes, uint8_t* *mat_U_bytes, uint8_t* *message_bytes, const nmod_mat_t *commits, const nmod_mat_t mat_U){
    *commits_bytes = malloc(sizeof(uint8_t) * COMMIT_SIZE);
    *mat_U_bytes   = malloc(sizeof(uint8_t) * PUB_KEY_SIZE);
    *message_bytes = malloc(sizeof(uint8_t) * MESSAGE_SIZE); 
    *hash          = malloc(sizeof(uint8_t) * 32);
    uint8_t* concat = malloc(sizeof(uint8_t) * HASH_SIZE);

	for(int i = 0; i < 4; i++){
		uint32_t number = arc4random();
		for(int j = 0; j < 4; j++)
			(*message_bytes)[i * 4 + j] = (number >> j * 8) & ((1 << 8) - 1);
	}

	list_mat_to_bytes(*commits_bytes, commits, T);
    mat_to_bytes(*mat_U_bytes, mat_U, T);

    for(int i = 0; i < COMMIT_SIZE; i++)
		concat[i] = (*commits_bytes)[i];
	for(int i = COMMIT_SIZE; i < PUB_KEY_SIZE; i++)
		concat[i] = (*mat_U_bytes)[i];
	for(int i = PUB_KEY_SIZE; i < MESSAGE_SIZE; i++)
		concat[i] = (*message_bytes)[i];
	
    sha256ctx state;
    sha256_inc_init(&state);
    sha256(*hash, concat, HASH_SIZE);
    free(concat);
}
