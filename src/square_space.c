#include "square_space.h"
#include "arc4random.h"

u_int32_t 
arc4random(void)
{
        u_int32_t val;
        _ARC4_LOCK();
        arc4_count -= 4;
        if (arc4_count <= 0 || !rs_initialized || arc4_stir_pid != getpid())
                arc4_stir();
        val = arc4_getword();
        _ARC4_UNLOCK();
        return val;
}

u_int32_t
arc4random_uniform(u_int32_t upper_bound)
{
        u_int32_t r, min;
        if (upper_bound < 2)
                return 0;
#if (ULONG_MAX > 0xffffffffUL)
        min = 0x100000000UL % upper_bound;
#else
        /* Calculate (2**32 % upper_bound) avoiding 64-bit math */
        if (upper_bound > 0x80000000)
                min = 1 + ~upper_bound;                /* 2**32 - upper_bound */
        else {
                /* (2**32 - (x * 2)) % x == 2**32 % x when x <= 2**31 */
                min = ((0xffffffff - (upper_bound * 2)) + 1) % upper_bound;
        }
#endif
        /*
         * This could theoretically loop forever but each retry has
         * p > 0.5 (worst case, usually far better) of selecting a
         * number inside the range we need, so it should rarely need
         * to re-roll.
         */
        for (;;) {
                r = arc4random();
                if (r >= min)
                        break;
        }
        return r % upper_bound;
}
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

