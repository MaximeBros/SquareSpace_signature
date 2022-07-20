#include "flint_utils.h"

/**
 * Initialize context over Fp and (Fp)^m, with its modulus
 * @param fmpz_mod_poly_t*  modulus :  Modulus over Fp^m
 * @param fmpz_mod_ctx_t*   ctx_fp  :  Context over Fp
 * @param fq_ctx_t*         ctx_fpm :  Context over Fp^m (modulus)
 * @param fmpz_t            p       :  Characteristic
 */
void init_modulus(fmpz_mod_poly_t* modulus, fmpz_mod_ctx_t* ctx_fp, fq_ctx_t* ctx_fpm, fmpz_t p){
	fmpz_mod_ctx_init(*ctx_fp, p);
    fmpz_mod_poly_init(*modulus, *ctx_fp);
    
    fmpz_mod_poly_set_coeff_si(*modulus, D1, C1, *ctx_fp);
    fmpz_mod_poly_set_coeff_si(*modulus, D2, C2, *ctx_fp);
    fmpz_mod_poly_set_coeff_si(*modulus, D3, C3, *ctx_fp);

    fq_ctx_init_modulus(*ctx_fpm, *modulus, *ctx_fp, "x");
    fq_ctx_print(*ctx_fpm);
}

/**
 * Allocating memory for an array of elements, initializing each of them
 * @param fq_t**    array   :  List initialized
 * @param slong     size    :  List size
 * @param fq_ctx_t  ctx_fpm :  Context over (Fp)^m (modulus)
 */
void init_array(fq_t* *array, slong size, fq_ctx_t ctx_fpm){
    *array = malloc(sizeof(fq_t) * size);
	for(int i = 0; i < size; i++)
		fq_init((*array)[i], ctx_fpm);
}


/**
 * Free memory for an array of elements, clearing each of them
 * @param fq_t**    array   :  List initialized
 * @param slong     size    :  List size
 * @param fq_ctx_t  ctx_fpm :  Context over (Fp)^m (modulus)
 */
void clear_array(fq_t* *array, slong size, fq_ctx_t ctx_fpm){
	for(int i = 0; i < size; i++)
		fq_clear((*array)[i], ctx_fpm);
	free(*array);
}

/**
 * Free memory for every tabs, context and modulus
 */
void clear_vars(fq_t* *array_E, fq_t* *array_U, fq_t* *array_X, fq_t* *array_K, fq_t* *array_XE, fq_t* *array_Y2, fmpz_mod_poly_t* modulus, fmpz_mod_ctx_t* ctx_fp, fq_ctx_t* ctx_fpm, fmpz_t* p){
	clear_array(array_E, R, *ctx_fpm);
	clear_array(array_U, T, *ctx_fpm);
	clear_array(array_X, R, *ctx_fpm);
	clear_array(array_K, T, *ctx_fpm);
	clear_array(array_XE, R * R, *ctx_fpm);
	clear_array(array_Y2, T2, *ctx_fpm);
	
    fmpz_mod_poly_clear(*modulus, *ctx_fp);
    fq_ctx_clear(*ctx_fpm);
    fmpz_mod_ctx_clear(*ctx_fp);
    fmpz_clear(*p);
}

/**
 * Converting a list of elements over (Fp)^m to its matrix representation
 * @param nmod_mat_t  mat       :  Matrix returned
 * @param fmpz_t      p         :  Characteristic
 * @param fq_t*       fqm_list  :  List of elements over (Fp)^m
 * @param slong       length    :  List size
 * @param fq_ctx_t    ctx_fpm   :  Context over (Fp)^m (modulus)
 */
long fqm_list_to_mat(nmod_mat_t mat, fmpz_t p, const fq_t* fqm_list, slong length, fq_ctx_t ctx_fpm){
	nmod_mat_init(mat, M, length, P);
    fmpz_mod_mat_t mat_buffer;
    fmpz_mod_mat_init(mat_buffer, M, 1, p);
    
    for(int j = 0; j < length; j++){
        fq_get_fmpz_mod_mat(mat_buffer, fqm_list[j], ctx_fpm);
		for(int i = 0; i < M; i++){
		    nmod_mat_entry(mat,i,j) = *fmpz_mod_mat_entry(mat_buffer,i,0);
		}
    }
    fmpz_mod_mat_clear(mat_buffer);
    
    long rank = nmod_mat_rref(mat);
    return rank;
}

/**
 * Converting matrix of elements over (Fp)^m to its list representation
 * @param nmod_mat_t       mat      :  Matrix of elements over (Fp)^m
 * @param fq_t*            fqm_list :  List returned
 * @param slong            length   :  List size
 * @param fq_ctx_t         ctx_fpm  :  Context over (Fp)^m (modulus)
 * @param fmpz_mod_ctx_t*  ctx_fp   :  Context over Fp
 */
void mat_to_fqm_list(const nmod_mat_t mat, fq_t* fqm_list, slong length, fq_ctx_t ctx_fpm, fmpz_mod_ctx_t ctx_fp){
    fmpz_mod_poly_t poly; 
	fq_t elt_fqm; 
    fq_init(elt_fqm, ctx_fpm); 
    fmpz_mod_poly_init(poly, ctx_fp); 
     
	for(int i = 0; i < length; i++){
		fmpz_mod_poly_zero(poly, ctx_fp);
        for(int j = 0; j < M; j++){
    		fmpz_mod_poly_set_coeff_ui(poly, j, nmod_mat_get_entry(mat, j, i), ctx_fp);
    	}
    	fq_set_fmpz_mod_poly(elt_fqm, poly, ctx_fpm);
		fq_set(fqm_list[i], elt_fqm, ctx_fpm);
    }
    fmpz_mod_poly_clear(poly, ctx_fp);
    fq_clear(elt_fqm, ctx_fpm);
}

/**
 * Write matrix of elements over (Fp)^m into file
 * @param nmod_mat_t  mat     :  Matrix of elements over (Fp)^m
 * @param slong       size    :  Number of elements
 * @param fq_ctx_t    ctx_fpm :  Context over (Fp)^m (modulus)
 * @param FILE*       file    :  Output file
 */
void write_file_mat(nmod_mat_t* mat, slong size, fq_ctx_t ctx_fpm, FILE* file){
	 uint64_t bitBuffer = 0;  
     uint8_t count = 0;
     for(int i = 0; i < size; i++)
    	for(int j = 0; j < M; j++){
			bitBuffer |= nmod_mat_entry(*mat,i,j);
			count = count + 3;
			if(count >= 8){
				fwrite(&bitBuffer, 1, 1, file);
				bitBuffer >>= 8;
				count = count - 8;
			}
			fwrite(&bitBuffer, 1, 1, file);
			bitBuffer >>= 8;
    	}
}

/**
 * Write list of elements over (Fp)^m into file
 * @param fq_t*     fqm_list :  List of elements over (Fp)^m
 * @param slong     size     :  Number of elements
 * @param fq_ctx_t  ctx_fpm  :  Context over (Fp)^m (modulus)
 * @param FILE*     file     :  Output file
 */
void write_file_array(fq_t* *fqm_list, slong size, fq_ctx_t ctx_fpm, FILE* file){
     uint64_t bitBuffer = 0;  
     uint8_t count = 0;
     for(int i = 0; i < size; i++)
    	for(int j = 0; j < M; j++){
			bitBuffer |= (*fqm_list[i])->coeffs[j];
			count = count + 3;

			if(count >= 8){
				fwrite(&bitBuffer, 1, 1, file);
				bitBuffer >>= 8;
				count = count - 8;
			}
			fwrite(&bitBuffer, 1, 1, file);
			bitBuffer >>= 8;
    	}
}
