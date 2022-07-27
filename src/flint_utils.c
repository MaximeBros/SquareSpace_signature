#include "flint_utils.h"

uint8_t mask = ((1 << 8) - 1);

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
void clear_vars(fq_t* *array_E, fq_t* *array_U, fq_t* *array_X, fq_t* *array_K, fq_t* *array_XE, fmpz_mod_poly_t* modulus, fmpz_mod_ctx_t* ctx_fp, fq_ctx_t* ctx_fpm, fmpz_t* p, uint8_t* *hash, uint8_t* *commits_bytes, uint8_t* *mat_U_bytes, uint8_t* *message_bytes){
	clear_array(array_E, R, *ctx_fpm);
	clear_array(array_U, T, *ctx_fpm);
    for(int z = 0; z < 128; z++){
    	clear_array(&array_X[z], R, *ctx_fpm);
    }
	clear_array(array_K, T, *ctx_fpm);
	clear_array(array_XE, R * R, *ctx_fpm);
	
    fmpz_mod_poly_clear(*modulus, *ctx_fp);
    fq_ctx_clear(*ctx_fpm);
    fmpz_mod_ctx_clear(*ctx_fp);
    fmpz_clear(*p);
    
    free(*hash);
    free(*commits_bytes);
    free(*mat_U_bytes);
    free(*message_bytes);
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
	nmod_mat_init(mat, length, M, P);
	
    for(int i = 0; i < length; i++){
		for(int j = 0; j < M; j++){
		    nmod_mat_entry(mat,i,j) = fqm_list[i]->coeffs[j];
		}
    }
	
    return nmod_mat_rref(mat);
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
 * Converting matrix over (Fp)^m into 8-byte list representation
 * @param uint8_t*    tab_bytes  :  Matrix in 8-byte list returned
 * @param nmod_mat_t  mat        :  Matrix of elements over (Fp)^m
 * @param slong       size       :  # of columns
 */
void mat_to_bytes(uint8_t *tab_bytes, const nmod_mat_t mat, slong size){
	 uint64_t bitBuffer = 0;  
     uint8_t count = 0;
     int k = 0;
     	 
     for(int i = 0; i < size; i++)
    	for(int j = 0; j < M; j++){
			bitBuffer |= nmod_mat_entry(mat,i,j);
			count = count + 3;
			if(count >= 8){
				tab_bytes[k] = (bitBuffer >> count) & mask;
				bitBuffer >>= 8;
				count = count - 8;
				k++;
			}
			tab_bytes[k] = (bitBuffer >> count) & mask;
			bitBuffer >>= 8;
			k++;
    	}
    if(count > 0){
    	tab_bytes[k] = (bitBuffer >> count) & mask;
		bitBuffer >>= 8;
    }
    	
}


/**
 * Converting list of matrices over (Fp)^m into 8-byte list representation
 * @param uint8_t*    tab_bytes   :  List of matrices in 8-byte list returned
 * @param nmod_mat_t*  mat        :  List of matrices of elements over (Fp)^m
 * @param slong       size        :  # of columns
 */
void list_mat_to_bytes(uint8_t *tab_bytes, const nmod_mat_t* mat, slong size){
	 uint64_t bitBuffer = 0;  
     uint8_t count = 0;
	 int k = 0;
	 
     for(int h = 0; h < 128; h++)
		 for(int i = 0; i < size; i++)
			for(int j = 0; j < M; j++){
				bitBuffer |= nmod_mat_entry(mat[h],i,j);
				count = count + 3;
				if(count >= 8){
					tab_bytes[k] = (bitBuffer >> count) & mask;
					bitBuffer >>= 8;
					count = count - 8;
					k++;
				}
				tab_bytes[k] = (bitBuffer >> count) & mask;
				bitBuffer >>= 8;
				k++;
			}
	if(count > 0){
    	tab_bytes[k] = (bitBuffer >> count) & mask;
		bitBuffer >>= 8;
    }
}


/**
 * Converting commitments (list of matrices) into 8-byte list representation
 * @param uint8_t*     tab_bytes   :  List of commitments in 8-byte list returned
 * @param nmod_mat_t*  mat         :  List of commitments
 * @param uint8_t*     hash        :  Hash value to determine the # of commitment columns
 */
void responses_to_bytes(uint8_t *tab_bytes, const nmod_mat_t* mat, uint8_t* hash){
	 uint64_t bitBuffer = 0;  
     uint8_t count = 0;
	 int k = 0;
	 int size;
     for(int h = 0; h < 128; h++){
     	if((hash[h / 8] >> h % 8) & 1)
     	    size = R;
     	else
     	 	size = R2;
    	for(int i = 0; i < size; i++)
			for(int j = 0; j < M; j++){
				bitBuffer |= nmod_mat_entry(mat[h],i,j);
				count = count + 3;
				if(count >= 8){
					tab_bytes[k] = (bitBuffer >> count) & mask;
					bitBuffer >>= 8;
					count = count - 8;
					k++;
				}
				tab_bytes[k] = (bitBuffer >> count) & mask;
				bitBuffer >>= 8;
				k++;
			}
     }
		 
	if(count > 0){
    	tab_bytes[k] = (bitBuffer >> count) & mask;
		bitBuffer >>= 8;
		k++;
    }
}


/**
 * Write signature (commits and responses) into file
 * @param uint8_t*  commits_bytes    :  List of matrix K commited in 8-byte list
 * @param uint8_t*  responses_bytes  :  List of responses (X or XE) in 8-byte list
 * @param slong     response_size    :  Size of responses list
 */
void write_signature(uint8_t *commits_bytes, uint8_t *responses_bytes, slong response_size){
     FILE* sign = fopen("signature", "wb");
	 fwrite(commits_bytes, 1, COMMIT_SIZE, sign);
	 fwrite(responses_bytes, 1, response_size, sign);
	 fclose(sign);
}
