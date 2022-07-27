#ifndef FLINT_UTILS_H
#define FLINT_UTILS_H

#include <flint/fq.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include "parameters.h"

void init_array(fq_t* *array, slong size, fq_ctx_t ctx_fqm);
void clear_array(fq_t* *array, slong size, fq_ctx_t ctx_fqm);
void init_modulus(fmpz_mod_poly_t* modulus, fmpz_mod_ctx_t* ctx_fq, fq_ctx_t* ctx_fqm, fmpz_t p);
void clear_vars(fq_t* *array_E, fq_t* *array_U, fq_t* *array_X, fq_t* *array_K, fq_t* *array_XE, fmpz_mod_poly_t* modulus, fmpz_mod_ctx_t* ctx_fq, fq_ctx_t* ctx_fqm, fmpz_t* p, uint8_t* *hash, uint8_t* *commits_bytes, uint8_t* *mat_U_bytes, uint8_t* *message_bytes);

long fqm_list_to_mat(nmod_mat_t mat, fmpz_t p, const fq_t* fqm_list, slong length, fq_ctx_t ctx_fqm);
void mat_to_fqm_list(const nmod_mat_t mat, fq_t* fqm_list, slong length, fq_ctx_t ctx_fqm, fmpz_mod_ctx_t ctx_fq);

void list_mat_to_bytes(uint8_t *tab_bytes, const nmod_mat_t *mat, slong size);
void mat_to_bytes(uint8_t *tab_bytes, const nmod_mat_t mat, slong size);
void responses_to_bytes(uint8_t *tab_bytes, const nmod_mat_t* mat, uint8_t* hash);
void write_signature(uint8_t *commits_bytes, uint8_t *responses_bytes, slong size_responses);

#endif //FLINT_UTILS_H
