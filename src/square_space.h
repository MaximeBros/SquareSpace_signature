#ifndef SQUARE_SPACE_H
#define SQUARE_SPACE_H

#include <flint/fq.h>
#include <stdlib.h>
#include <stdint.h>
#include "parameters.h"

void random_fqm_list(fq_t* array_poly, fq_ctx_t ctx_fpm, fmpz_mod_ctx_t ctx_fp);
void square(fq_t* array_square, const fq_t* array_elt, slong nb_elts, fq_ctx_t ctx_fpm);
void multiply(fq_t* array_res, const fq_t* array_elt_1, const fq_t* array_elt_2, slong nb_elts, fq_ctx_t ctx_fpm);
void generate_hash(uint8_t* *hash, const nmod_mat_t *commits, const nmod_mat_t mat_U);
u_int32_t arc4random(void);

#endif //SQUARE_SPACE_H
