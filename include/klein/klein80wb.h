/*
 * @Author: Weijie Li 
 * @Date: 2017-11-17 00:50:46 
 * @Last Modified by: Weijie Li
 * @Last Modified time: 2017-11-19 11:35:34
 */
#ifndef KLEIN80WB_H_
#define KLEIN80WB_H_


// #include <stdint.h>

#include "kleinSbox.h"

#include <NTL/mat_GF2.h>
#include <affine/MatrixUtils.h>

//rounds of Klein-80;
#define ROUNDS_80 16

#define AFFINE_IMPL 
// #define AFFINE_IMPL_DEBUG  

// #define MAT_BIJECTIONS

typedef struct klein80_wb_t {
    int rounds;
    // uint8_t T[8][256];     // final xor
    uint32_t Te[ROUNDS_80][8][256];
    uint8_t SE[2][8][256]; //0 for initial encoding, 1 for final encoding
    

    #ifdef MAT_BIJECTIONS
    uint8_t fi[8][256],g[8][256];
    #endif //MAT_BIJECTIONS

    #ifdef AFFINE_IMPL_DEBUG 
        affine_transform_t xorInv[ROUNDS_80+1][2];
    #endif

} klein80_wb_t;

#define gen_klein80wb_table(key, ctx) gen_klein80wb_table_rounds(key, ctx, ROUNDS_80)

int gen_klein80wb_table_rounds(const uint8_t *key, klein80_wb_t * klein80wb_ctx, int rounds);

void klein80wb_encrypt(const uint8_t *plain, const klein80_wb_t * klein80wb_ctx,  uint8_t *cipher);

#endif /* KLEIN80WB_H_ */