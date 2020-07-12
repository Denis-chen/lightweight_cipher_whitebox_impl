/*
 * @Author: Weijie Li 
 * @Date: 2018-01-23 23:35:52 
 * @Last Modified by: Weijie Li
 * @Last Modified time: 2018-01-29 22:52:37
 */


#ifndef __PRESENT80WB_H__
#define __PRESENT80WB_H__

#include <NTL/mat_GF2.h>
#include <affine/MatrixUtils.h>

/*********
 * 
 * present config 
 * 
 * ***************/

#define PRESENT_WB_DEBUG 0

// comment this out if this is used on PC
#define __UINT_T__

#ifndef __UINT_T__
#define __UINT_T__
typedef unsigned char uint8_t;
typedef unsigned short uint16_t;
typedef unsigned long uint32_t;
typedef unsigned long long uint64_t;
#endif /* __UINT_T__ */

#define PRESENT_ROUNDS 31

typedef struct _present_wb_ctx {
    uint8_t rk[PRESENT_ROUNDS+1][8][256];
    matrix_transform_t pLayer[PRESENT_ROUNDS][8];
    uint8_t stmp[8][256]; // start map

#if PRESENT_WB_DEBUG
	matrix_transform_t fc_inv[PRESENT_ROUNDS][2];
	matrix_transform_t f_inv[PRESENT_ROUNDS][8];
    matrix_transform_t g_inv[PRESENT_ROUNDS+1][8];
#endif //PRESENT_WB_DEBUG
    
    uint8_t rounds;
} present_wb_ctx;

// // full-round PRESENT block cipher
// #define present_wb(plain, key, cipher) present_wb_rounds((plain), (key), 31, (cipher))
// // full-round PRESENT whitebox block cipher init
// #define present_wb_init(key, ctx) present_wb_init_rounds((key), 31, (ctx))

// full-round should be 31, i.e. rounds = 31
void present_wb_init(const uint8_t *key, present_wb_ctx &ctx);

// recycling memory space
void present_wb_release(present_wb_ctx ctx);

// plain and cipher can overlap, so do key and cipher
void present_wb_enc(const uint8_t *plain, const present_wb_ctx &ctx, uint8_t *cipher);


#endif /* __PRESENT80WB_H__ */
