/*
 * @Author: Weijie Li 
 * @Date: 2017-11-06 22:04:17 
 * @Last Modified by: Weijie Li
 * @Last Modified time: 2017-11-26 22:04:41
 */


#ifndef  SMS4_WB_H
#define  SMS4_WB_H


#include <stdint.h>
#include <sms4/sms4.h>

#include <affine/MatrixUtils.h>


#define SMS4_ENCODING
// #define DEBUG_INFO


#define SMS4_ENCRYPT 1
#define SMS4_DECRYPT 0

#define SMS4_NUM_STATES (SMS4_NUM_ROUNDS+5)

typedef struct sms4_wb_t{
    uint32_t rk;
    // affine_transform_t Q[SMS4_NUM_STATES];

#ifdef SMS4_ENCODING 
    uint32_t    ssbox_enc[SMS4_NUM_ROUNDS+1][4][256];
    affine_transform_t M[SMS4_NUM_STATES][3];
    affine_transform_t C[SMS4_NUM_ROUNDS];
    affine_transform_t D[SMS4_NUM_ROUNDS];
    // start encoding
    affine_transform_t SE[4];
    affine_transform_t FE[4];
#else 
    uint8_t     ssbox_enc[SMS4_NUM_ROUNDS+1][4][256];
#endif /* SMS4_ENCODING */


    #ifdef DEBUG_INFO
    //debug to be del
    affine_transform_t P[SMS4_NUM_STATES][2];
    #endif /*  DEBUG_INFO */   

} sms4_wb_t;



void sms4_wb_enc(const uint8_t *in, uint8_t *out, sms4_wb_t *sms4_wb_ctx);


#endif //SMS4_WB_H