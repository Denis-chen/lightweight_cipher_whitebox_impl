/*
 * @Author: Weijie Li 
 * @Date: 2017-11-08 19:35:45 
 * @Last Modified by: Weijie Li
 * @Last Modified time: 2017-11-26 22:00:24
 */

 #include <sms4wb/sms4_wb.h>
 #include "../sms4/sms4_lcl.h"

#ifdef DEBUG_INFO 
    #include <internal/DebugUtils.h>
#endif /* DEBUG_INFO */

  


#ifdef SMS4_ENCODING

    #ifdef DEBUG_INFO
        #define PRINT_INFO(x0, x1, x2, x3, x4, i) \
        printf("Round %d:\n", i); \
        xt0 = applyAffineToU32(sms4_wb_ctx->P[i][1], x0); \
        printf("\txt0\t%08X\n", xt0);\
        xt1 = applyAffineToU32(sms4_wb_ctx->P[i+1][1], x1); \
        printf("\txt1\t%08X\n", xt1); \    
        xt2 = applyAffineToU32(sms4_wb_ctx->P[i+2][1], x2); \
        printf("\txt2\t%08X\n", xt2); \
        xt3 = applyAffineToU32(sms4_wb_ctx->P[i+3][1], x3); \
        printf("\txt3\t%08X\n", xt3);\
        xt4 = applyAffineToU32(sms4_wb_ctx->P[i+4][1], x4); \    
        printf("\txt4\t%08X\n", xt4)
    #else 
        #define PRINT_INFO(x0, x1, x2, x3, x4, i)
    #endif /* DEBUG_INFO */

    #define SK32(A)						\
        ((skbox[0   +  (((A) >> 24) & 0xff)]) ^		\
        (skbox[256 +  (((A) >> 16) & 0xff)]) ^		\
        (skbox[512 +  (((A) >>  8) & 0xff)]) ^		\
        (skbox[768 +  (((A))       & 0xff)]      ))


    #define ROUND(x0, x1, x2, x3, x4, i)			\
    xt1 = applyAffineToU32(sms4_wb_ctx->M[i][0], x1); \
    xt2 = applyAffineToU32(sms4_wb_ctx->M[i][1], x2); \
    xt3 = applyAffineToU32(sms4_wb_ctx->M[i][2], x3); \
    x4 = xt1 ^ xt2 ^ xt3;			\
    x4 = SK32(x4);					        \
    skbox += 1024;                         \
    xt0 = applyAffineToU32(sms4_wb_ctx->C[i], x0); \
    xt4 = applyAffineToU32(sms4_wb_ctx->D[i], x4); \
    x4 = xt0 ^ xt4;\
    PRINT_INFO(x0, x1, x2, x3, x4, i)

#else 

    #define L32(x)						\
            ((x) ^						\
            ROT32((x),  2) ^				\
            ROT32((x), 10) ^				\
            ROT32((x), 18) ^				\
            ROT32((x), 24))

    #define SK32(A)						\
        ((skbox[0   +  (((A) >> 24) & 0xff)] << 24) ^		\
        (skbox[256 +  (((A) >> 16) & 0xff)] << 16) ^		\
        (skbox[512 +  (((A) >>  8) & 0xff)] <<  8) ^		\
        (skbox[768 +  (((A))       & 0xff)]      ))

    #define ROUND(x0, x1, x2, x3, x4, i)			\
    x4 = x1 ^ x2 ^ x3;			\
    x4 = SK32(x4);					        \
    skbox += 1024;                         \
    x4 = x0 ^ L32(x4)


#endif /* SMS4_ENCODING */


void sms4_wb_enc(const uint8_t *in, uint8_t *out, sms4_wb_t *sms4_wb_ctx) {
    // const uint32_t *rk = key->rk;

    uint32_t x0, x1, x2, x3, x4;
    
    uint32_t xt0, xt1, xt2, xt3, xt4;

	x0 = GET32(in     );
	x1 = GET32(in +  4);
	x2 = GET32(in +  8);
    x3 = GET32(in + 12);
    
    

    #ifdef SMS4_ENCODING
        uint32_t *skbox;
        x0 = applyAffineToU32(sms4_wb_ctx->SE[0], x0);
        x1 = applyAffineToU32(sms4_wb_ctx->SE[1], x1);
        x2 = applyAffineToU32(sms4_wb_ctx->SE[2], x2);
        x3 = applyAffineToU32(sms4_wb_ctx->SE[3], x3);

        skbox = (uint32_t*)sms4_wb_ctx->ssbox_enc;
    #else
        uint8_t *skbox;
        skbox = (uint8_t*)sms4_wb_ctx->ssbox_enc;
    #endif /* SMS4_ENCODING */


    ROUNDS(x0, x1, x2, x3, x4);   

    #ifdef SMS4_ENCODING
        x0 = applyAffineToU32(sms4_wb_ctx->FE[3], x0);
        x4 = applyAffineToU32(sms4_wb_ctx->FE[2], x4);
        x3 = applyAffineToU32(sms4_wb_ctx->FE[1], x3);
        x2 = applyAffineToU32(sms4_wb_ctx->FE[0], x2);
    #endif /* SMS4_ENCODING */

    PUT32(x0, out     );
	PUT32(x4, out +  4);
	PUT32(x3, out +  8);
	PUT32(x2, out + 12);

	x0 = x1 = x2 = x3 = x4 = 0;

    // DUMPPINFO(out, 16, "encrypt");
}

