/*
 * @Author: Weijie Li 
 * @Date: 2017-11-08 16:39:13 
 * @Last Modified by: Weijie Li
 * @Last Modified time: 2017-11-26 23:48:04
 */

#include <sms4wb/sms4_wb_generator.h>
#include <affine/MatrixUtils.h>
#include <affine/NTLUtils.h>

#include <assert.h>



#include <NTL/mat_GF2.h>
#include <NTL/vec_GF2.h>


#include <internal/DebugUtils.h>

extern uint8_t SBOX[256];



typedef struct sms4_table_assistant_t{
    affine_transform_t P[SMS4_NUM_STATES][2];
    affine_transform_t E[SMS4_NUM_STATES][5][2];
    // affine_transform_t M[SMS4_NUM_STATES];
    affine_transform_t Q[SMS4_NUM_STATES][2];

    uint8_t skbox_enc[SMS4_NUM_ROUNDS+1][4][256];

} sms4_table_assistant_t;

long sms4_csl_xor_matrix[32] = {0xA0202080, 0x50101040, 0x28080820, 0x14040410, 0xA020208, 0x5010104, 0x2808082, 0x1404041, 0x80A02020, 0x40501010, 0x20280808, 0x10140404, 0x80A0202, 0x4050101, 0x82028080, 0x41014040, 0x2080A020, 0x10405010, 0x8202808, 0x4101404, 0x2080A02, 0x1040501, 0x80820280, 0x40410140, 0x202080A0, 0x10104050, 0x8082028, 0x4041014, 0x202080A, 0x1010405, 0x80808202, 0x40404101};

void sms4_wb_gen_secrect_sbox(sms4_key_t *sms4_key, sms4_table_assistant_t *assistant) {
    uint8_t *sst;
    uint32_t *rk = sms4_key->rk;
    sst = (uint8_t*)assistant->skbox_enc;
    
    for (int i=0; i< SMS4_NUM_ROUNDS; i++) {
        for (int j=0; j<256; j++) {
            sst[   0 + j ] = SBOX[ j ^ ((rk[i] >> 24) & 0xff) ];
            sst[ 256 + j ] = SBOX[ j ^ ((rk[i] >> 16) & 0xff) ];
            sst[ 512 + j ] = SBOX[ j ^ ((rk[i] >>  8) & 0xff) ];
            sst[ 768 + j ] = SBOX[ j ^ ((rk[i]      ) & 0xff) ];
        }
        sst += 1024;
    }
}


#ifdef  SMS4_ENCODING


#define RANDOM_AFFINE_MAT(x, xi, y, yi, d)   genRandomAffineMatrix(x, xi, y, d);\
        yi = xi*y


#define RANDOM_DEBUG_AFFINE_MAT(x, xi, y, yi, d)   genIndMatrix(x, d);    \
        genIndMatrix(xi, d);                                         \
        genZeroVec(y, d);                                   \
        genZeroVec(yi, d)


void sms4_wb_gen_affine(sms4_wb_t *sms4_wb_ctx, sms4_table_assistant_t *assistant) {
    int i,j;
    for (i=0; i<SMS4_NUM_STATES; i++) {
        //gen P affine matrix
        RANDOM_AFFINE_MAT(assistant->P[i][0].linearMap, assistant->P[i][1].linearMap, assistant->P[i][0].vectorTranslation, assistant->P[i][1].vectorTranslation, 32);

        //gen E affine matrix
        for (j=0; j<4; j++) {
            RANDOM_AFFINE_MAT(assistant->E[i][j][0].linearMap, assistant->E[i][j][1].linearMap, assistant->E[i][j][0].vectorTranslation, assistant->E[i][j][1].vectorTranslation, 8);

        }
        NTL::mat_GF2 T,S;
        combineDiagMat(T, assistant->E[i][0][0].linearMap, assistant->E[i][1][0].linearMap);
        combineDiagMat(S, T, assistant->E[i][2][0].linearMap);
        combineDiagMat(assistant->E[i][4][0].linearMap, S, assistant->E[i][3][0].linearMap);
        inv(assistant->E[i][4][1].linearMap, assistant->E[i][4][0].linearMap);
        assistant->E[i][4][0].vectorTranslation.SetLength(32);
        for (j=0; j<32; j++) {
            assistant->E[i][4][0].vectorTranslation.put(j, assistant->E[i][j/8][0].vectorTranslation.get(j%8));
        }
        assistant->E[i][4][1].vectorTranslation = assistant->E[i][4][1].linearMap * assistant->E[i][4][0].vectorTranslation;

        //gen Q affine matrix
        RANDOM_AFFINE_MAT(assistant->Q[i][0].linearMap, assistant->Q[i][1].linearMap, assistant->Q[i][0].vectorTranslation, assistant->Q[i][1].vectorTranslation, 32);
    }

    
}

//gen CSL matrix
// x + (x<<<2) + (x<<<10) + (x<<<18) + (x<<<24)
void sms4_wb_gen_csl_matrix(mat_GF2 &csl_matrix) {
    csl_matrix.SetDims(32, 32);
    initMatrixFromBit(csl_matrix, sms4_csl_xor_matrix);
}

int sms4_wb_table_check(sms4_table_assistant_t *assistant, mat_GF2 &csl_matrix) {
    mat_GF2 m32a, m32b;
    vec_GF2 v32a, v32b, v8a, v8b;
    random(m32a, 32, 32);
    random(v32a, 32);
    random(v8a, 8);
    int i,j;
    for (i=0; i<SMS4_NUM_STATES; i++) {
        
        v32b = assistant->P[i][0].linearMap * v32a + assistant->P[i][0].vectorTranslation;
        v32b = assistant->P[i][1].linearMap * (v32b + assistant->P[i][0].vectorTranslation);
        assert(v32a==v32b);

        v32b = assistant->P[i][0].linearMap * v32a + assistant->P[i][0].vectorTranslation;
        v32b = assistant->P[i][1].linearMap * v32b + assistant->P[i][1].vectorTranslation;
        assert(v32a==v32b);

        v32b = assistant->Q[i][0].linearMap * v32a + assistant->Q[i][0].vectorTranslation;
        v32b = assistant->Q[i][1].linearMap * (v32b + assistant->Q[i][0].vectorTranslation);
        assert(v32a==v32b);

        v32b = assistant->Q[i][0].linearMap * v32a + assistant->Q[i][0].vectorTranslation;
        v32b = assistant->Q[i][1].linearMap * v32b + assistant->Q[i][1].vectorTranslation;
        assert(v32a==v32b);

        for (j=0; j<4; j++) {
            v8b = assistant->E[i][j][0].linearMap * v8a + assistant->E[i][j][0].vectorTranslation;
            v8b = assistant->E[i][j][1].linearMap * (v8b + assistant->E[i][j][0].vectorTranslation);
            assert(v8a==v8b);

            v8b = assistant->E[i][j][0].linearMap * v8a + assistant->E[i][j][0].vectorTranslation;
            v8b = assistant->E[i][j][1].linearMap * v8b + assistant->E[i][j][1].vectorTranslation;
            assert(v8a==v8b);
        }
        
        v32b = assistant->E[i][4][0].linearMap * v32a + assistant->E[i][4][0].vectorTranslation;
        v32b = assistant->E[i][4][1].linearMap * (v32b + assistant->E[i][4][0].vectorTranslation);
        assert(v32a==v32b);

        v32b = assistant->E[i][4][0].linearMap * v32a + assistant->E[i][4][0].vectorTranslation;
        v32b = assistant->E[i][4][1].linearMap * v32b + assistant->E[i][4][1].vectorTranslation;
        assert(v32a==v32b);
    }
    return 0;
}

int sms4_combine_affine_table(sms4_wb_t *sms4_wb_ctx, sms4_table_assistant_t *assistant, mat_GF2 &csl_matrix) {
    int i,j;
    for (i=0; i<SMS4_NUM_ROUNDS; i++) {

        //part 1. gen M affine matrix
        sms4_wb_ctx->M[i][0].linearMap = assistant->E[i][4][1].linearMap * assistant->P[i + 1][1].linearMap;
        sms4_wb_ctx->M[i][0].vectorTranslation =    assistant->E[i][4][1].linearMap * (assistant->P[i + 1][1].vectorTranslation) + \
                                                    assistant->E[i][4][1].vectorTranslation ;

        sms4_wb_ctx->M[i][1].linearMap = assistant->E[i][4][1].linearMap * assistant->P[i + 2][1].linearMap;
        sms4_wb_ctx->M[i][1].vectorTranslation =    assistant->E[i][4][1].linearMap * (assistant->P[i + 2][1].vectorTranslation) + \
                                                    assistant->E[i][4][1].vectorTranslation ;

        sms4_wb_ctx->M[i][2].linearMap = assistant->E[i][4][1].linearMap * assistant->P[i + 3][1].linearMap;
        sms4_wb_ctx->M[i][2].vectorTranslation =    assistant->E[i][4][1].linearMap * (assistant->P[i + 3][1].vectorTranslation) + \
                                                    assistant->E[i][4][1].vectorTranslation ;


        //part 2. gen Q combine L into 4 matrix
        mat_GF2 QL;

        QL = assistant->Q[i][0].linearMap * csl_matrix;
        mat_GF2 QLi[4];
        for (j=0; j<4; j++) {
            QLi[j].SetDims(32, 8);
            for (int ii=0; ii<32; ii++) {
                for (int jj=0; jj<8; jj++) {
                    QLi[j].put(ii, jj, QL.get(ii, j*8+jj));
                }
            }
        }

        int k,d;
        uint32_t r = get32FromVec(assistant->Q[i][0].vectorTranslation);

        for (k=0; k<256; k++) {

            for (int d=0; d<4; d++) {
                vec_GF2 kv;
                int kd; 
                initVecFromBit(kv, k, 8);
                kv = assistant->E[i][d][0].linearMap * kv + assistant->E[i][d][0].vectorTranslation;
                kd = assistant->skbox_enc[i][d][getDigitalFromVec(kv)];
                initVecFromBit(kv, kd, 8);
                
                vec_GF2 dest;
                dest = QLi[d]*kv;
                uint32_t result = get32FromVec(dest);
                
                sms4_wb_ctx->ssbox_enc[i][d][k] = result;
            }
            sms4_wb_ctx->ssbox_enc[i][3][k] = sms4_wb_ctx->ssbox_enc[i][3][k] ^ r;
            
        }

        //part 3. gen C D matrix, C for Xi0, D for T(Xi1+Xi2+Xi3+rk)
        
        sms4_wb_ctx->C[i].linearMap =   assistant->P[i+4][0].linearMap * \
                                        assistant->P[i][1].linearMap;
        sms4_wb_ctx->C[i].vectorTranslation =   assistant->P[i+4][0].vectorTranslation + \
                                                assistant->P[i+4][0].linearMap * \
                                                assistant->P[i][1].vectorTranslation;

        sms4_wb_ctx->D[i].linearMap =   assistant->P[i+4][0].linearMap * \
                                        assistant->Q[i][1].linearMap;
        sms4_wb_ctx->D[i].vectorTranslation =   assistant->P[i+4][0].vectorTranslation + \
                                                assistant->P[i+4][0].linearMap * \
                                                assistant->Q[i][1].vectorTranslation;

        //Important: A super pit that wasted my three days
        sms4_wb_ctx->D[i].vectorTranslation =   assistant->P[i+4][0].vectorTranslation + \
                                                sms4_wb_ctx->D[i].vectorTranslation;

    }


    //external encoding

    for (int i=0; i<4; i++) {
        sms4_wb_ctx->SE[i].linearMap = assistant->P[i][0].linearMap;
        sms4_wb_ctx->SE[i].vectorTranslation = assistant->P[i][0].vectorTranslation;

        sms4_wb_ctx->FE[i].linearMap = assistant->P[SMS4_NUM_ROUNDS+i][1].linearMap;
        sms4_wb_ctx->FE[i].vectorTranslation = assistant->P[SMS4_NUM_ROUNDS+i][1].vectorTranslation;
    }

    


    return 0;
}
 
#endif /* SMS4_ENCODING  */



int sms4_wb_gen_tables(const uint8_t *key, sms4_wb_t *sms4_wb_ctx, int enc) {
    sms4_key_t sms4_enc_key;
    int ret = 0;
    
    if (enc==SMS4_ENCRYPT)
        sms4_set_encrypt_key(&sms4_enc_key, key);
    else 
        sms4_set_decrypt_key(&sms4_enc_key, key);

    sms4_table_assistant_t assistant;

    sms4_wb_gen_secrect_sbox(&sms4_enc_key, &assistant);

    // DUMP(sms4_wb_ctx->skbox_enc, sizeof (sms4_wb_ctx->skbox_enc));
    
    #ifdef SMS4_ENCODING

        sms4_wb_gen_affine(sms4_wb_ctx, &assistant);

        mat_GF2 csl_matrix;
        sms4_wb_gen_csl_matrix(csl_matrix);
    
        

        if (sms4_wb_table_check(&assistant, csl_matrix)!=0)
            return 1;

        ret = sms4_combine_affine_table(sms4_wb_ctx, &assistant, csl_matrix);

        #ifdef DEBUG_INFO
        // memcpy(sms4_wb_ctx->P, assistant.P, SMS4_NUM_STATES*2*sizeof(affine_transform_t));
        for (int i=0; i<SMS4_NUM_STATES; i++) {
            sms4_wb_ctx->P[i][0] = assistant.P[i][0];
            sms4_wb_ctx->P[i][1] = assistant.P[i][1];
        }
        #endif /*  DEBUG_INFO */   


    #else 

        memcpy(sms4_wb_ctx->ssbox_enc, assistant.skbox_enc, sizeof(assistant.skbox_enc));
        
    #endif /* SMS4_ENCODING */


    
    return ret;

}
