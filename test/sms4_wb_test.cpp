extern "C" {

#include "utils/testutil.h"

}

#include <sms4/sms4.h>
#include <sms4wb/sms4_wb_generator.h>

#include<internal/DebugUtils.h>


//#include <NTL/mat.h>


typedef struct {
    size_t size;
    const unsigned char *data;
}  SIZED_DATA;

typedef struct {
    size_t size;
    unsigned char *data;
}  MUTABLE_SIZED_DATA;

#define u8 unsigned char

/**********************************************************************
 *
 * Test of cbc128
 * void sms4_cbc_encrypt(const unsigned char *in, unsigned char *out,
    size_t len, const sms4_key_t *key, unsigned char *iv, int enc);
 ***/

/* Test Case 1 */
static const u8 SMS4_CBC_K1[16] = {0};
static const u8 SMS4_CBC_P1[] = { 0 };
static const u8 SMS4_CBC_IV1[12] = {0} ;
static const u8 SMS4_CBC_C1[] = { 0 };

/* Test Case 2 */
static const u8 SMS4_CBC_K2[16] = {
	0x2B, 0x7E, 0x15, 0x16, 0x28, 0xAE, 0xD2, 0xA6,
	0xAB, 0xF7, 0x15, 0x88, 0x09, 0xCF, 0x4F, 0x3C
};

static const u8 SMS4_CBC_IV2[] = {
    0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07,
 	0x08, 0x09, 0x0A, 0x0B, 0x0C, 0x0D, 0x0E, 0x0F
};

static const u8 SMS4_CBC_P2[] = {
	0x6B, 0xC1, 0xBE, 0xE2, 0x2E, 0x40, 0x9F, 0x96,
 	0xE9, 0x3D, 0x7E, 0x11, 0x73, 0x93, 0x17, 0x2A
};

static const u8 SMS4_CBC_C2[] = {
	0x76, 0x49, 0xAB, 0xAC, 0x81, 0x19, 0xB2, 0x46,
 	0xCE, 0xE9, 0x8E, 0x9B, 0x12, 0xE9, 0x19, 0x7D
};

/* Test Case 3 */
# define SMS4_CBC_K3 SMS4_CBC_K2
# define SMS4_CBC_IV3 SMS4_CBC_C2
static const u8 SMS4_CBC_P3[] = {
	0xAE, 0x2D, 0x8A, 0x57, 0x1E, 0x03, 0xAC, 0x9C,
 	0x9E, 0xB7, 0x6F, 0xAC, 0x45, 0xAF, 0x8E, 0x51
 };

static const u8 SMS4_CBC_C3[] = {
	0x50, 0x86, 0xCB, 0x9B, 0x50, 0x72, 0x19, 0xEE,
 	0x95, 0xDB, 0x11, 0x3A, 0x91, 0x76, 0x78, 0xB2
};

/* Test Case 4 */
static const u8 SMS4_CBC_K4[16] = {
    0x01, 0x23, 0x45, 0x67, 0x89, 0xab, 0xcd, 0xef,
    0xfe, 0xdc, 0xba, 0x98, 0x76, 0x54, 0x32, 0x10
};

# define SMS4_CBC_IV4 SMS4_CBC_C2

static const u8 SMS4_CBC_P4[] = {
	0x01, 0x23, 0x45, 0x67, 0x89, 0xab, 0xcd, 0xef,
    0xfe, 0xdc, 0xba, 0x98, 0x76, 0x54, 0x32, 0x10
 };

static const u8 SMS4_CBC_C4[] = {
    0x68, 0x1E, 0xDF, 0x34, 0xD2, 0x06, 0x96, 0x5E,
    0x86, 0xB3, 0xE9, 0x4F, 0x53, 0x6E, 0x42, 0x46
};




#define CBC128_TEST_VECTOR(n)                   \
    {                                           \
        {sizeof(SMS4_CBC_K##n), SMS4_CBC_K##n},                   \
        {sizeof(SMS4_CBC_IV##n), SMS4_CBC_IV##n},                 \
        {sizeof(SMS4_CBC_P##n), SMS4_CBC_P##n},                   \
        {sizeof(SMS4_CBC_C##n), SMS4_CBC_C##n},                   \
    }
static struct cbc128_data {
    const SIZED_DATA K;
    const SIZED_DATA IV;
    const SIZED_DATA P;
    const SIZED_DATA C;
} cbc128_vectors[] = {
    CBC128_TEST_VECTOR(1),
    CBC128_TEST_VECTOR(2),
    CBC128_TEST_VECTOR(3),
    CBC128_TEST_VECTOR(4)    
};

/* only check the encrypted date can decrypt correct*/
int test_sms4_wb128(int idx)
{
    unsigned char out[512];
    SIZED_DATA K = cbc128_vectors[idx].K;
    SIZED_DATA C_IV = cbc128_vectors[idx].IV;
    SIZED_DATA P = cbc128_vectors[idx].P;
    sms4_key_t key;
    sms4_wb_t sms4_wb_enc_ctx, sms4_wb_dec_ctx;

    MUTABLE_SIZED_DATA C;
    MUTABLE_SIZED_DATA IV;

    /* Size 1 inputs are special-cased to signal NULL. */
    if (P.size == 1)
        P.data = NULL;

    memset(out, 0, P.size);
    IV.size = C_IV.size;
    IV.data = (uint8_t*)malloc(IV.size);
    memcpy(IV.data, C_IV.data, IV.size);
    
    sms4_wb_gen_tables(K.data, &sms4_wb_enc_ctx, SMS4_ENCRYPT);
    sms4_wb_gen_tables(K.data, &sms4_wb_dec_ctx, SMS4_DECRYPT);
//    sms4_set_encrypt_key(&key, K.data);

    if (P.data != NULL)
        sms4_wb_enc(P.data, out, &sms4_wb_enc_ctx);
//sms4_encrypt(P.data, out, &key);

    // if (C.data != NULL
    //     	&& !TEST_mem_eq(out, P.size, C.data, P.size))
    //     return 0;

    // if (P.data != NULL)
    //     DUMPPINFO(P.data, 16, "p.data");
    // DUMPPINFO(out, 16, "cipher");

    IV.size = C_IV.size;
    memcpy(IV.data, C_IV.data, IV.size);
    // sms4_set_decrypt_key(&key, K.data);

    if (P.data != NULL) {
        C.size = P.size;
        C.data = (uint8_t*)malloc(C.size);
        memcpy(C.data, out, C.size);
        memset(out, 0, P.size);
//        sms4_decrypt(C.data, out, &key);
        sms4_wb_enc(C.data, out, &sms4_wb_dec_ctx);
        // DUMPPINFO(out, 16, "decrypt");
    }

    if (P.data != NULL
            && !TEST_mem_eq(out, P.size, P.data, P.size))
        return 1;

    return 0;
}
