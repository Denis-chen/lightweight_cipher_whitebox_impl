/* ====================================================================
 * Copyright (c) 2014 - 2017 The GmSSL Project.  All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in
 *    the documentation and/or other materials provided with the
 *    distribution.
 *
 * 3. All advertising materials mentioning features or use of this
 *    software must display the following acknowledgment:
 *    "This product includes software developed by the GmSSL Project.
 *    (http://gmssl.org/)"
 *
 * 4. The name "GmSSL Project" must not be used to endorse or promote
 *    products derived from this software without prior written
 *    permission. For written permission, please contact
 *    guanzhi1980@gmail.com.
 *
 * 5. Products derived from this software may not be called "GmSSL"
 *    nor may "GmSSL" appear in their names without prior written
 *    permission of the GmSSL Project.
 *
 * 6. Redistributions of any form whatsoever must retain the following
 *    acknowledgment:
 *    "This product includes software developed by the GmSSL Project
 *    (http://gmssl.org/)"
 *
 * THIS SOFTWARE IS PROVIDED BY THE GmSSL PROJECT ``AS IS'' AND ANY
 * EXPRESSED OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE GmSSL PROJECT OR
 * ITS CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
 * NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
 * STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 * OF THE POSSIBILITY OF SUCH DAMAGE.
 * ====================================================================
 */

#ifndef HEADER_SMS4_H
#define HEADER_SMS4_H

// #include <openssl/opensslconf.h>
#ifndef NO_GMSSL

# define SMS4_ENCRYPT     1
# define SMS4_DECRYPT     0

#define SMS4_KEY_LENGTH		16
#define SMS4_BLOCK_SIZE		16
#define SMS4_IV_LENGTH		(SMS4_BLOCK_SIZE)
#define SMS4_NUM_ROUNDS		32

#include <sys/types.h>
#include <stdint.h>
#include <string.h>
#include <openssl/modes.h>
#include <stdlib.h>
#include <stdio.h>

#define SIMPLE_EXTERNAL_ENCODINGS

#ifdef __cplusplus
extern "C" {
#endif

typedef struct sms4_key_t {
	uint32_t rk[SMS4_NUM_ROUNDS];
#ifdef SIMPLE_EXTERNAL_ENCODINGS
    uint8_t G[SMS4_NUM_ROUNDS+1][256];
    uint8_t F[SMS4_NUM_ROUNDS+1][256];
#endif //SIMPLE_EXTERNAL_ENCODINGS
} sms4_key_t;

typedef struct sms4_key_t SMS4_KEY;

/**
 sms4_set_encrypt_key

 @param key SMS4_KEY
 @param user_key key
 */
void sms4_set_encrypt_key(sms4_key_t *key, const unsigned char *user_key);

/**
 sms4_set_decrypt_key

 @param key SMS4_KEY
 @param user_key key
 */
void sms4_set_decrypt_key(sms4_key_t *key, const unsigned char *user_key);

/**
 sms4_encrypt

 @param in in
 @param out out
 @param key SMS4_KEY
 */
void sms4_encrypt(const unsigned char *in, unsigned char *out, const sms4_key_t *key);
#define sms4_decrypt(in,out,key)  sms4_encrypt(in,out,key)

/**
 sms4_encrypt_init

 @param key SMS4_KEY
 */
void sms4_encrypt_init(sms4_key_t *key);
void sms4_encrypt_8blocks(const unsigned char *in, unsigned char *out, const sms4_key_t *key);
void sms4_encrypt_16blocks(const unsigned char *in, unsigned char *out, const sms4_key_t *key);

/**
 sms4_ecb_encrypt

 @param in in
 @param out out
 @param key key
 @param enc 1 to SMS4_ENCRYPT, 0 to SMS4_DECRYPT
 */
void sms4_ecb_encrypt(const unsigned char *in, unsigned char *out,
	const sms4_key_t *key, int enc);

/**
 sms4_cbc_encrypt

 @param in in
 @param out out
 @param len byte size of in
 @param key key
 @param iv iv
 @param enc 1 to SMS4_ENCRYPT, 0 to SMS4_DECRYPT
 */
void sms4_cbc_encrypt(const unsigned char *in, unsigned char *out,
	size_t len, const sms4_key_t *key, unsigned char *iv, int enc);
void sms4_cfb128_encrypt(const unsigned char *in, unsigned char *out,
	size_t len, const sms4_key_t *key, unsigned char *iv, int *num, int enc);
void sms4_ofb128_encrypt(const unsigned char *in, unsigned char *out,
	size_t len, const sms4_key_t *key, unsigned char *iv, int *num);
/**
 sms4_ctr128_encrypt
 The input encrypted as though 128bit counter mode is being used.  The
 extra state information to record how much of the 128bit block we have
 used is contained in *num, and the encrypted counter is kept in
 ecount_buf.  Both *num and ecount_buf must be initialised with zeros
 before the first call to sms4_ctr128_encrypt(). This algorithm assumes
 that the counter is in the x lower bits of the IV (ivec), and that the
 application has full control over overflow and the rest of the IV.  This
 implementation takes NO responsibility for checking that the counter
 doesn't overflow into the rest of the IV when incremented.

 @param in in
 @param out out
 @param len byte size of in
 @param key key
 @param iv iv
 @param ecount_buf extra state, must be initialised with zeros before the first call
 @param num extra state, must be initialised with zeros before the first call
 */
void sms4_ctr128_encrypt(const unsigned char *in, unsigned char *out,
	size_t len, const sms4_key_t *key, unsigned char *iv,
	unsigned char ecount_buf[SMS4_BLOCK_SIZE], unsigned int *num);

/**
 sms4_wrap_key

 @param key key
 @param iv iv
 @param out out
 @param in in
 @param inlen byte size of in
 @return 1 to successful, otherwise fault
 */
int sms4_wrap_key(sms4_key_t *key, const unsigned char *iv,
	unsigned char *out, const unsigned char *in, unsigned int inlen);
/**
 sms4_unwrap_key

 @param key key
 @param iv iv
 @param out out
 @param in in
 @param inlen byte size of in
 @return 1 to successful, otherwise fault
 */
int sms4_unwrap_key(sms4_key_t *key, const unsigned char *iv,
	unsigned char *out, const unsigned char *in, unsigned int inlen);

typedef GCM128_CONTEXT SMS4_GCM128_CONTEXT;

/**
 sms4_gcm128_init

 @param ctx SMS4_GCM128_CONTEXT
 @param key SMS4_KEY
 */
void sms4_gcm128_init(SMS4_GCM128_CONTEXT *ctx, SMS4_KEY *key);

/**
 sms4_gcm128_setiv

 @param ctx SMS4_GCM128_CONTEXT
 @param ivec iv
 @param len byte size of iv
 */
void sms4_gcm128_setiv(SMS4_GCM128_CONTEXT *ctx, const unsigned char *ivec,
                      size_t len);

/**
 addition message of gcm

 @param ctx SMS4_GCM128_CONTEXT
 @param aad addition message
 @param len byte size of aad
 @return 1 to successful, otherwises fault
 */
int sms4_gcm128_aad(SMS4_GCM128_CONTEXT *ctx, const unsigned char *aad,
                    size_t len);

/**
 sms4_gcm128_encrypt

 @param in in
 @param out out
 @param length byte size of in
 @param ctx SMS4_GCM128_CONTEXT
 @param enc 1 to SMS4_ENCRYPT, 0 to SMS4_DECRYPT
 @return 1 to successful, otherwises fault
 */
int sms4_gcm128_encrypt(const unsigned char *in, unsigned char *out,
                        size_t length, SMS4_GCM128_CONTEXT *ctx, const int enc);

/**
 get tag of sms4_gcm128

 @param ctx SMS4_GCM128_CONTEXT
 @param tag memory for storage tag
 @param len byte size of tag
 */
void sms4_gcm128_tag(SMS4_GCM128_CONTEXT *ctx, unsigned char *tag,
                    size_t len);

/**
 sms4_gcm128_finish

 @param ctx SMS4_GCM128_CONTEXT
 @param tag memory for storage tag
 @param len byte size of tag
 @return 1 to successful, otherwises fault
 */
int sms4_gcm128_finish(SMS4_GCM128_CONTEXT *ctx, const unsigned char *tag,
                      size_t len);

/**
 release SMS4_GCM128_CONTEXT

 @param ctx SMS4_GCM128_CONTEXT
 */
void sms4_gcm128_release(SMS4_GCM128_CONTEXT *ctx);

/*
void sms4_avx2_encrypt_init(sms4_key_t *key);
void sms4_avx2_encrypt_8blocks(const unsigned char *in, unsigned char *out, const sms4_key_t *key);
void sms4_avx2_encrypt_16blocks(const unsigned char *in, unsigned char *out, const sms4_key_t *key);

void sms4_knc_encrypt_init(sms4_key_t *key);
void sms4_knc_encrypt_8blocks(const unsigned char *in, unsigned char *out, const sms4_key_t *key);
void sms4_knc_encrypt_16blocks(const unsigned char *in, unsigned char *out, const sms4_key_t *key);

#define SMS4_EDE_KEY_LENGTH	32

typedef struct {
	sms4_key_t k1;
	sms4_key_t k2;
} sms4_ede_key_t;

void sms4_ede_set_encrypt_key(sms4_ede_key_t *key, const unsigned char *user_key);
void sms4_ede_set_decrypt_key(sms4_ede_key_t *key, const unsigned char *user_key);
void sms4_ede_encrypt(sms4_ede_key_t *key, const unsigned char *in, unsigned char *out);
void sms4_ede_encrypt_8blocks(sms4_ede_key_t *key, const unsigned char *in, unsigned char *out);
void sms4_ede_encrypt_16blocks(sms4_ede_key_t *key, const unsigned char *in, unsigned char *out);
void sms4_ede_decrypt(sms4_ede_key_t *key, const unsigned char *in, unsigned char *out);
void sms4_ede_decrypt_8blocks(sms4_ede_key_t *key, const unsigned char *in, unsigned char *out);
void sms4_ede_decrypt_16blocks(sms4_ede_key_t *key, const unsigned char *in, unsigned char *out);
*/

/*
*  Created on: 2017.10.18
*      Author: lzj
*
*  encrypted file format:
*
*  MAGIC_TAG  : AYCF-SM4-GCM
*  VERSION	  : 1
*  TAG		    : cipher tag
*  SIZE	      : plaintext size
*  CONTENT	  : cipher text
*
*	how to use?  see the test function under below .check out!
*/
#define SMS4_GCM_FILE_MAGIC_TAG "AYCF-SM4-GCM"
#define SMS4_GCM_FILE_MAGIC_TAG_LEN 12
#define SMS4_GCM_FILE_VERSION 1
#define SMS4_GCM_FILE_TAG_LEN 16
#define SM4_GCM_FILE_MAX_BLOCK_LEN 128

/**
 * sm4 gcm file context
 */
typedef struct {
		SMS4_GCM128_CONTEXT *gcm;
		unsigned char tag[SMS4_GCM_FILE_TAG_LEN];
} sm4_gcmf_context;

/**
 * init the sm4 gcm file context
 *
 * @param  ctx [in]		gcm file context
 *
 * @return     [flag]		if successful o,otherwise failed
 */
int sm4_gcmf_init(sm4_gcmf_context *ctx);

/**
 * gcm file context free
 *
 * @param  ctx [in]		gcm file context
 *
 * @return     [flag]		if successful o,otherwise failed
 */
int sm4_gcmf_free(sm4_gcmf_context *ctx);

/**
 * set sm4 iv param
 *
 * @param  ctx [in]		gcm file context
 *
 * @param  iv  [iv]		iv array
 *
 * @param  len [in]		iv array length
 *
 * @return     [flag]		if successful o,otherwise failed
 */
int sm4_gcmf_set_iv(sm4_gcmf_context *ctx, const unsigned char * iv, size_t len);

/**
 * set sm4 key
 *
 * @param  ctx [in]		gcm file context
 *
 * @param  key [in]		key array
 *
 * @param  len [in]		key array length
 *
 * @return     [flag]		if successful o,otherwise failed
 */
int sm4_gcmf_set_key(sm4_gcmf_context *ctx, const unsigned char * key, size_t len);


/**
 * encrypte file
 *
 * @param  ctx      [in]		gcm file context
 *
 * @param  infpath  [in]		plaintext file input path
 *
 * @param  outfpath [in]		cipher file output path
 *
 * @return          [fage]		if successful o,otherwise failed
 */
int sm4_gcmf_encrypt_file(sm4_gcmf_context * ctx, char *infpath, char *outfpath);


/**
 * decrypt file
 *
 * @param  ctx      [in]		gcm file context
 *
 * @param  infpath  [in]		cipher file input path
 *
 * @param  outfpath [in]		plaintext file output path
 *
 * @return          [flag]		if successful o,otherwise failed
 */
int sm4_gcmf_decrypt_file(sm4_gcmf_context * ctx, char *infpath, char *outfpath);

#ifdef __cplusplus
}
#endif
#endif
#endif
