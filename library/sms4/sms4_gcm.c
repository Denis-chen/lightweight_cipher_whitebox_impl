/*
    created by Weijie Li
 */
#include <sms4/sms4.h>
#include <openssl/modes.h>

void sms4_gcm128_init(SMS4_GCM128_CONTEXT *ctx, SMS4_KEY *key) {
    CRYPTO_gcm128_init(ctx, key, (block128_f)sms4_encrypt);
}

void sms4_gcm128_setiv(SMS4_GCM128_CONTEXT *ctx, const unsigned char *ivec,
                      size_t len) {
    CRYPTO_gcm128_setiv(ctx, ivec, len);
}

int sms4_gcm128_aad(SMS4_GCM128_CONTEXT *ctx, const unsigned char *aad,
                    size_t len) {
    return CRYPTO_gcm128_aad(ctx, aad, len);
}

int sms4_gcm128_encrypt(const unsigned char *in, unsigned char *out,
                        size_t length, SMS4_GCM128_CONTEXT *ctx, const int enc) {
    int ret;
    if (enc)
        ret = CRYPTO_gcm128_encrypt(ctx, in, out, length);
    else
        ret = CRYPTO_gcm128_decrypt(ctx, in, out, length);
    return ret;
}

void sms4_gcm128_tag(SMS4_GCM128_CONTEXT *ctx, unsigned char *tag,
                    size_t len) {
    CRYPTO_gcm128_tag(ctx, tag, len);
}

int sms4_gcm128_finish(SMS4_GCM128_CONTEXT *ctx, const unsigned char *tag,
                      size_t len) {
    return CRYPTO_gcm128_finish(ctx, tag, len);
}

void sms4_gcm128_release(SMS4_GCM128_CONTEXT *ctx) {
    CRYPTO_gcm128_release(ctx);
}
