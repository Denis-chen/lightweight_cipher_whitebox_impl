/*
* sm4gcmf.c
*
*  Created on: 2017-10-18
*      Author: lzj
*      Updated: Payne
*
*/

#include <sms4/sms4.h>

int sm4_gcmf_init(sm4_gcmf_context *ctx) {
	memset(ctx, 0, sizeof(sm4_gcmf_context));
	ctx->gcm = (SMS4_GCM128_CONTEXT *)malloc(sizeof(SMS4_GCM128_CONTEXT));
	return 0;
}

int sm4_gcmf_free(sm4_gcmf_context *ctx) {
	int ret = sms4_gcm128_finish(ctx->gcm, ctx->tag, SMS4_GCM_FILE_TAG_LEN);
	if (ret != 0) return ret;
	free(ctx->gcm);
	memset(ctx, 0, sizeof(sm4_gcmf_context));
	return 0;
}

int sm4_gcmf_set_iv(sm4_gcmf_context *ctx, const unsigned char * iv, size_t len) {
	sms4_gcm128_setiv(ctx->gcm, iv, len);
	return 0;
}

int sm4_gcmf_set_key(sm4_gcmf_context *ctx, const unsigned char * key,
		size_t len) {
	SMS4_KEY sms4_key;
	sms4_set_encrypt_key(&sms4_key, key);
	sms4_gcm128_init(ctx->gcm, &sms4_key);
	return 0;
}

int sm4_gcmf_encrypt_file(sm4_gcmf_context * ctx, char *infpath, char *outfpath) {
	int len_outfpath, file_size, block_size;
	FILE *infp, *tmpfp, *outfp;

	// #0 Open file
	len_outfpath = strlen(outfpath);
	char *tmp_file_path = (char *)malloc((len_outfpath + 4) * sizeof(char));
	memcpy(tmp_file_path, outfpath, len_outfpath);
	memcpy(tmp_file_path + len_outfpath, ".tmp", 4);

	if ((infp = fopen(infpath, "rb")) && infp != NULL) return -1;
	if ((tmpfp = fopen(tmp_file_path, "wb+")) && tmpfp != NULL) return -2;
	if ((outfp = fopen(outfpath, "wb+")) && outfp != NULL) return -3;

	// #1 Encrypt file (slice blocks)

	unsigned char buf[SM4_GCM_FILE_MAX_BLOCK_LEN];
	unsigned char out[SM4_GCM_FILE_MAX_BLOCK_LEN];
	unsigned char size_buf[sizeof(size_t)];
	memset(buf, 0, SM4_GCM_FILE_MAX_BLOCK_LEN);
	memset(out, 0, SM4_GCM_FILE_MAX_BLOCK_LEN);

	// Get size of file
	fseek(infp, 0, SEEK_END);
	file_size = ftell(infp);
	memset(size_buf, 0, sizeof(size_t));
	memcpy(size_buf, &file_size, sizeof(size_t));
	fseek(infp, 0, SEEK_SET);

	// GCM Block encrypt
	while (
		(block_size = fread(buf, sizeof(unsigned char), SM4_GCM_FILE_MAX_BLOCK_LEN, infp)) &&
		block_size != 0
	) {
		sms4_gcm128_encrypt(buf, out, block_size, ctx->gcm, 1);
		fwrite(out, sizeof(unsigned char), block_size, tmpfp);

		memset(buf, 0, SM4_GCM_FILE_MAX_BLOCK_LEN);
		memset(out, 0, SM4_GCM_FILE_MAX_BLOCK_LEN);
	}

	// Get Tag of GCM encrypt
	sms4_gcm128_tag(ctx->gcm, ctx->tag, SMS4_GCM_FILE_TAG_LEN);

	// Finish encrypt
	fflush(tmpfp);
	fseek(tmpfp, 0, SEEK_SET);
	fclose(infp);

	// #2 Write cipher file

	// add file flag
	fputs(SMS4_GCM_FILE_MAGIC_TAG, outfp);
	// add version (default: 1)
	int ver = SMS4_GCM_FILE_VERSION;
	fwrite(&ver, sizeof(int), 1, outfp);
	// add tag
	fwrite(&ctx->tag, sizeof(unsigned char), SMS4_GCM_FILE_TAG_LEN, outfp);
	// add file len
	fwrite(size_buf, sizeof(unsigned char), sizeof(size_t), outfp);
	// copy cipher
	memset(buf, 0, SM4_GCM_FILE_MAX_BLOCK_LEN);
	while (
		(block_size = fread(buf, sizeof(unsigned char), SM4_GCM_FILE_MAX_BLOCK_LEN, tmpfp)) &&
		block_size != 0
	) {
		fwrite(buf, sizeof(unsigned char), block_size, outfp);
		memset(buf, 0, SM4_GCM_FILE_MAX_BLOCK_LEN);
	}

	fclose(tmpfp);
	fclose(outfp);
	remove(tmp_file_path);
	free(tmp_file_path);

	return 0;
}

int sm4_gcmf_decrypt_file(sm4_gcmf_context * ctx, char *infpath, char *outfpath) {
	int ret = 0;
	int file_size = 0, read_size = 0, block_size;
	unsigned char size_buf[sizeof(size_t)];
	FILE *infp, *outfp;
	if ((infp = fopen(infpath, "rb")) && infp != NULL) {
		return -1;
	}

	char flag[SMS4_GCM_FILE_TAG_LEN + 1];
	memset(flag, 0, sizeof(flag));
	//read file flag
	fgets(flag, SMS4_GCM_FILE_MAGIC_TAG_LEN + 1, infp);
	if (strcmp(flag, SMS4_GCM_FILE_MAGIC_TAG)) {
		ret = -2;
		goto end;
	}

	//read version
	int version;
	if (fread(&version, sizeof(version), 1, infp) == 0) {
		ret = -3;
		goto end;
	}
	if (version != SMS4_GCM_FILE_VERSION) {
		ret = -3;
		goto end;
	}

	//read tag
	if (fread(ctx->tag, sizeof(unsigned char), SMS4_GCM_FILE_TAG_LEN, infp) == 0) {
		ret = -4;
		goto end;
	}

	memset(size_buf, 0, sizeof(size_t));
	// read real length
	if (fread(size_buf, sizeof(unsigned char), sizeof(size_t), infp) == 0) {
		ret = -5;
		goto end;
	}
	file_size = (int) ((((size_buf[3] & 0xff) << 24)
			| ((size_buf[2] & 0xff) << 16) | ((size_buf[1] & 0xff) << 8)
			| ((size_buf[0] & 0xff) << 0)));

	// Write Plaintext
	if ((outfp = fopen(outfpath, "wb+")) && outfp != NULL) {
		ret = -6;
		goto end;
	}

	unsigned char buf[SM4_GCM_FILE_MAX_BLOCK_LEN];
	unsigned char out[SM4_GCM_FILE_MAX_BLOCK_LEN];
	memset(buf, 0, SM4_GCM_FILE_MAX_BLOCK_LEN);
	memset(out, 0, SM4_GCM_FILE_MAX_BLOCK_LEN);
	sms4_gcm128_tag(ctx->gcm, ctx->tag, SMS4_GCM_FILE_TAG_LEN);

	while (
		(block_size = fread(buf, sizeof(unsigned char), SM4_GCM_FILE_MAX_BLOCK_LEN, infp)) &&
		block_size != 0
	) {
		sms4_gcm128_encrypt(buf, out, block_size, ctx->gcm, 0);
		read_size += block_size;
		fwrite(out, sizeof(unsigned char), block_size, outfp);

		memset(buf, 0, SM4_GCM_FILE_MAX_BLOCK_LEN);
		memset(out, 0, SM4_GCM_FILE_MAX_BLOCK_LEN);
	}

	fflush(outfp);
	// Invalid file
	if (read_size != file_size) {
		ret = -7;
		goto end;
	}

	end:
	fclose(infp);
	fclose(outfp);

	return ret;
}
