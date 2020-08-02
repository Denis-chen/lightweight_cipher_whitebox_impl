/*
 * @Author: Weijie Li
 * @Date: 2018-01-24 00:26:02
 * @Last Modified by: Weijie Li
 * @Last Modified time: 2018-01-29 22:54:37
 */
#ifdef _cplusplus
extern "C" {
#endif
#include "random.h"
#include "WBMatrix.h"
#include "structure.h"
#include "inverse.h"


#ifdef _cplusplus
}
#endif
#include <present/present80wb.h>
#include <internal/ayssl_random.h>
 // actually is (sbox[] << 4)
static const uint8_t sbox[16] = {
	0xC0, 0x50, 0x60, 0xB0, 0x90, 0x00, 0xA0, 0xD0, 0x30, 0xE0, 0xF0, 0x80, 0x40, 0x70, 0x10, 0x20,
};

typedef struct _wb_helper {
	M8 f[PRESENT_ROUNDS][8];//matrix_transform_t
	M8 f_inv[PRESENT_ROUNDS][8];

	// combine every 4 f matrix
	M32 fc_inv[PRESENT_ROUNDS][2];

	M8 g[PRESENT_ROUNDS + 1][8];
	M8 g_inv[PRESENT_ROUNDS + 1][8];
} wb_helper;

void present_wb_helper_init(int rounds, wb_helper &wbh) {
	int i, k;

	for (i = 0; i < rounds; i++) {
		for (k = 0; k < 8; k++) {
			genMatpairM8(&wbh.f[i][k], &wbh.f_inv[i][k]);//Éú³É¿ÉÄæ¾ØÕó¼°ÆäÄæ¾ØÕó
			//identityM8(&wbh.f[i][k]);
			//identityM8(&wbh.f_inv[i][k]);
		}

		//Ã¿ËÄ¸ö¾ØÕó×é³É¶Ô½ÇÏß¾ØÕó
		MatrixcomM8to32(wbh.f_inv[i][0], wbh.f_inv[i][1], wbh.f_inv[i][2], wbh.f_inv[i][3], &wbh.fc_inv[i][0]);
		MatrixcomM8to32(wbh.f_inv[i][4], wbh.f_inv[i][5], wbh.f_inv[i][6], wbh.f_inv[i][7], &wbh.fc_inv[i][1]);

	}
	for (i = 0; i <= rounds; i++) {
		for (k = 0; k < 8; k++) {
			genMatpairM8(&wbh.g[i][k], &wbh.g_inv[i][k]);//¿ÉÄæ¾ØÕóÒÔ¼°Äæ¾ØÕóÉú³É
			//identityM8(&wbh.g[i][k]);
			//identityM8(&wbh.g_inv[i][k]);
		}
	}
}

void Mat8MulMatM32(M8 a, M32 b, M32 *c) {

	// divide matrix b to 4 8x8 matrixs
	M8 b1, b2, b3, b4;
	for (int i = 0; i < 8; i++) {
		b1.M[i] = (b.M[i] >> 24) & 0xff;
		b2.M[i] = (b.M[i] >> 16) & 0xff;
		b3.M[i] = (b.M[i] >> 8) & 0xff;
		b4.M[i] = (b.M[i]) & 0xff;
	}
	//a * bi respectively;
	MatMulMatM8(a, b1, &b1);
	MatMulMatM8(a, b2, &b2);
	MatMulMatM8(a, b3, &b3);
	MatMulMatM8(a, b4, &b4);
	//combine bi to c
	for (int i = 0; i < 8; i++) {
		c->M[i] = (b1.M[i] << 24) ^ (b2.M[i] << 16) ^ (b3.M[i] << 8) ^ (b4.M[i]);
	}
	//return 0;
}//8*8¾ØÕó³Ë8*32¾ØÕó

void initMatrixFromBit(M32 *Mat, long* data) {
	long i, j;
	M32* mat = Mat;
	for (i = 0; i <8;i++) {

		mat->M[i] = data[i];
	}
}

uint8_t applyMatToU8(const M8 mat, uint8_t data) {

	return MatMulNumM8(mat, data);
}
uint32_t applyMatToU32(const M32 mat, uint32_t data) {

	return MatMulNumM32(mat, data);
}
uint8_t present_sbox8(uint8_t x) {
	uint8_t y = 0;
	y = sbox[x >> 4 & 0x0F] | (sbox[x & 0x0F] >> 4);
	return y;
}//SºÐ±ä»»£¿ÓÒÒÆËÄÎ»

/**
 * @brief combine addRoundKey and sbox, then confuse it
 *
 * @param key
 * @param wbh
 * @param ctx
 */
void present_wb_init_rkAsbox(const uint8_t *key, const wb_helper &wbh, present_wb_ctx &ctx) {
	const uint8_t rounds = ctx.rounds;
	uint8_t round_counter = 1;//ÂÖ´ÎÊý³õÊ¼»¯
	int i, j;

	uint8_t state[3];
	uint8_t round_key[10];
	// combine addRoundKey and sbox



	for (i = 0; i < 8; i++) {
		for (j = 0; j < 256; j++) {
			uint8_t n_int;
			n_int = applyMatToU8(wbh.g_inv[0][i], j);
			n_int = present_sbox8(n_int ^ key[i]);//ÂÖÃÜÔ¿Òì»ò²¢½øÐÐSºÐ´ú»»
			ctx.rk[0][i][j] = applyMatToU8(wbh.f[0][i], n_int);
			n_int = applyMatToU8(wbh.g[0][i], j);
			ctx.stmp[i][j] = n_int;
		}
	}

	// update keyÃÜÔ¿À©Õ¹£¨Ñ­»·ÓÒÒÆ18Î»£©
	round_key[9] = key[6] << 5 | key[7] >> 3;
	round_key[8] = key[5] << 5 | key[6] >> 3;
	round_key[7] = key[4] << 5 | key[5] >> 3;
	round_key[6] = key[3] << 5 | key[4] >> 3;
	round_key[5] = key[2] << 5 | key[3] >> 3;
	round_key[4] = key[1] << 5 | key[2] >> 3;
	round_key[3] = key[0] << 5 | key[1] >> 3;
	round_key[2] = key[9] << 5 | key[0] >> 3;
	round_key[1] = key[8] << 5 | key[9] >> 3;
	round_key[0] = key[7] << 5 | key[8] >> 3;

	round_key[0] = (round_key[0] & 0x0F) | sbox[round_key[0] >> 4];

	round_key[7] ^= round_counter >> 1;
	round_key[8] ^= round_counter << 7;

	for (round_counter = 2; round_counter <= rounds; round_counter++) {
		for (i = 0; i < 8; i++) {
			for (j = 0; j < 256; j++) {
				uint8_t n_int;
				n_int = applyMatToU8(wbh.g_inv[round_counter - 1][i], j);
				n_int = present_sbox8(n_int ^ round_key[i]);
				ctx.rk[round_counter - 1][i][j] = applyMatToU8(wbh.f[round_counter - 1][i], n_int);
			}
		}
		round_key[5] ^= round_counter << 2; // do this first, which may be faster

		// use state[] for temporary storage
		state[2] = round_key[9];
		state[1] = round_key[8];
		state[0] = round_key[7];

		round_key[9] = round_key[6] << 5 | round_key[7] >> 3;
		round_key[8] = round_key[5] << 5 | round_key[6] >> 3;
		round_key[7] = round_key[4] << 5 | round_key[5] >> 3;
		round_key[6] = round_key[3] << 5 | round_key[4] >> 3;
		round_key[5] = round_key[2] << 5 | round_key[3] >> 3;
		round_key[4] = round_key[1] << 5 | round_key[2] >> 3;
		round_key[3] = round_key[0] << 5 | round_key[1] >> 3;
		round_key[2] = state[2] << 5 | round_key[0] >> 3;
		round_key[1] = state[1] << 5 | state[2] >> 3;
		round_key[0] = state[0] << 5 | state[1] >> 3;

		round_key[0] = (round_key[0] & 0x0F) | sbox[round_key[0] >> 4];
	}

	for (i = 0; i < 8; i++) {
		for (j = 0; j < 256; j++) {
			uint8_t n_int;
			n_int = applyMatToU8(wbh.g_inv[rounds][i], j);
			n_int = (n_int ^ round_key[i]);
			ctx.rk[rounds][i][j] = n_int;
		}
	}

}

void present_wb_init_player(const wb_helper &wbh, present_wb_ctx &ctx) {
	// permutation£¬PÖÃ»»
	uint8_t round_counter;

	static long pLayer0_hex[] = { 0x80000000,
									0x08000000,
									0x00800000,
									0x00080000,
									0x00008000,
									0x00000800,
									0x00000080,
									0x00000008 };
	static long pLayer1_hex[] = { 0x40000000,
									0x04000000,
									0x00400000,
									0x00040000,
									0x00004000,
									0x00000400,
									0x00000040,
									0x00000004 };
	static long pLayer2_hex[] = { 0x20000000,
									0x02000000,
									0x00200000,
									0x00020000,
									0x00002000,
									0x00000200,
									0x00000020,
									0x00000002 };
	static long pLayer3_hex[] = { 0x10000000,
									0x01000000,
									0x00100000,
									0x00010000,
									0x00001000,
									0x00000100,
									0x00000010,
									0x00000001 };

	static M32 p0, p1, p2, p3;
	initMatrixFromBit(&p0, pLayer0_hex);
	initMatrixFromBit(&p1, pLayer1_hex);
	initMatrixFromBit(&p2, pLayer2_hex);
	initMatrixFromBit(&p3, pLayer3_hex);

	for (round_counter = 0; round_counter < ctx.rounds; round_counter++) {

		M32 Mat;
		initM32(&Mat);
		MatMulMatM32(p0, wbh.fc_inv[round_counter][0], &Mat);//µÃMat8*32
		Mat8MulMatM32(wbh.g[round_counter + 1][0], Mat, &ctx.pLayer[round_counter][0]);
		MatMulMatM32(p0, wbh.fc_inv[round_counter][1], &Mat);
		Mat8MulMatM32(wbh.g[round_counter + 1][1], Mat, &ctx.pLayer[round_counter][1]);
		MatMulMatM32(p1, wbh.fc_inv[round_counter][0], &Mat);
		Mat8MulMatM32(wbh.g[round_counter + 1][2], Mat, &ctx.pLayer[round_counter][2]);
		MatMulMatM32(p1, wbh.fc_inv[round_counter][1], &Mat);
		Mat8MulMatM32(wbh.g[round_counter + 1][3], Mat, &ctx.pLayer[round_counter][3]);
		MatMulMatM32(p2, wbh.fc_inv[round_counter][0], &Mat);
		Mat8MulMatM32(wbh.g[round_counter + 1][4], Mat, &ctx.pLayer[round_counter][4]);
		MatMulMatM32(p2, wbh.fc_inv[round_counter][1], &Mat);
		Mat8MulMatM32(wbh.g[round_counter + 1][5], Mat, &ctx.pLayer[round_counter][5]);
		MatMulMatM32(p3, wbh.fc_inv[round_counter][0], &Mat);
		Mat8MulMatM32(wbh.g[round_counter + 1][6], Mat, &ctx.pLayer[round_counter][6]);   
		MatMulMatM32(p3, wbh.fc_inv[round_counter][1], &Mat);
		Mat8MulMatM32(wbh.g[round_counter + 1][7], Mat, &ctx.pLayer[round_counter][7]);

	}
}

#if PRESENT_WB_DEBUG

void 	present_wb_debug_prepare(const wb_helper &wbh, present_wb_ctx &ctx) {
	int i, k;
	int rounds = ctx.rounds;
	for (i = 0; i < rounds; i++) {
		ctx.fc_inv[i][0] = wbh.fc_inv[i][0];
		ctx.fc_inv[i][1] = wbh.fc_inv[i][1];
		for (k = 0; k < 8; k++) {
			ctx.f_inv[i][k] = wbh.f_inv[i][k];
		}
	}
	for (i = 0; i <= rounds; i++) {
		for (k = 0; k < 8; k++) {
			ctx.g_inv[i][k] = wbh.g_inv[i][k];
		}
	}
}

#endif //PRESENT_WB_DEBUG

void present_wb_init(const uint8_t *key, present_wb_ctx &ctx) {
	ctx.rounds = PRESENT_ROUNDS;
	wb_helper wbh;
	present_wb_helper_init(ctx.rounds, wbh);
	present_wb_init_rkAsbox(key, wbh, ctx);
	present_wb_init_player(wbh, ctx);
#if PRESENT_WB_DEBUG
	present_wb_debug_prepare(wbh, ctx);
#endif //PRESENT_WB_DEBUG

}

void present_wb_release(present_wb_ctx ctx) {}

// full-round should be 31, i.e. rounds = 31
// plain and cipher can overlap, so do key and cipher
void present_wb_enc(const uint8_t *plain, const present_wb_ctx &ctx, uint8_t *cipher)
{
	uint8_t round_counter = 1;

	uint8_t state[8];

	int rounds = ctx.rounds;

	int i;
	for (i = 0; i < 8; i++) {
		cipher[i] = ctx.stmp[i][plain[i]];
	}

	for (round_counter = 0; round_counter < rounds; round_counter++) {

#if PRESENT_WB_DEBUG
		printf("Round %d:\n", round_counter + 1);
		printf("State : ");
		for (i = 0; i < 8; i++) {
			printf("%02X", applyMatToU8(ctx.g_inv[round_counter][i], cipher[i]));
		}
		printf("\n");
#endif //PRESENT_WB_DEBUG

		for (i = 0; i < 8; i++) {
			state[i] = ctx.rk[round_counter][i][cipher[i]];
		}

#if PRESENT_WB_DEBUG
		printf("AddRk and sbox : ");
		for (i = 0; i < 8; i++) {
			printf("%02X", applyMatToU8(ctx.f_inv[round_counter][i], state[i]));
		}
		printf("\n");
#endif //PRESENT_WB_DEBUG

		uint32_t x0, x1;
		x0 = state[0] << 24 |
			state[1] << 16 |
			state[2] << 8 |
			state[3];
		x1 = state[4] << 24 |
			state[5] << 16 |
			state[6] << 8 |
			state[7];
		for (i = 0; i < 8; i += 2) {
			cipher[i + 0] = (uint8_t)((applyMatToU32(ctx.pLayer[round_counter][i + 0], x0))>>24)&0xff;
			cipher[i + 1] = (uint8_t)((applyMatToU32(ctx.pLayer[round_counter][i + 1], x1))>>24)&0xff;
		}

#if PRESENT_WB_DEBUG
		printf("Cipher : ");
		for (i = 0; i < 8; i++) {
			printf("%02X", applyMatToU8(ctx.g_inv[round_counter + 1][i], cipher[i]));
		}
		printf("\n");
#endif //PRESENT_WB_DEBUG

	}

	// for ( i=0; i<8; i++) {
	// 	printf("%02X", cipher[i]);
	// }
	// printf("\n");

	for (i = 0; i < 8; i++) {
		cipher[i] = ctx.rk[round_counter][i][cipher[i]];
	}
}
