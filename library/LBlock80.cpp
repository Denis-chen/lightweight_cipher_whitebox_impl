/*
 * @Author: Weijie Li 
 * @Date: 2018-01-25 16:12:15 
 * @Last Modified by: Weijie Li
 * @Last Modified time: 2018-01-26 17:40:57
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

#include <LBlock/LBlock80wb.h>

#include <cstdio>

#include<NTL/mat_GF2.h>
#include<affine/MatrixUtils.h>
using namespace std;
typedef struct LBlock80wb_assistant {
	uint8_t f[LBLOCK80_ROUNDS+3][8][16];
	uint8_t fi[LBLOCK80_ROUNDS+3][8][16];
}LBlock80wb_assistant;

int embedded_sbox_key(LBlock80wb_ctx *ctx, const uint8_t *key) {
	uint8_t K[20];
	uint8_t temp[8];
	//每个4bits
	K[0]=key[0];   K[1]=key[1];   K[2]=key[2];   K[3]=key[3];
	K[4]=key[4];   K[5]=key[5];   K[6]=key[6];   K[7]=key[7];
	K[8]=key[8];   K[9]=key[9];   K[10]=key[10]; K[11]=key[11];
	K[12]=key[12]; K[13]=key[13]; K[14]=key[14]; K[15]=key[15];
	K[16]=key[16]; K[17]=key[17]; K[18]=key[18]; K[19]=key[19]; 
	//20个key
	int i,d;
	int rounds = LBLOCK80_ROUNDS;
	ctx->rounds = rounds; //32rounds
	for (i=1; i<=rounds; i++) {
		// test OK
		// printf("Round %i key : ", i);
		// for(j=19; j>=12; j--)printf("%x ",K[j]);  一轮里有8个密钥
		// printf("\n");
		for (d=0; d<16; d++) {
			ctx->Te[i-1][0][d] = sbox0[ d^K[12] ]; //这里他妈是3维数组
			ctx->Te[i-1][1][d] = sbox1[ d^K[13] ];//SubKey的生成？
			ctx->Te[i-1][2][d] = sbox2[ d^K[14] ];
			ctx->Te[i-1][3][d] = sbox3[ d^K[15] ];
			ctx->Te[i-1][4][d] = sbox4[ d^K[16] ];
			ctx->Te[i-1][5][d] = sbox5[ d^K[17] ];
			ctx->Te[i-1][6][d] = sbox6[ d^K[18] ];
			ctx->Te[i-1][7][d] = sbox7[ d^K[19] ];
		}

		if (i!=rounds) { //i取1-31才能进来，80bits的密钥再次被打乱，然后取最左边32bits
			//Key schedule
			temp[7] = K[19]; temp[6] = K[18]; temp[5] = K[17]; 
			temp[4] = K[16]; temp[3] = K[15]; temp[2] = K[14];
			temp[1] = K[13]; temp[0] = K[12]; 

			K[19] = ((K[12]<<1)|(K[11]>>3))     & 0x0F;  K[18] = ((K[11]<<1)|(K[10]>>3))     & 0x0F;
			K[17] = ((K[10]<<1)|(K[9]>>3))      & 0x0F;  K[16] = ((K[9]<<1)|(K[8]>>3))       & 0x0F;
			K[15] = ((K[8]<<1) |(K[7]>>3))      & 0x0F;  K[14] = ((K[7]<<1)|(K[6]>>3))       & 0x0F;
			K[13] = ((K[6]<<1) |(K[5]>>3))      & 0x0F;  K[12] = ((K[5]<<1)|(K[4]>>3))       & 0x0F;
			K[11] = ((K[4]<<1) |(K[3]>>3))      & 0x0F;  K[10] = ((K[3]<<1)|(K[2]>>3))       & 0x0F;
			K[9]  = ((K[2]<<1) |(K[1]>>3))      & 0x0F;  K[8]  = ((K[1]<<1)|(K[0]>>3))       & 0x0F;
			K[7]  = ((K[0]<<1) |(temp[7]>>3))   & 0x0F;  K[6]  = ((temp[7]<<1)|(temp[6]>>3)) & 0x0F;
			K[5]  = ((temp[6]<<1) |(temp[5]>>3))& 0x0F;  K[4]  = ((temp[5]<<1)|(temp[4]>>3)) & 0x0F;
			K[3]  = ((temp[4]<<1) |(temp[3]>>3))& 0x0F;  K[2]  = ((temp[3]<<1)|(temp[2]>>3)) & 0x0F;
			K[1]  = ((temp[2]<<1) |(temp[1]>>3))& 0x0F;  K[0]  = ((temp[1]<<1)|(temp[0]>>3)) & 0x0F;
			//所有的K值发生替换，有点像轮函数F
			K[19] = sbox9[K[19]];
			K[18] = sbox8[K[18]];

			K[11] = (K[11] ^ (i<<2)) & 0x0F;
			K[12] =  K[12] ^  (i>>2);
		}
	}
	return 0;
}

void initVecFromBit(V4 *vec, long data) {
	vec->V = 0;
	for (int i = 0; i < 4; i++) {
		
		vec->V = vec->V ^ ((data % 2) << i);
		data = data >> 1;
	}
}//初始化vec


int gen_bijections_table(LBlock80wb_ctx *ctx, LBlock80wb_assistant *assistant) {
	int i,j,k;
	int rounds = ctx->rounds;
	//NTL::mat_GF2 X,Y;
	M4 Mat4[3];
	V4 ans;
	for (i=0; i<rounds+3; i++) {
		for (k=0; k<8; k++) {
			//genRandomInvMatrix(X, Y, 4);// 4x4的矩阵
			//initinvbaseM8(initM4_max);
			genMatpairM4(&Mat4[0], &Mat4[1]);//生成一对可逆矩阵

			for (j=0; j<16; j++) {
				//NTL::vec_GF2 tv;
				V4 vec4;
				//initVecFromBit(tv, j, 4);
				initVecFromBit(&vec4, j);

				MatMulVecM4(Mat4[0],vec4,&ans);
				//tv = X*tv; //向量乘矩阵              
				int t = (uint32_t) ans.V;
				//int t = getDigitalFromVec(tv);//将向量抽数据出来
				assistant->f[i][k][j] = t;
				assistant->fi[i][k][t] = j;//可逆的encoding之F函数？
			}
		}
	}
	for (k=0; k<8; k++) {
		for (j=0; j<16; j++) {
			ctx->fi[0][k][j] = assistant->fi[0][k][j];
			ctx->fi[1][k][j] = assistant->fi[1][k][j];
			ctx->g[1][k][j] = assistant->f[rounds+1][k][j];
			ctx->g[0][k][j] = assistant->f[rounds][k][j];
		}
	}	

	#if LBLOCK_WB_DEBUG
		for (i=0; i<rounds+2; i++) {
			for (k=0; k<8; k++) {
				for (j=0; j<16; j++) {
					ctx->restore[i][k][j] = assistant->f[i][k][j];
				}
			}
		}
	#endif //LBLOCK_WB_DEBUG

	return 0;
}

int apply_bijections_table(LBlock80wb_ctx *ctx, LBlock80wb_assistant *assistant) {
	//看名字应该是和应用混淆有关
	int i,d;
	int rounds = ctx->rounds;//32轮
	for (i=0; i<rounds; i++) {	
		uint8_t Te_T[8][16];
		for (d=0; d<16; d++) { 
			Te_T[0][d] = 	assistant->fi[i+2][2][ \
									ctx->Te[i][0][assistant->f[i+1][0][d]]] ;
			Te_T[1][d] = 	assistant->fi[i+2][0][ \
									ctx->Te[i][1][assistant->f[i+1][1][d]]] ;
			Te_T[2][d] = 	assistant->fi[i+2][3][ \
									ctx->Te[i][2][assistant->f[i+1][2][d]]] ;
			Te_T[3][d] = 	assistant->fi[i+2][1][ \
									ctx->Te[i][3][assistant->f[i+1][3][d]]] ;
			Te_T[4][d] = 	assistant->fi[i+2][6][ \
									ctx->Te[i][4][assistant->f[i+1][4][d]]] ;
			Te_T[5][d] = 	assistant->fi[i+2][4][ \
									ctx->Te[i][5][assistant->f[i+1][5][d]]] ;
			Te_T[6][d] = 	assistant->fi[i+2][7][ \
									ctx->Te[i][6][assistant->f[i+1][6][d]]] ;
			Te_T[7][d] = 	assistant->fi[i+2][5][ \
									ctx->Te[i][7][assistant->f[i+1][7][d]]] ;
								
		}
		for (d=0; d<16; d++) {
			ctx->Te[i][0][d] = Te_T[0][d];
			ctx->Te[i][1][d] = Te_T[1][d];
			ctx->Te[i][2][d] = Te_T[2][d];
			ctx->Te[i][3][d] = Te_T[3][d];
			ctx->Te[i][4][d] = Te_T[4][d];
			ctx->Te[i][5][d] = Te_T[5][d];
			ctx->Te[i][6][d] = Te_T[6][d];
			ctx->Te[i][7][d] = Te_T[7][d];

			ctx->remap[i][0][d] = 	assistant->fi[i+2][0][ \
									assistant->f [i   ][6][d]];
			ctx->remap[i][1][d] = 	assistant->fi[i+2][1][ \
									assistant->f [i   ][7][d]];
			ctx->remap[i][2][d] = 	assistant->fi[i+2][2][ \
									assistant->f [i   ][0][d]];
			ctx->remap[i][3][d] = 	assistant->fi[i+2][3][ \
									assistant->f [i   ][1][d]];
			ctx->remap[i][4][d] = 	assistant->fi[i+2][4][ \
									assistant->f [i   ][2][d]];
			ctx->remap[i][5][d] = 	assistant->fi[i+2][5][ \
									assistant->f [i   ][3][d]];
			ctx->remap[i][6][d] = 	assistant->fi[i+2][6][ \
									assistant->f [i   ][4][d]];
			ctx->remap[i][7][d] = 	assistant->fi[i+2][7][ \
									assistant->f [i   ][5][d]];
		}
		//一看这个remap就知道是重新映射的意思 。。。
	}
	return 0;
}

int gen_LBlock80_wb_ctx(LBlock80wb_ctx *ctx, const uint8_t *key) {
	//有个generate，生成的意思，应该是生成ctx这个对象的函数
	ctx->rounds = LBLOCK80_ROUNDS;
	embedded_sbox_key(ctx, key);//调用了嵌入S盒密钥

	LBlock80wb_assistant assistant;
	gen_bijections_table(ctx, &assistant);
	apply_bijections_table(ctx, &assistant);//制造了一次混淆，对ctx进行混淆

	// correctnesss check out
	// int temptrya =3;
	// int temptryb =8;
	// printf("%X ^ %X : %X ^ %X\n%X : %X : %X\n", temptrya, temptryb, assistant.f[1][2][temptrya], assistant.f[1][2][temptryb], \
	// 		temptrya ^ temptryb, assistant.f[1][2][temptrya] ^assistant.f[1][2][temptryb], assistant.fi[1][2][assistant.f[1][2][temptrya] ^assistant.f[1][2][temptryb]]);
	
	return 0;
}

// 64 bits plain
void LBlock80wb_encrypte_algorithm(const uint8_t *plain, const LBlock80wb_ctx *ctx, uint8_t* cipher){
	//using 4 bits each uint8 果然这个B也是和我采用同样方法
	//重头戏加密算法来了
	
	uint8_t state[16];
	uint8_t temp[8];
	uint8_t *stateR = state;
	uint8_t *stateL = state + 8;
	uint8_t *tempPiont;
	uint8_t i;

	stateR[0] = plain[0]; stateR[1] = plain[1];
	stateR[2] = plain[2]; stateR[3] = plain[3];
	stateR[4] = plain[4]; stateR[5] = plain[5];
	stateR[6] = plain[6]; stateR[7] = plain[7];

	stateL[0] = plain[8];  stateL[1] = plain[9];
	stateL[2] = plain[10]; stateL[3] = plain[11];
	stateL[4] = plain[12]; stateL[5] = plain[13];
	stateL[6] = plain[14]; stateL[7] = plain[15];
	//明文也是4bits的，神奇

	for (i=0; i<8; i++)  {
		stateL[i] = ctx->fi[1][i][stateL[i]];
		stateR[i] = ctx->fi[0][i][stateR[i]];
	}//开始加密了么

       
#if LBLOCK_WB_DEBUG
	printf("Round 0:\t");
	for(j=7; j>=0; j--) printf("%x ", ctx->restore[1][j][stateL[j]]);
	for(j=7; j>=0; j--) printf("%x ", ctx->restore[0][j][stateR[j]]);
	printf("\n");
#endif //LBLOCK_WB_DEBUG
	

	for(i=0; i < ctx->rounds; i++)
	{
		//one round of LBlock
		temp[0] = stateR[0]; temp[1]=stateR[1];

		#if LBLOCK_WB_DEBUG
			printf("Round %d:\t\n", i);

			printf("Left:\t");
			printf("%X%X%X%X%X%X%X%X\n",  \
			 														ctx->restore[i+2][7][ctx->Te[i][6][stateL[6]]], \
																	ctx->restore[i+2][6][ctx->Te[i][4][stateL[4]]], \
																	ctx->restore[i+2][5][ctx->Te[i][7][stateL[7]]], \
																	ctx->restore[i+2][4][ctx->Te[i][5][stateL[5]]], \
																	ctx->restore[i+2][3][ctx->Te[i][2][stateL[2]]], \
																	ctx->restore[i+2][2][ctx->Te[i][0][stateL[0]]], \
																	ctx->restore[i+2][1][ctx->Te[i][3][stateL[3]]], \
																	ctx->restore[i+2][0][ctx->Te[i][1][stateL[1]]] );
			printf("Right:\t");
			printf("%X%X%X%X%X%X%X%X\n",  \
																	ctx->restore[i+2][7][ctx->remap[i][7][stateR[5]]], \
																	ctx->restore[i+2][6][ctx->remap[i][6][stateR[4]]], \
																	ctx->restore[i+2][5][ctx->remap[i][5][stateR[3]]], \
																	ctx->restore[i+2][4][ctx->remap[i][4][stateR[2]]], \
																	ctx->restore[i+2][3][ctx->remap[i][3][temp[1]  ]],  \
																	ctx->restore[i+2][2][ctx->remap[i][2][temp[0]  ]], \
																	ctx->restore[i+2][1][ctx->remap[i][1][stateR[7]]], \
																	ctx->restore[i+2][0][ctx->remap[i][0][stateR[6]]] );
		#endif //LBLOCK_WB_DEBUG

		stateR[0] = ctx->remap[i][0][stateR[6]] ^ ctx->Te[i][1][stateL[1]];
		stateR[1] = ctx->remap[i][1][stateR[7]] ^ ctx->Te[i][3][stateL[3]];
		stateR[6] = ctx->remap[i][6][stateR[4]] ^ ctx->Te[i][4][stateL[4]];
		stateR[7] = ctx->remap[i][7][stateR[5]] ^ ctx->Te[i][6][stateL[6]];
		stateR[4] = ctx->remap[i][4][stateR[2]] ^ ctx->Te[i][5][stateL[5]];
		stateR[5] = ctx->remap[i][5][stateR[3]] ^ ctx->Te[i][7][stateL[7]];
		stateR[2] = ctx->remap[i][2][temp[0]  ] ^ ctx->Te[i][0][stateL[0]];
		stateR[3] = ctx->remap[i][3][temp[1]  ] ^ ctx->Te[i][2][stateL[2]];

		tempPiont = stateR;
		stateR = stateL;
		stateL = tempPiont;

		        
#if LBLOCK_WB_DEBUG
		for(j=7; j>=0; j--) printf("%x ", ctx->restore[i+2][j][stateL[j]]);
        for(j=7; j>=0; j--) printf("%x ", ctx->restore[i+1][j][stateR[j]]);
        printf("\n");
#endif //LBLOCK_WB_DEBUG

	}
	tempPiont = stateR;
	stateR = stateL;
	stateL = tempPiont;

	for (i=0; i<8; i++) {
		stateR[i] = ctx->g[1][i][stateR[i]];  
		stateL[i] = ctx->g[0][i][stateL[i]];
	}

	//ouput cipher text
	cipher[0] = stateR[0]; cipher[1] = stateR[1]; 
	cipher[2] = stateR[2]; cipher[3] = stateR[3]; 
	cipher[4] = stateR[4]; cipher[5] = stateR[5]; 
	cipher[6] = stateR[6]; cipher[7] = stateR[7];

	cipher[8]  = stateL[0]; cipher[9]  = stateL[1];
	cipher[10] = stateL[2]; cipher[11] = stateL[3];
	cipher[12] = stateL[4]; cipher[13] = stateL[5];
	cipher[14] = stateL[6]; cipher[15] = stateL[7];
}

