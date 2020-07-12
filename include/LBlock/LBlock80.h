/*************************************
@berif  : this file is the implemention of Lblock in the way of 4-bit S-box. 
		  It is writed to test the speed of LBlock80.
@author : Shusheng Liu
@time   : 11.4 2011
@Copyright (c) 2011, south china normal university
**************************************/

#ifndef _LBLOCK_H_
#define _LBLOCK_H_


#ifdef __cplusplus 
extern "C" { 
#endif

#include "LblocSbox.h"

void LBlock80_encrypte_algorithm(uint8 *plain, const uint8 rounds, uint8 *key, uint8* cipher)
{
	//uint8 round_key[8];
	uint8 K[20];
	uint8 state[16];
	uint8 temp[8];
	uint8 *stateR = state;
	uint8 *stateL = state + 8;
	uint8 *tempPiont;
	uint8 i;

	stateR[0] = plain[0]; stateR[1] = plain[1];
	stateR[2] = plain[2]; stateR[3] = plain[3];
	stateR[4] = plain[4]; stateR[5] = plain[5];
	stateR[6] = plain[6]; stateR[7] = plain[7];

	stateL[0] = plain[8];  stateL[1] = plain[9];
	stateL[2] = plain[10]; stateL[3] = plain[11];
	stateL[4] = plain[12]; stateL[5] = plain[13];
	stateL[6] = plain[14]; stateL[7] = plain[15];

	K[0]=key[0];   K[1]=key[1];   K[2]=key[2];   K[3]=key[3];
	K[4]=key[4];   K[5]=key[5];   K[6]=key[6];   K[7]=key[7];
	K[8]=key[8];   K[9]=key[9];   K[10]=key[10]; K[11]=key[11];
	K[12]=key[12]; K[13]=key[13]; K[14]=key[14]; K[15]=key[15];
	K[16]=key[16]; K[17]=key[17]; K[18]=key[18]; K[19]=key[19];    

	// printf("Round 0:\t");
	// for(j=15; j>=0; j--) printf("%x ",state[j]);
	// printf("\n");

	for(i=1; i < rounds; i++)
	{
		// printf("Round %i key : ", i);
		// for(j=19; j>=12; j--)printf("%x ",K[j]);
		// printf("\n");
		
		//one round of LBlock
		temp[0] = stateR[0]; temp[1]=stateR[1];
		
		stateR[0] = stateR[6] ^ sbox1[K[13] ^ stateL[1]];
		stateR[1] = stateR[7] ^ sbox3[K[15] ^ stateL[3]];
		stateR[6] = stateR[4] ^ sbox4[K[16] ^ stateL[4]];
		stateR[7] = stateR[5] ^ sbox6[K[18] ^ stateL[6]];
		stateR[4] = stateR[2] ^ sbox5[K[17] ^ stateL[5]];
		stateR[5] = stateR[3] ^ sbox7[K[19] ^ stateL[7]];
		stateR[2] = temp[0]   ^ sbox0[K[12] ^ stateL[0]];
		stateR[3] = temp[1]   ^ sbox2[K[14] ^ stateL[2]];

		// printf("Round %d:\t", i);
		// for(j=7; j>=0; j--) printf("%x ",stateR[j]);
		// printf("\n");

		tempPiont = stateR;
		stateR = stateL;
		stateL = tempPiont;

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

		K[19] = sbox9[K[19]];
		K[18] = sbox8[K[18]];

		K[11] = (K[11] ^ (i<<2)) & 0x0F;
		K[12] =  K[12] ^  (i>>2);

		
	}

	//last round of LBlock
	temp[0] = stateR[0]; temp[1]=stateR[1];
	stateR[0] = stateR[6] ^ sbox1[K[13] ^ stateL[1]];
	stateR[1] = stateR[7] ^ sbox3[K[15] ^ stateL[3]];
	stateR[6] = stateR[4] ^ sbox4[K[16] ^ stateL[4]];
	stateR[7] = stateR[5] ^ sbox6[K[18] ^ stateL[6]];
	stateR[4] = stateR[2] ^ sbox5[K[17] ^ stateL[5]];
	stateR[5] = stateR[3] ^ sbox7[K[19] ^ stateL[7]];
	stateR[2] = temp[0]   ^ sbox0[K[12] ^ stateL[0]];
	stateR[3] = temp[1]   ^ sbox2[K[14] ^ stateL[2]];

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


#ifdef __cplusplus 
}
#endif

#endif