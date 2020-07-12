/*************************************
@berif  : this file is the implemention of Lblock in the way of 8-bit S-box. 
		  It is writed to test the speed of LBlock80.
@author : Shusheng Liu
@time   : 11.4 2011
@Copyright (c) 2011, south china normal university
**************************************/

#ifndef _LBLOCK80_8BIT_H_
#define _LBLOCK80_8BIT_H_

#include "LblockSbox_8Bit.h"

void LBlock80_encrypte_function(const uint8 *plain,const uint8 *key, const int rounds, uint8 *cipher)
{
	uint8 K[10];
	uint8 state[8];
	uint8 temp[4];
	uint8 t;

	uint8 *stateR = state;
	uint8 *stateL = state + 4;
	uint8 *tempPiont;
	uint8 i;

	state[0] = plain[0];  
	state[1] = plain[1]; 
	state[2] = plain[2]; 
	state[3] = plain[3];
	state[4] = plain[4]; 
	state[5] = plain[5];
	state[6] = plain[6]; 
	state[7] = plain[7];

	K[0] = key[0];  
	K[1] = key[1];  
	K[2] = key[2];
	K[3] = key[3];
	K[4] = key[4];
	K[5] = key[5];
	K[6] = key[6]; 
	K[7] = key[7]; 
	K[8] = key[8];
	K[9] = key[9];

	for(i=1; i<rounds; i++)
	{
		//add key and sbox

		temp[0] = sbox0[stateL[0] ^ K[6]];
		temp[1] = sbox1[stateL[1] ^ K[7]];
		temp[2] = sbox2[stateL[2] ^ K[8]];
		temp[3] = sbox3[stateL[3] ^ K[9]];

		//rotation
		t = stateR[0];

		stateR[0] = stateR[3] ^ ((temp[0]>>4) | (temp[1] & 0xF0));	  
		stateR[3] = stateR[2] ^ ((temp[3]<<4) | (temp[2] & 0x0F));
		stateR[2] = stateR[1] ^ ((temp[2]>>4) | (temp[3] & 0xF0));
		stateR[1] = t         ^ ((temp[1]<<4) | (temp[0] & 0x0F));

		tempPiont = stateL;
		stateL = stateR;
		stateR = tempPiont;

		//key schedule
		temp[3] = K[9]; temp[2] = K[8];
		temp[1] = K[7]; temp[0] = K[6];

		K[9] = (K[6]<<5)|(K[5]>>3); K[9]=sboxk[K[9]];
		K[8] = (K[5]<<5)|(K[4]>>3);
		K[7] = (K[4]<<5)|(K[3]>>3);
		K[6] = ((K[3]<<5)|(K[2]>>3)) ^ (i>>2);
		K[5] = ((K[2]<<5)|(K[1]>>3)) ^ (i<<6);
		K[4] = (K[1]<<5)|(K[0]>>3);

		K[3] = (K[0]<<5)|(temp[3]>>3);
		K[2] = (temp[3]<<5) | (temp[2]>>3);
		K[1] = (temp[2]<<5) | (temp[1]>>3);
		K[0] = (temp[1]<<5) | (temp[0]>>3);
	}
	//Last round
	//add key and sbox
	temp[0] = sbox0[stateL[0] ^ K[6]];
	temp[1] = sbox1[stateL[1] ^ K[7]];
	temp[2] = sbox2[stateL[2] ^ K[8]];
	temp[3] = sbox3[stateL[3] ^ K[9]];

	//rotation
	t=stateR[0];

	stateR[0] = stateR[3] ^ ((temp[0]>>4) | (temp[1] & 0xF0));	  
	stateR[3] = stateR[2] ^ ((temp[3]<<4) | (temp[2] & 0x0F));
	stateR[2] = stateR[1] ^ ((temp[2]>>4) | (temp[3] & 0xF0));
	stateR[1] = t         ^ ((temp[1]<<4) | (temp[0] & 0x0F));

	//output cipher
	cipher[0]=stateR[0]; cipher[1]=stateR[1];
	cipher[2]=stateR[2]; cipher[3]=stateR[3];

	cipher[4]=stateL[0]; cipher[5]=stateL[1];
	cipher[6]=stateL[2]; cipher[7]=stateL[3];
}

#endif