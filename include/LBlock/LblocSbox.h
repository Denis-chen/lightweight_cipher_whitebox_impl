/*************************************
@berif  : this file is the 4-bit S-box of Lblock. 
		  It is writed to test the speed of LBlock80.
@author : Shusheng Liu
@time   : 11.4 2011
@Copyright (c) 2010, south china normal university
**************************************/

#ifndef _LBLOCKSBOX_H_
#define _LBLOCKSBOX_H_

typedef unsigned char uint8;
typedef unsigned short uint16;
typedef unsigned int uint32;

static const uint8 sbox0[16] = {14,9,15,0,13,4,10,11,1,2,8,3,7,6,12,5};

static const uint8 sbox1[16] = {4,11,14,9,15,13,0,10,7,12,5,6,2,8,1,3};

static const uint8 sbox2[16] = {1,14,7,12,15,13,0,6,11,5,9,3,2,4,8,10};

static const uint8 sbox3[16] = {7,6,8,11,0,15,3,14,9,10,12,13,5,2,4,1};

static const uint8 sbox4[16] = {14,5,15,0,7,2,12,13,1,8,4,9,11,10,6,3};

static const uint8 sbox5[16] = {2,13,11,12,15,14,0,9,7,10,6,3,1,8,4,5};

static const uint8 sbox6[16] = {11,9,4,14,0,15,10,13,6,12,5,7,3,8,1,2};

static const uint8 sbox7[16] = {13,10,15,0,14,4,9,11,2,1,8,3,7,5,12,6};

static const uint8 sbox8[16] = {8,7,14,5,15,13,0,6,11,12,9,10,2,4,1,3};

static const uint8 sbox9[16] = {11,5,15,0,7,2,9,13,4,8,1,12,14,10,3,6};

#endif