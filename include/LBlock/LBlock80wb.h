/*
 * @Author: Weijie Li 
 * @Date: 2018-01-25 16:12:11 
 * @Last Modified by: Weijie Li
 * @Last Modified time: 2018-01-26 11:01:42
 */


#ifndef _LBLOCK80WB_H_
#define _LBLOCK80WB_H_

#include "LblocSbox.h"
#include <cstdint>

/*********
 * 
 * LBlock config 
 * 
 * ***************/

#define LBLOCK80_ROUNDS 32
#define LBLOCK_WB_DEBUG 0

typedef struct LBlock80wb_ctx {
	int rounds ;
	uint8_t Te[LBLOCK80_ROUNDS][8][16];
	uint8_t remap[LBLOCK80_ROUNDS][8][16];
	uint8_t fi[2][8][16];
	uint8_t g[2][8][16];
	
	#if LBLOCK_WB_DEBUG
		uint8_t restore[LBLOCK80_ROUNDS+2][8][16];
	#endif //LBLOCK_WB_DEBUG
	
} LBlock80wb_ctx;

int gen_LBlock80_wb_ctx(LBlock80wb_ctx *ctx, const uint8_t *key);

void LBlock80wb_encrypte_algorithm(const uint8_t *plain, const LBlock80wb_ctx *ctx, uint8_t* cipher);

#endif //_LBLOCK80WB_H_