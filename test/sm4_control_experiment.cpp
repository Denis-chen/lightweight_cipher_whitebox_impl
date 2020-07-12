/*
 * @Author: Weijie Li 
 * @Date: 2018-01-29 23:02:23 
 * @Last Modified by: Weijie Li
 * @Last Modified time: 2018-02-01 00:08:11
 */

#include <stdio.h>
#include <time.h>

#include <sms4/sms4.h>
#include <sms4wb/sms4_wb_generator.h>
#include <sms4wb/sms4_wb.h>

#include "count_cycles.h"

#define EXPERIMENT_TIMES 1


const uint8_t key16[16] = {
    0
    //0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF,
    //0x12, 0x34, 0x56, 0x78, 0x90, 0xAB, 0xCD, 0xEF, 0x12, 0x34, 0x56, 0x78, 0x90, 0xAB, 0xCD, 0xEF,
};

uint8_t message[16] = {
    //0
    // 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF,
    0x12, 0x34, 0x56, 0x78, 0x90, 0xAB, 0xCD, 0xEF, 0x12, 0x34, 0x56, 0x78, 0x90, 0xAB, 0xCD, 0xEF,
};

uint8_t cipher[16] = {0};

int sm4_stand()
{
	

	int i;

    sms4_key_t sm4_key;
    sms4_set_encrypt_key(&sm4_key, key16);
	
	clock_t program_start, program_end = 0;
	program_start = clock();
	set_cycles_start();

	for(i = 0; i < EXPERIMENT_TIMES; i++)
	{
	    sms4_encrypt(message, cipher, &sm4_key);
    }

	set_cycles_ends();
	program_end = clock();

	double ts = program_end-program_start;
	ts = ts*1000/CLOCKS_PER_SEC;

	printf("The %d times encryptions spend %lf ms\n", i, ts);
	printf("The %d times encryptions spend %llu cycles\n", i, get_cycles_elapsed());

	for (i = 0; i < 8; i++) {
	    printf("%02X ", cipher[i]);
	}

	printf("\n");

    sms4_set_decrypt_key(&sm4_key, key16);

	program_start = clock();
	set_cycles_start();

	for(i = 0; i < EXPERIMENT_TIMES; i++)
	{
	    sms4_decrypt(cipher, message, &sm4_key);
	}

	set_cycles_ends();
	program_end = clock();
	ts = program_end-program_start;
	ts = ts*1000/CLOCKS_PER_SEC;

	printf("The %d times decryptions spend %lf ms\n", i, ts);
	printf("The %d times decryptions spend %llu cycles\n", i, get_cycles_elapsed());

	for (i = 0; i < 8; i++) {
		printf("%02X ", message[i]);
	}

	printf("\n");

	return 0;
}

int sm4_wb()
{

	sms4_wb_t wb_ctx;

	int i;
	clock_t program_start, program_end = 0;
	program_start = clock();
	set_cycles_start();

	sms4_wb_gen_tables(key16, &wb_ctx, SMS4_ENCRYPT);

	set_cycles_ends();
	program_end = clock();
	double ts = program_end - program_start;
	ts = ts*1000/CLOCKS_PER_SEC;

	printf("The generating of encryptions table spend %lf ms\n", ts);
	printf("The generating of encryptions table spend %llu cycles\n", get_cycles_elapsed());
	
	program_start = clock();
	set_cycles_start();

	for(i = 0; i < EXPERIMENT_TIMES; i++)
	{
        //sm464_encrypt(message, key8, cipher);
	    sms4_wb_enc(message, cipher, &wb_ctx);
	    //sm496_encrypt(message, key12, cipher);
    }

	set_cycles_ends();
	program_end = clock();

	ts = program_end - program_start;
	ts = ts*1000/CLOCKS_PER_SEC;

	printf("The %d times encryptions spend %lf ms\n", i, ts);
	printf("The %d times encryptions spend %llu cycles\n", i, get_cycles_elapsed());

	for (i = 0; i < 8; i++) {
	    printf("%02X ", cipher[i]);
	}

	printf("\n");

	// program_start = clock();

	// for(i = 0; i < EXPERIMENT_TIMES; i++)
	// {
	// 	//sm464_decrypt(cipher, key8, message);
	//     sm480_decrypt(cipher, key10, message);
	// 	//sm496_decrypt(cipher, key12, message);
	// }

	// program_end = clock();

	// printf("The %d times decryptions spend %ld ms\n", i, program_end - program_start);

	// for (i = 0; i < 8; i++) {
	// 	printf("%02X ", message[i]);
	// }

	// printf("\n");

	return 0;
}

int sm4_test_main() {
	printf("sm4 standard test:\n");
	sm4_stand();
	printf("sm4 wb test:\n");
	sm4_wb();
	printf("\n");
	return 0;
}