#include <stdio.h>
#include <time.h>
#include <klein/klein64.h>
#include <klein/klein80.h>
#include <klein/klein96.h>

#include <klein/klein80wb.h>

#include "count_cycles.h"

#define EXPERIMENT_TIMES 1

int klein_stand()
{

    //test values
	const uint8_t key8[8] = {
		0
	    //0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF,
		//0x12, 0x34, 0x56, 0x78, 0x90, 0xAB, 0xCD, 0xEF,
	};


	const uint8_t key10[10] = {
		0
		//0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF,
		//0x12, 0x34, 0x56, 0x78, 0x90, 0xAB, 0xCD, 0xEF, 0x12, 0x34,
	};


	const uint8_t key12[12] = {
		//0
		0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF,
		//0x12, 0x34, 0x56, 0x78, 0x90, 0xAB, 0xCD, 0xEF, 0x12, 0x34, 0x56, 0x78,
	};

	uint8_t message[8] = {
		//0
		// 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF,
		0x12, 0x34, 0x56, 0x78, 0x90, 0xAB, 0xCD, 0xEF,
	};

	uint8_t cipher[8] = {0};

	int i;
	
	clock_t program_start, program_end = 0;
	program_start = clock();
	set_cycles_start();

	for(i = 0; i < EXPERIMENT_TIMES; i++)
	{
	    klein80_encrypt(message, key10, cipher);
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

	program_start = clock();
	set_cycles_start();
	

	for(i = 0; i < EXPERIMENT_TIMES; i++)
	{
	    klein80_decrypt(cipher, key10, message);
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

int klein_wb()
{

	//klein_testvectors();

    //test values
	const uint8_t key8[8] = {
		0
	    //0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF,
		//0x12, 0x34, 0x56, 0x78, 0x90, 0xAB, 0xCD, 0xEF,
	};


	const uint8_t key10[10] = {
		0
		//0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF,
		//0x12, 0x34, 0x56, 0x78, 0x90, 0xAB, 0xCD, 0xEF, 0x12, 0x34,
	};


	const uint8_t key12[12] = {
		//0
		0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF,
		//0x12, 0x34, 0x56, 0x78, 0x90, 0xAB, 0xCD, 0xEF, 0x12, 0x34, 0x56, 0x78,
	};

	uint8_t message[8] = {
		//0
		// 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF,
		0x12, 0x34, 0x56, 0x78, 0x90, 0xAB, 0xCD, 0xEF,
	};

	uint8_t cipher[8] = {0};

	klein80_wb_t wb_ctx;



	int i;
	clock_t program_start, program_end = 0;
	program_start = clock();
	set_cycles_start();

	gen_klein80wb_table(key10, &wb_ctx);
	
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
        //klein64_encrypt(message, key8, cipher);
	    klein80wb_encrypt(message, &wb_ctx, cipher);
	    //klein96_encrypt(message, key12, cipher);
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
	// 	//klein64_decrypt(cipher, key8, message);
	//     klein80_decrypt(cipher, key10, message);
	// 	//klein96_decrypt(cipher, key12, message);
	// }

	// program_end = clock();

	// printf("The %d times decryptions spend %ld ms\n", i, program_end - program_start);

	// for (i = 0; i < 8; i++) {
	// 	printf("%02X ", message[i]);
	// }

	// printf("\n");

	return 0;
}

int klein_test_main() {
	printf("klein standard test:\n");
	klein_stand();
	printf("\nklein wb test:\n");
	klein_wb();
	printf("\n");
	return 0;
}