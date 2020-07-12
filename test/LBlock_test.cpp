/*
 * @Author: Weijie Li 
 * @Date: 2018-01-31 15:15:39 
 * @Last Modified by: Weijie Li
 * @Last Modified time: 2018-01-31 18:27:10
 */

#include <stdio.h>
#include <time.h>
#include <LBlock/LBlock80.h>
#include <LBlock/LBlock80wb.h>
const int Experiment_Times = 1;

#include "count_cycles.h"

int lblock_standard_test()
{
	uint8 plaintext[16];
	uint8 key[20];
	uint8 ciphertext[16];
	clock_t begin,end; 
	
	//test the implementation is right or not.
	int i;
	uint8 j=0;
	for(i=15; i>=0; i--) 
	{
		plaintext[i]=0;
		j++;
	}
	j=0;
	for(i=19; i>3; i--) 
	{
		key[i]=0;
		j++;
	}
	//key[3]=15; key[2]=14; key[1]=13; key[0]=12;
	key[3]=0; key[2]=0; key[1]=0; key[0]=0;
	
	begin = clock();
	set_cycles_start();

	for(i=0; i<Experiment_Times; i++)
	{
		LBlock80_encrypte_algorithm(plaintext,32,key,ciphertext);
	}

	set_cycles_ends();
	end = clock();
	double ts = end-begin;
	ts = ts * 1000 / CLOCKS_PER_SEC;
	printf("Encryption %d using %lf ms\n",Experiment_Times, ts);
	printf("Encryption %d using %llu cycles\n",Experiment_Times, get_cycles_elapsed());
	for(i=15; i>=0; i--)
	{
		printf("%x ",ciphertext[i]);		
	}
	printf("\n");

	return 0;
}

int lblock_wb_test()
{
	uint8 plaintext[16];
	uint8 key[20];
	uint8 ciphertext[16];
	clock_t program_start,program_end; 
	
	//test the implementation is right or not.
	int i;
	uint8 j=0;
	for(i=15; i>=0; i--) 
	{
		plaintext[i]=0;
		j++;
	}
	// j = 1;
	// for(i=15; i>=8; i--) 
	// {
	// 	plaintext[i]=j;
	// 	j++;
	// }
	j=0;
	for(i=19; i>3; i--) 
	{
		key[i]=0;
		j++;
	}
	//key[3]=15; key[2]=14; key[1]=13; key[0]=12;
	key[3]=0; key[2]=0; key[1]=0; key[0]=0;
	LBlock80wb_ctx ctx;

	program_start = clock();
	set_cycles_start();

	gen_LBlock80_wb_ctx(&ctx, key);

	set_cycles_ends();
	program_end = clock();

	double ts = program_end - program_start;
	ts = ts*1000/CLOCKS_PER_SEC;

	printf("The generating of encryptions table spend %lf ms\n", ts);
	printf("The generating of encryptions table spend %llu cycles\n", get_cycles_elapsed());

	program_start = clock();
	set_cycles_start();

	for(i=0; i<Experiment_Times; i++)
	{
		LBlock80wb_encrypte_algorithm(plaintext, &ctx, ciphertext);
	}
	
	set_cycles_ends();
	program_end = clock();

	ts = program_end - program_start;
	ts = ts*1000/CLOCKS_PER_SEC;

	printf("Encryption %d using %lf ms\n",Experiment_Times, ts);
	printf("Encryption %d using %llu cycles\n",Experiment_Times, get_cycles_elapsed());
	for(i=15; i>=0; i--)
	{
		printf("%x ",ciphertext[i]);		
	}
	printf("\n");

	return 0;
}

int lblock_test_main() {
	printf("LBlock80 standard test:\n");
	lblock_standard_test();
	printf("LBlock80 wb test:\n");
	lblock_wb_test();
	printf("\n");
	return 0;
}