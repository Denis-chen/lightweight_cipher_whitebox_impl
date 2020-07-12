


#include <affine/NTLUtils.h>

#include <NTL/mat_GF2.h>
#include <NTL/vec_GF2.h>

// #include <LBlock/LBlock80.h>

void print_u8(uint8_t *out, int len) {
    int ct = 8;
    for (int i=0; i<len; i++) {
        if (i%ct==0) printf("\n");
        printf("%02X", *(out+i));
    }
    printf("\n");
}

// int test_sms4_wb128(int idx);
int NTLTest();

int klein_test_main();

int lblock_test_main();
int present80_test_main();
int sm4_test_main();


int main(int argc, char **argv) 
{
    
    printf("-------------klein test-------------\n");
    klein_test_main();


    printf("-------------LBlock test-------------\n");
    lblock_test_main();

    printf("-------------present80 test-------------\n");
    present80_test_main();

    printf("-------------sm4 test-------------\n");
    sm4_test_main();

    return 0;
}
