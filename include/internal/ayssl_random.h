#ifndef HEADER_AYSSL_RANDOM_H
#define HEADER_AYSSL_RANDOM_H

#ifdef __cplusplus 
extern "C" { 
#endif

typedef unsigned int u32;
typedef unsigned short u16;
typedef unsigned char u8;

int ayssl_random_number(int upper);
int ayssl_random_shutter_u32(u32* list, int len);
int ayssl_random_shutter_u8(u8* list, int len);

#ifdef __cplusplus 
}
#endif


#endif
