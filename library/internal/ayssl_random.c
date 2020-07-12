#include <internal/ayssl_random.h>

#include <time.h>
#include <stdlib.h>
#include <stdbool.h>


static bool isInintial = false;

void initial() {
   if (isInintial)
       return;
   isInintial = true;
   srand(time(NULL));   // should only be called once
}

int ayssl_random_number(int upper) {
    #ifndef DUBUG
    initial();
    #endif
    return rand()%upper; 
}

int ayssl_random_shutter_u32(u32* list, int len) {
    #ifndef DUBUG
    initial();
    #endif
    int t;
    while(len) {
        int r = ayssl_random_number(len);
        t = *(list+len);
        *(list+len) = *(list+r);
        *(list+r) = t;
        len--;
    }
    return 0;
}

int ayssl_random_shutter_u8(u8* list, int len) {
#ifndef DUBUG
    initial();
#endif
    int t;
    while(len) {
        u32 r = ayssl_random_number(len--);
        t = *(list+len);
        *(list+len) = *(list+r);
        *(list+r) = t;
    }
    return 0;
}
