/*
 * @Author: Weijie Li 
 * @Date: 2017-11-24 10:05:18 
 * @Last Modified by: Weijie Li
 * @Last Modified time: 2017-11-24 16:28:07
 */



#ifndef _DEBUG_UTILS_H_
#define _DEBUG_UTILS_H_

#include <stdio.h>


static char hexchars[] ="0123456789ABCDEF";

#define DUMP(d, len) \
        dumpPointer(d, len, __FILE__, __LINE__)

#define DUMPPINFO(d, len, info) \
        dumpPointerInfo(d, len, __FILE__, __LINE__, info)

static void dumpPointer(const uint8_t *d, int len, char _file[], int _line) {
    printf("dump <%s,%d>: %d\n", _file, _line, len);
    size_t i;
    for (i = 0; i < len; i++) {
        printf("%02X", (*d++));
    }
    printf("\n");
}

static void dumpPointerInfo(const uint8_t *d, int len, char _file[], int _line, char info[]) {
    printf("[%s] dump <%s,%d>: %d\n", info, _file, _line, len);
    size_t i;
    for (i = 0; i < len; i++) {
        printf("%02X", (*d++));
    }
    printf("\n");
}


/**
 * hexify(in, out, len):
 * Convert ${len} bytes from ${in} into hexadecimal, writing the resulting
 * 2 * ${len} bytes to ${out}; and append a NUL byte.
 */
static void hexify(const uint8_t *in, char *out, size_t len) {
    char * p = out;
    size_t i;
    
    for (i = 0; i < len; i++) {
        *p++ = hexchars[in[i] >> 4];
        *p++ = hexchars[in[i] & 0x0f];
    }
    *p = '\0';
}

/**
 * unhexify(in, out, len):
 * Convert 2 * ${len} hexadecimal characters from ${in} to ${len} bytes
 * and write them to ${out}.  This function will only fail if the input is
 * not a sequence of hexadecimal characters.
 */
static int unhexify(const char *in, uint8_t *out, size_t len) {
    size_t i;
    
    /* Make sure we have at least 2 * ${len} hex characters. */
    for (i = 0; i < 2 * len; i++) {
        if ((in[i] == '\0') || (strchr(hexchars, in[i]) == NULL))
            goto err0;
    }
    
    for (i = 0; i < len; i++) {
        out[i] = (strchr(hexchars, in[2 * i]) - hexchars) & 0x0f;
        out[i] <<= 4;
        out[i] += (strchr(hexchars, in[2 * i + 1]) - hexchars) & 0x0f;
    }
    
    /* Success! */
    return (0);
    
err0:
    /* Bad input string. */
    return (-1);
}



#endif /* _DEBUG_UTILS_H_ */