/*
 * @Author: Weijie Li 
 * @Date: 2017-11-07 19:24:54 
 * @Last Modified by: Weijie Li
 * @Last Modified time: 2017-11-24 11:27:29
 */



 #include <stdint.h>
 #include <sms4/sms4.h>

 #include <sms4wb/sms4_wb.h>


void sms4_wb_gen_rk(uint8_t *key);

void sms4_wb_gen_coding_map();

void sms4_wb_gen_mixing_bijections();

void sms4_wb_gen_ext_encoding();

int sms4_wb_gen_tables(const uint8_t *key, sms4_wb_t *sms4_wb_ctx, int enc);

