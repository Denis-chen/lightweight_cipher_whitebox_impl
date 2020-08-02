/*
 * @Author: Weijie Li
 * @Date: 2017-11-17 01:27:14
 * @Last Modified by: Weijie Li
 * @Last Modified time: 2017-11-19 12:27:56
 */

#include <klein/klein80wb.h>

//#include <affine/MatrixUtils.h>
//#include <affine/NTLUtils.h>
//#include <NTL/mat_GF2.h>

#include <wbmatrix/WBMatrix.h>


#if (defined(_WIN32) || defined(_WIN64)) && !defined(__MINGW32__)
typedef unsigned __int64 u64;
#define U64(C)  C##UI64
#elif defined(__arch64__)
typedef unsigned long u64;
#define U64(C)  C##UL
#else
typedef unsigned long long u64;
#define U64(C)  C##ULL
#endif

#define u64 uint64_t
#define u32 uint32_t
#define u16 uint16_t
#define u8  uint8_t


/*-
Te [x] = [x].[02, 01, 01, 03, 02, 01, 01, 03];
Te0[x] = [x].[02, 01, 01, 03];
Te1[x] = [x].[03, 02, 01, 01];
Te2[x] = [x].[01, 03, 02, 01];
Te3[x] = [x].[01, 01, 03, 02];
*/
#define Te0 (u32)((u64*)((u8*)Te+0))
#define Te1 (u32)((u64*)((u8*)Te+3))
#define Te2 (u32)((u64*)((u8*)Te+2))
#define Te3 (u32)((u64*)((u8*)Te+1))

static const u64 Te[256] = {
    U64(0x0000000000000000), U64(0x0301010203010102),
    U64(0x0602020406020204), U64(0x0503030605030306),
    U64(0x0C0404080C040408), U64(0x0F05050A0F05050A),
    U64(0x0A06060C0A06060C), U64(0x0907070E0907070E),
    U64(0x1808081018080810), U64(0x1B0909121B090912),
    U64(0x1E0A0A141E0A0A14), U64(0x1D0B0B161D0B0B16),
    U64(0x140C0C18140C0C18), U64(0x170D0D1A170D0D1A),
    U64(0x120E0E1C120E0E1C), U64(0x110F0F1E110F0F1E),
    U64(0x3010102030101020), U64(0x3311112233111122),
    U64(0x3612122436121224), U64(0x3513132635131326),
    U64(0x3C1414283C141428), U64(0x3F15152A3F15152A),
    U64(0x3A16162C3A16162C), U64(0x3917172E3917172E),
    U64(0x2818183028181830), U64(0x2B1919322B191932),
    U64(0x2E1A1A342E1A1A34), U64(0x2D1B1B362D1B1B36),
    U64(0x241C1C38241C1C38), U64(0x271D1D3A271D1D3A),
    U64(0x221E1E3C221E1E3C), U64(0x211F1F3E211F1F3E),
    U64(0x6020204060202040), U64(0x6321214263212142),
    U64(0x6622224466222244), U64(0x6523234665232346),
    U64(0x6C2424486C242448), U64(0x6F25254A6F25254A),
    U64(0x6A26264C6A26264C), U64(0x6927274E6927274E),
    U64(0x7828285078282850), U64(0x7B2929527B292952),
    U64(0x7E2A2A547E2A2A54), U64(0x7D2B2B567D2B2B56),
    U64(0x742C2C58742C2C58), U64(0x772D2D5A772D2D5A),
    U64(0x722E2E5C722E2E5C), U64(0x712F2F5E712F2F5E),
    U64(0x5030306050303060), U64(0x5331316253313162),
    U64(0x5632326456323264), U64(0x5533336655333366),
    U64(0x5C3434685C343468), U64(0x5F35356A5F35356A),
    U64(0x5A36366C5A36366C), U64(0x5937376E5937376E),
    U64(0x4838387048383870), U64(0x4B3939724B393972),
    U64(0x4E3A3A744E3A3A74), U64(0x4D3B3B764D3B3B76),
    U64(0x443C3C78443C3C78), U64(0x473D3D7A473D3D7A),
    U64(0x423E3E7C423E3E7C), U64(0x413F3F7E413F3F7E),
    U64(0xC0404080C0404080), U64(0xC3414182C3414182),
    U64(0xC6424284C6424284), U64(0xC5434386C5434386),
    U64(0xCC444488CC444488), U64(0xCF45458ACF45458A),
    U64(0xCA46468CCA46468C), U64(0xC947478EC947478E),
    U64(0xD8484890D8484890), U64(0xDB494992DB494992),
    U64(0xDE4A4A94DE4A4A94), U64(0xDD4B4B96DD4B4B96),
    U64(0xD44C4C98D44C4C98), U64(0xD74D4D9AD74D4D9A),
    U64(0xD24E4E9CD24E4E9C), U64(0xD14F4F9ED14F4F9E),
    U64(0xF05050A0F05050A0), U64(0xF35151A2F35151A2),
    U64(0xF65252A4F65252A4), U64(0xF55353A6F55353A6),
    U64(0xFC5454A8FC5454A8), U64(0xFF5555AAFF5555AA),
    U64(0xFA5656ACFA5656AC), U64(0xF95757AEF95757AE),
    U64(0xE85858B0E85858B0), U64(0xEB5959B2EB5959B2),
    U64(0xEE5A5AB4EE5A5AB4), U64(0xED5B5BB6ED5B5BB6),
    U64(0xE45C5CB8E45C5CB8), U64(0xE75D5DBAE75D5DBA),
    U64(0xE25E5EBCE25E5EBC), U64(0xE15F5FBEE15F5FBE),
    U64(0xA06060C0A06060C0), U64(0xA36161C2A36161C2),
    U64(0xA66262C4A66262C4), U64(0xA56363C6A56363C6),
    U64(0xAC6464C8AC6464C8), U64(0xAF6565CAAF6565CA),
    U64(0xAA6666CCAA6666CC), U64(0xA96767CEA96767CE),
    U64(0xB86868D0B86868D0), U64(0xBB6969D2BB6969D2),
    U64(0xBE6A6AD4BE6A6AD4), U64(0xBD6B6BD6BD6B6BD6),
    U64(0xB46C6CD8B46C6CD8), U64(0xB76D6DDAB76D6DDA),
    U64(0xB26E6EDCB26E6EDC), U64(0xB16F6FDEB16F6FDE),
    U64(0x907070E0907070E0), U64(0x937171E2937171E2),
    U64(0x967272E4967272E4), U64(0x957373E6957373E6),
    U64(0x9C7474E89C7474E8), U64(0x9F7575EA9F7575EA),
    U64(0x9A7676EC9A7676EC), U64(0x997777EE997777EE),
    U64(0x887878F0887878F0), U64(0x8B7979F28B7979F2),
    U64(0x8E7A7AF48E7A7AF4), U64(0x8D7B7BF68D7B7BF6),
    U64(0x847C7CF8847C7CF8), U64(0x877D7DFA877D7DFA),
    U64(0x827E7EFC827E7EFC), U64(0x817F7FFE817F7FFE),
    U64(0x9B80801B9B80801B), U64(0x9881811998818119),
    U64(0x9D82821F9D82821F), U64(0x9E83831D9E83831D),
    U64(0x9784841397848413), U64(0x9485851194858511),
    U64(0x9186861791868617), U64(0x9287871592878715),
    U64(0x8388880B8388880B), U64(0x8089890980898909),
    U64(0x858A8A0F858A8A0F), U64(0x868B8B0D868B8B0D),
    U64(0x8F8C8C038F8C8C03), U64(0x8C8D8D018C8D8D01),
    U64(0x898E8E07898E8E07), U64(0x8A8F8F058A8F8F05),
    U64(0xAB90903BAB90903B), U64(0xA8919139A8919139),
    U64(0xAD92923FAD92923F), U64(0xAE93933DAE93933D),
    U64(0xA7949433A7949433), U64(0xA4959531A4959531),
    U64(0xA1969637A1969637), U64(0xA2979735A2979735),
    U64(0xB398982BB398982B), U64(0xB0999929B0999929),
    U64(0xB59A9A2FB59A9A2F), U64(0xB69B9B2DB69B9B2D),
    U64(0xBF9C9C23BF9C9C23), U64(0xBC9D9D21BC9D9D21),
    U64(0xB99E9E27B99E9E27), U64(0xBA9F9F25BA9F9F25),
    U64(0xFBA0A05BFBA0A05B), U64(0xF8A1A159F8A1A159),
    U64(0xFDA2A25FFDA2A25F), U64(0xFEA3A35DFEA3A35D),
    U64(0xF7A4A453F7A4A453), U64(0xF4A5A551F4A5A551),
    U64(0xF1A6A657F1A6A657), U64(0xF2A7A755F2A7A755),
    U64(0xE3A8A84BE3A8A84B), U64(0xE0A9A949E0A9A949),
    U64(0xE5AAAA4FE5AAAA4F), U64(0xE6ABAB4DE6ABAB4D),
    U64(0xEFACAC43EFACAC43), U64(0xECADAD41ECADAD41),
    U64(0xE9AEAE47E9AEAE47), U64(0xEAAFAF45EAAFAF45),
    U64(0xCBB0B07BCBB0B07B), U64(0xC8B1B179C8B1B179),
    U64(0xCDB2B27FCDB2B27F), U64(0xCEB3B37DCEB3B37D),
    U64(0xC7B4B473C7B4B473), U64(0xC4B5B571C4B5B571),
    U64(0xC1B6B677C1B6B677), U64(0xC2B7B775C2B7B775),
    U64(0xD3B8B86BD3B8B86B), U64(0xD0B9B969D0B9B969),
    U64(0xD5BABA6FD5BABA6F), U64(0xD6BBBB6DD6BBBB6D),
    U64(0xDFBCBC63DFBCBC63), U64(0xDCBDBD61DCBDBD61),
    U64(0xD9BEBE67D9BEBE67), U64(0xDABFBF65DABFBF65),
    U64(0x5BC0C09B5BC0C09B), U64(0x58C1C19958C1C199),
    U64(0x5DC2C29F5DC2C29F), U64(0x5EC3C39D5EC3C39D),
    U64(0x57C4C49357C4C493), U64(0x54C5C59154C5C591),
    U64(0x51C6C69751C6C697), U64(0x52C7C79552C7C795),
    U64(0x43C8C88B43C8C88B), U64(0x40C9C98940C9C989),
    U64(0x45CACA8F45CACA8F), U64(0x46CBCB8D46CBCB8D),
    U64(0x4FCCCC834FCCCC83), U64(0x4CCDCD814CCDCD81),
    U64(0x49CECE8749CECE87), U64(0x4ACFCF854ACFCF85),
    U64(0x6BD0D0BB6BD0D0BB), U64(0x68D1D1B968D1D1B9),
    U64(0x6DD2D2BF6DD2D2BF), U64(0x6ED3D3BD6ED3D3BD),
    U64(0x67D4D4B367D4D4B3), U64(0x64D5D5B164D5D5B1),
    U64(0x61D6D6B761D6D6B7), U64(0x62D7D7B562D7D7B5),
    U64(0x73D8D8AB73D8D8AB), U64(0x70D9D9A970D9D9A9),
    U64(0x75DADAAF75DADAAF), U64(0x76DBDBAD76DBDBAD),
    U64(0x7FDCDCA37FDCDCA3), U64(0x7CDDDDA17CDDDDA1),
    U64(0x79DEDEA779DEDEA7), U64(0x7ADFDFA57ADFDFA5),
    U64(0x3BE0E0DB3BE0E0DB), U64(0x38E1E1D938E1E1D9),
    U64(0x3DE2E2DF3DE2E2DF), U64(0x3EE3E3DD3EE3E3DD),
    U64(0x37E4E4D337E4E4D3), U64(0x34E5E5D134E5E5D1),
    U64(0x31E6E6D731E6E6D7), U64(0x32E7E7D532E7E7D5),
    U64(0x23E8E8CB23E8E8CB), U64(0x20E9E9C920E9E9C9),
    U64(0x25EAEACF25EAEACF), U64(0x26EBEBCD26EBEBCD),
    U64(0x2FECECC32FECECC3), U64(0x2CEDEDC12CEDEDC1),
    U64(0x29EEEEC729EEEEC7), U64(0x2AEFEFC52AEFEFC5),
    U64(0x0BF0F0FB0BF0F0FB), U64(0x08F1F1F908F1F1F9),
    U64(0x0DF2F2FF0DF2F2FF), U64(0x0EF3F3FD0EF3F3FD),
    U64(0x07F4F4F307F4F4F3), U64(0x04F5F5F104F5F5F1),
    U64(0x01F6F6F701F6F6F7), U64(0x02F7F7F502F7F7F5),
    U64(0x13F8F8EB13F8F8EB), U64(0x10F9F9E910F9F9E9),
    U64(0x15FAFAEF15FAFAEF), U64(0x16FBFBED16FBFBED),
    U64(0x1FFCFCE31FFCFCE3), U64(0x1CFDFDE11CFDFDE1),
    U64(0x19FEFEE719FEFEE7), U64(0x1AFFFFE51AFFFFE5)
};

typedef struct klein_wb_assistant_t {
    int rounds;

#ifdef MAT_BIJECTIONS
    uint8_t fi[8][256][2], g[8][256][2];
#endif 

#ifdef AFFINE_IMPL
    Aff8 P[ROUNDS_80 + 1][8][2];
    Aff32 PC[ROUNDS_80 + 1][2][2];
    Aff32 xorMap[ROUNDS_80][2][2];

#endif
} klein_wb_assistant_t;


u32 inverse_u32(u32 a) {
    u32 b;
    u8* p = (u8*)&b;
    int i;
    for (i = 0; i < 4; i++) {     
        *p = (u8)(a >> (24 - 8 * i));
        p++;
    }
    return b;
}

/*void initVecFromBit8(V8 *vec, long data){   //init Vector with data
    vec->V = 0;
    for (int i = 0; i < 8; i++) {
        vec->V = vec->V ^ ((data % 2) << i);
        data = data >> 1;
    }
}*/

uint8_t applyMatToU8(const M8 mat, uint8_t data) {  //apply data in matrix
    //V8 vec8, ans;
    //uint8_t n_int;
    //initVecFromBit8(&vec8, data);
    //MatMulVecM8(mat, vec8, &ans);
    //n_int = (uint8_t)ans.V;
    return MatMulNumM8(mat, data);
}


int gen_klein80wb_ctx_table(const uint8_t* key, klein80_wb_t* klein80wb_ctx, int rounds) {
    int i, j, k;

    uint8_t round_key[10];
    uint8_t temp_state[8];

    klein80wb_ctx->rounds = rounds;

    round_key[0] = key[0];
    round_key[1] = key[1];
    round_key[2] = key[2];
    round_key[3] = key[3];
    round_key[4] = key[4];
    round_key[5] = key[5];
    round_key[6] = key[6];
    round_key[7] = key[7];
    round_key[8] = key[8];
    round_key[9] = key[9];


    uint8_t*** T = new uint8_t * *[rounds];   //Dynamically generate arrays T[rounds][8][256]
    for (i = 0; i < rounds; i++)
        T[i] = new uint8_t * [8];
    for (i = 0; i < rounds; i++)
        for (j = 0; j < 8; j++)
            T[i][j] = new uint8_t[256];
    //uint8_t T[rounds][8][256];


    for (i = 1; i <= rounds; i++) {
        // printf("round %d:\n", i);
        // for (j=0; j<8; j++) 
        //     printf("%02X", round_key[j]);
        // printf("\n");

        for (j = 0; j < 8; j++) {
            // printf("\tj:%d\tk:%02X\n", j, round_key[j]);
            for (k = 0; k < 256; k++) {
                T[i - 1][j][k] = sbox8[k ^ round_key[j]];

                // klein80wb_ctx->T[i-1][j][k] = k^ round_key[j];
                // printf("%02X:%02X\t", k, klein80wb_ctx->T[i-1][j][k]);
            }
            // printf("\n");
        }

        for (k = 0; k < 256; k++) {
            klein80wb_ctx->Te[i - 1][0][k] = inverse_u32(Te0[T[i - 1][2][k]]);
            klein80wb_ctx->Te[i - 1][1][k] = inverse_u32(Te1[T[i - 1][3][k]]);
            klein80wb_ctx->Te[i - 1][2][k] = inverse_u32(Te2[T[i - 1][4][k]]);
            klein80wb_ctx->Te[i - 1][3][k] = inverse_u32(Te3[T[i - 1][5][k]]);

            klein80wb_ctx->Te[i - 1][4][k] = inverse_u32(Te0[T[i - 1][6][k]]);
            klein80wb_ctx->Te[i - 1][5][k] = inverse_u32(Te1[T[i - 1][7][k]]);
            klein80wb_ctx->Te[i - 1][6][k] = inverse_u32(Te2[T[i - 1][0][k]]);
            klein80wb_ctx->Te[i - 1][7][k] = inverse_u32(Te3[T[i - 1][1][k]]);
        }

        temp_state[0] = round_key[0];
        temp_state[1] = round_key[1];
        temp_state[2] = round_key[2];
        temp_state[3] = round_key[3];
        temp_state[4] = round_key[4];
        temp_state[5] = round_key[5];

        round_key[0] = round_key[6];
        round_key[1] = round_key[7];
        round_key[2] = round_key[8] ^ i;
        round_key[3] = round_key[9];
        round_key[4] = round_key[5];

        round_key[5] = temp_state[1] ^ round_key[6];
        round_key[6] = sbox8[temp_state[2] ^ round_key[7]];
        round_key[7] = sbox8[temp_state[3] ^ round_key[8]];
        round_key[8] = temp_state[4] ^ round_key[9];
        round_key[9] = temp_state[0] ^ temp_state[5];
    }

    for (j = 0; j < 8; j++) {
        for (k = 0; k < 256; k++) {
            klein80wb_ctx->SE[0][j][k] = k ^ round_key[j];
        }
    }

    return 0;
}
    

#ifdef AFFINE_IMPL


int gen_affine_matrix(klein80_wb_t* ctx, klein_wb_assistant_t* assistant) {
    int i, j;
    int rounds = assistant->rounds;
    for (i = 0; i <= rounds; i++) {
        for (j = 0; j < 8; j++) {
            genaffinepairM8(&assistant->P[i][j][0], &assistant->P[i][j][1]);
        }

        MatrixcomM8to32(assistant->P[i][0][0].Mat, assistant->P[i][1][0].Mat, \
            assistant->P[i][2][0].Mat, assistant->P[i][3][0].Mat, &assistant->PC[i][0][0].Mat);

        invsM32(assistant->PC[i][0][0].Mat, &assistant->PC[i][0][1].Mat);

        VectorcomV8to32(assistant->P[i][0][0].Vec, assistant->P[i][1][0].Vec, \
            assistant->P[i][2][0].Vec, assistant->P[i][3][0].Vec, &assistant->PC[i][0][0].Vec);

        MatMulVecM32(assistant->PC[i][0][1].Mat, assistant->PC[i][0][0].Vec, &assistant->PC[i][0][1].Vec);

#ifdef AFFINE_IMPL_DEBUG
        ctx->xorInv[i][0].Mat = assistant->PC[i][0][0].Mat;
        ctx->xorInv[i][0].Vec = assistant->PC[i][0][0].Vec;
#endif 

        MatrixcomM8to32(assistant->P[i][4][0].Mat, assistant->P[i][5][0].Mat, \
            assistant->P[i][6][0].Mat, assistant->P[i][7][0].Mat, &assistant->PC[i][1][0].Mat);

        invsM32(assistant->PC[i][1][0].Mat, &assistant->PC[i][1][1].Mat);

        VectorcomV8to32(assistant->P[i][4][0].Vec, assistant->P[i][5][0].Vec, \
            assistant->P[i][6][0].Vec, assistant->P[i][7][0].Vec, &assistant->PC[i][1][0].Vec);

        MatMulVecM32(assistant->PC[i][1][1].Mat, assistant->PC[i][1][0].Vec, &assistant->PC[i][1][1].Vec);

#ifdef AFFINE_IMPL_DEBUG
        ctx->xorInv[i][1].Mat = assistant->PC[i][1][0].Mat;
        ctx->xorInv[i][1].Vec = assistant->PC[i][1][0].Vec;
#endif 
    }
    // for (i=0; i<8; i++) {
    //     genRandomAffineMatrix(  assistant->SE[0][i][0].linearMap, assistant->SE[0][i][1].linearMap, \
    //                             assistant->SE[0][i][0].vectorTranslation, assistant->SE[0][i][1].vectorTranslation, 8);
    //     genRandomAffineMatrix(  assistant->SE[1][i][0].linearMap, assistant->SE[1][i][1].linearMap, \
    //                             assistant->SE[1][i][0].vectorTranslation, assistant->SE[1][i][1].vectorTranslation, 8);
    // }
    return 0;
}


int apply_affine_matrix(klein80_wb_t* ctx, klein_wb_assistant_t* assistant) {
    int i, j, k;
    int rounds = ctx->rounds;

    uint32_t Te_T[256];

    for (i = 0; i < rounds; i++) {
        for (j = 0; j < 8; j++) {
            for (k = 0; k < 256; k++) {
                int pos = (j + 2) % 8;
                Te_T[k] = affineU32(assistant->PC[i + 1][j / 4][1], ctx->Te[i][j][affineU8(assistant->P[i][pos][0], k)]);
                if (j % 4 == 0) {
                    Te_T[k] = assistant->PC[i + 1][j / 4][1].Vec.V ^ Te_T[k];
                }
            }
            for (k = 0; k < 256; k++) {
                ctx->Te[i][j][k] = Te_T[k];
            }
        }
    }

    for (i = 0; i < 8; i++) {
        for (j = 0; j < 256; j++) {
            ctx->SE[1][i][j] = ctx->SE[0][i][affineU8(assistant->P[rounds][i][0], j)];
        }
    }

    for (i = 0; i < 8; i++) {
        for (j = 0; j < 256; j++) {
            ctx->SE[0][i][j] = affineU8(assistant->P[0][i][1], j);
        }
    }

    return 0;
}

#endif //AFFINE_IMPL

#ifdef MAT_BIJECTIONS

int    gen_gf_mat(klein80_wb_t* ctx, klein_wb_assistant_t* assistant) {
    M8 mat, mati;
    for (int i = 0; i < 8; i++) {
        initinvbaseM8(initM8_max);
        genMatpairM8(&mat, &mati);
        for (int j = 0; j < 256; j++) {
            int t = applyMatToU8(mat, j);
            assistant->fi[i][0][j] = ctx->fi[i][j] = t;
            assistant->fi[i][1][t] = j;
        }
        initinvbaseM8(initM8_max);
        genMatpairM8(&mat, &mati);
        for (int j = 0; j < 256; j++) {
            int t = applyMatToU8(mat, j);
            assistant->g[i][0][j] = ctx->g[i][j] = t;
            assistant->g[i][1][t] = j;
        }
    }
    return 0;
}

int    apply_mat_bijections(klein80_wb_t* ctx, klein_wb_assistant_t* assistant) {
    int i, j, k;
    int rounds = ctx->rounds;

    uint32_t Te_T[256];

    M8 p[8], pi[8], tp[8];
    uint32_t s[8];

    for (j = 0; j < 8; j++) {
        initinvbaseM8(initM8_max);
        genMatpairM8(&p[j], &tp[j]);
    }
    for (j = 0; j < 8; j++) {
        for (k = 0; k < 256; k++) {
            uint32_t s = ctx->Te[0][j][ctx->fi[j][k]];
            // uint32_t s = ctx->Te[0][j][ applyMatToU8(p[j], k) ];
            Te_T[k] = (applyMatToU8(p[j], (s >> 24) & 0xff) << 24) ^ \
                (applyMatToU8(p[j], (s >> 16) & 0xff) << 16) ^ \
                (applyMatToU8(p[j], (s >> 8) & 0xff) << 8) ^ \
                (applyMatToU8(p[j], (s) & 0xff));
        }
        for (k = 0; k < 256; k++) {
            ctx->Te[i][j][k] = Te_T[k];
        }
    }
    for (j = 0; j < 8; j++) {
        pi[j] = tp[j];
    }


    for (i = 1; i < rounds; i++) {

        for (j = 0; j < 8; j++) {
            initinvbaseM8(initM8_max);
            genMatpairM8(&p[j], &tp[j]);
        }
        for (j = 0; j < 8; j++) {
            for (k = 0; k < 256; k++) {
                uint32_t s = (applyMatToU8(pi[j], (s >> 24) & 0xff) << 24) ^ \
                    (applyMatToU8(pi[j], (s >> 16) & 0xff) << 16) ^ \
                    (applyMatToU8(pi[j], (s >> 8) & 0xff) << 8) ^ \
                    (applyMatToU8(pi[j], (s) & 0xff));
                // uint32_t s = ctx->Te[0][j][ applyMatToU8(p[j], k) ];
                Te_T[k] = (applyMatToU8(p[j], (s >> 24) & 0xff) << 24) ^ \
                    (applyMatToU8(p[j], (s >> 16) & 0xff) << 16) ^ \
                    (applyMatToU8(p[j], (s >> 8) & 0xff) << 8) ^ \
                    (applyMatToU8(p[j], (s) & 0xff));
            }
            for (k = 0; k < 256; k++) {
                ctx->Te[i][j][k] = Te_T[k];
            }
        }

        for (j = 0; j < 8; j++) {
            initinvbaseM8(initM8_max);
            genMatpairM8(&p[i], &tp[i]);
            for (k = 0; k < 256; k++) {
                Te_T[k] = affineU32(assistant->PC[i + 1][j / 4][1], ctx->Te[i][j][affineU8(assistant->P[i][j][0], k)]);
            }
            for (k = 0; k < 256; k++) {
                ctx->Te[i][j][k] = Te_T[k];
            }
        }

        for (j = 0; j < 8; j++) {
            pi[j] = tp[j];
        }
    }

    for (i = 0; i < 8; i++) {
        for (j = 0; j < 256; j++) {
            ctx->SE[0][i][j] = affineU8(assistant->P[0][i][1], j);
            ctx->SE[1][i][j] = affineU8(assistant->P[rounds][i][0], j);
        }
    }

    return 0;
}



#endif //MAT_BIJECTIONS


int gen_klein80wb_table_rounds(const uint8_t* key, klein80_wb_t* klein80wb_ctx, int rounds) {

    klein_wb_assistant_t assistant;
    assistant.rounds = rounds;
    gen_klein80wb_ctx_table(key, klein80wb_ctx, rounds);

#ifdef MAT_BIJECTIONS
    gen_gf_mat(klein80wb_ctx, &assistant);
    apply_mat_bijections(klein80wb_ctx, &assistant)
#endif //MAT_BIJECTIONS

#ifdef AFFINE_IMPL
        gen_affine_matrix(klein80wb_ctx, &assistant);
        apply_affine_matrix(klein80wb_ctx, &assistant);
#endif //AFFINE_IMPL

    return 0;
}


void klein80wb_encrypt(const uint8_t* plain, const klein80_wb_t* klein80wb_ctx, uint8_t* cipher) {
    int i;
    int rounds = klein80wb_ctx->rounds;
    uint8_t state[8];

#ifndef AFFINE_IMPL

    state[0] = plain[0];
    state[1] = plain[1];
    state[2] = plain[2];
    state[3] = plain[3];
    state[4] = plain[4];
    state[5] = plain[5];
    state[6] = plain[6];
    state[7] = plain[7];

#else

    //input encoding

    state[0] = klein80wb_ctx->SE[0][0][plain[0]];
    state[1] = klein80wb_ctx->SE[0][1][plain[1]];
    state[2] = klein80wb_ctx->SE[0][2][plain[2]];
    state[3] = klein80wb_ctx->SE[0][3][plain[3]];
    state[4] = klein80wb_ctx->SE[0][4][plain[4]];
    state[5] = klein80wb_ctx->SE[0][5][plain[5]];
    state[6] = klein80wb_ctx->SE[0][6][plain[6]];
    state[7] = klein80wb_ctx->SE[0][7][plain[7]];

#endif

    u32 t0, t1;



    for (i = 0; i < rounds; i++) {
        // temp_state[0] = klein80wb_ctx->T[i][2][state[2]];
        // temp_state[1] = klein80wb_ctx->T[i][3][state[3]];
        // temp_state[2] = klein80wb_ctx->T[i][4][state[4]];
        // temp_state[3] = klein80wb_ctx->T[i][5][state[5]];
        // temp_state[4] = klein80wb_ctx->T[i][6][state[6]];
        // temp_state[5] = klein80wb_ctx->T[i][7][state[7]];
        // temp_state[6] = klein80wb_ctx->T[i][0][state[0]];
        // temp_state[7] = klein80wb_ctx->T[i][1][state[1]];

        // printf("round %d:\n", i);
        // for (int j=0; j<8; j++) 
        //     printf("%02X", temp_state[j]);
        // printf("\n");

        // t0 = Te0[temp_state[0]] ^ Te1[temp_state[1]] ^ Te2[temp_state[2]] ^ Te3[temp_state[3]];
        // t1 = Te0[temp_state[4]] ^ Te1[temp_state[5]] ^ Te2[temp_state[6]] ^ Te3[temp_state[7]];

#ifdef AFFINE_IMPL_DEBUG
        printf("round %d:\n", i);
        uint32_t ttt;
        ttt = (state[0] << 24) | \
            (state[1] << 16) | \
            (state[2] << 8) | \
            (state[3]);
        ttt = affineU32(klein80wb_ctx->xorInv[i][0], ttt);
        printf("%08X", ttt);
        ttt = (state[4] << 24) | \
            (state[5] << 16) | \
            (state[6] << 8) | \
            (state[7]);
        ttt = affineU32((klein80wb_ctx->xorInv[i][1]), ttt);
        printf("%08X", ttt);
        printf("\n");

        u32 gg00, gg01, gg02, gg03;
        gg00 = affineU32((klein80wb_ctx->xorInv[i + 1][0]), klein80wb_ctx->Te[i][0][state[2]]);
        gg01 = affineU32((klein80wb_ctx->xorInv[i + 1][0]), klein80wb_ctx->Te[i][1][state[3]]);
        gg02 = affineU32((klein80wb_ctx->xorInv[i + 1][0]), klein80wb_ctx->Te[i][2][state[4]]);
        gg03 = affineU32((klein80wb_ctx->xorInv[i + 1][0]), klein80wb_ctx->Te[i][3][state[5]]);
        printf("t0:\t%08X ^ %08X ^ %08X ^ %08X = %08X\n", gg00, gg01, gg02, gg03, gg00 ^ gg01 ^ gg02 ^ gg03);
        gg00 = affineU32((klein80wb_ctx->xorInv[i + 1][1]), klein80wb_ctx->Te[i][4][state[6]]);
        gg01 = affineU32((klein80wb_ctx->xorInv[i + 1][1]), klein80wb_ctx->Te[i][5][state[7]]);
        gg02 = affineU32((klein80wb_ctx->xorInv[i + 1][1]), klein80wb_ctx->Te[i][6][state[0]]);
        gg03 = affineU32((klein80wb_ctx->xorInv[i + 1][1]), klein80wb_ctx->Te[i][7][state[1]]);
        printf("t1:\t%08X ^ %08X ^ %08X ^ %08X = %08X\n", gg00, gg01, gg02, gg03, gg00 ^ gg01 ^ gg02 ^ gg03);

#endif 

        // printf("round %d:\n", i);
        // uint32_t ttt;
        // ttt =   (state[0]<<24) | \
        //         (state[1]<<16) | \
        //         (state[2]<< 8) | \
        //         (state[3]    );
        // printf("%08X", ttt);
        // ttt =   (state[4]<<24) | \
        //         (state[5]<<16) | \
        //         (state[6]<< 8) | \
        //         (state[7]    );
        // printf("%08X", ttt);
        // printf("\n");

        // u32 gg00, gg01, gg02, gg03;
        // gg00 = klein80wb_ctx->Te[i][0][state[2]];
        // gg01 = klein80wb_ctx->Te[i][1][state[3]];
        // gg02 = klein80wb_ctx->Te[i][2][state[4]];
        // gg03 = klein80wb_ctx->Te[i][3][state[5]];
        // // klein80wb_ctx->Te[i][2][state[4]] ^ klein80wb_ctx->Te[i][3][state[5]]
        // printf("t0:\t%08X ^ %08X ^ %08X ^ %08X\n", gg00, gg01, gg02, gg03);
        // gg00 = klein80wb_ctx->Te[i][4][state[6]];
        // gg01 = klein80wb_ctx->Te[i][5][state[7]];
        // gg02 = klein80wb_ctx->Te[i][6][state[0]];
        // gg03 = klein80wb_ctx->Te[i][7][state[1]];
        // printf("t1:\t%08X ^ %08X ^ %08X ^ %08X\n", gg00, gg01, gg02, gg03);

        t0 = klein80wb_ctx->Te[i][0][state[2]] ^ klein80wb_ctx->Te[i][1][state[3]] ^ klein80wb_ctx->Te[i][2][state[4]] ^ klein80wb_ctx->Te[i][3][state[5]];
        t1 = klein80wb_ctx->Te[i][4][state[6]] ^ klein80wb_ctx->Te[i][5][state[7]] ^ klein80wb_ctx->Te[i][6][state[0]] ^ klein80wb_ctx->Te[i][7][state[1]];

        state[0] = t0 >> 24;
        state[1] = t0 >> 16;
        state[2] = t0 >> 8;
        state[3] = t0;
        state[4] = t1 >> 24;
        state[5] = t1 >> 16;
        state[6] = t1 >> 8;
        state[7] = t1;

    }

#ifdef AFFINE_IMPL

    //output encoding
    cipher[0] = klein80wb_ctx->SE[1][0][state[0]];
    cipher[1] = klein80wb_ctx->SE[1][1][state[1]];
    cipher[2] = klein80wb_ctx->SE[1][2][state[2]];
    cipher[3] = klein80wb_ctx->SE[1][3][state[3]];
    cipher[4] = klein80wb_ctx->SE[1][4][state[4]];
    cipher[5] = klein80wb_ctx->SE[1][5][state[5]];
    cipher[6] = klein80wb_ctx->SE[1][6][state[6]];
    cipher[7] = klein80wb_ctx->SE[1][7][state[7]];

#else
    cipher[0] = klein80wb_ctx->SE[0][0][state[0]];
    cipher[1] = klein80wb_ctx->SE[0][1][state[1]];
    cipher[2] = klein80wb_ctx->SE[0][2][state[2]];
    cipher[3] = klein80wb_ctx->SE[0][3][state[3]];
    cipher[4] = klein80wb_ctx->SE[0][4][state[4]];
    cipher[5] = klein80wb_ctx->SE[0][5][state[5]];
    cipher[6] = klein80wb_ctx->SE[0][6][state[6]];
    cipher[7] = klein80wb_ctx->SE[0][7][state[7]];

#endif

    return;
}



void klein80wb_decrypt(const uint8_t* cipher, const klein80_wb_t* klein80wb_ctx, uint8_t* plain) {

    return;
}
