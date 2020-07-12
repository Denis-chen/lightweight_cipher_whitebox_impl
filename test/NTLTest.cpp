/*
 * @Author: Weijie Li 
 * @Date: 2017-11-15 11:22:59 
 * @Last Modified by: Weijie Li
 * @Last Modified time: 2017-11-26 21:37:11
 */

#include <NTL/mat_gf2.h>
#include <NTL/GF2.h>
#include <NTL/vec_gf2.h>

#include <affine/NTLUtils.h>
#include <affine/MatrixUtils.h>

#define CSL32(x,i)					\
	(((x) << i) | ((x) >> (32-i)))

int NTLTest() {
    int i,j;
    printf("--------------NTL--------------\n");
    Mat<GF2> X;
    X.SetDims(2,2); 
    Mat<GF2> A,B;
    long data[] ={1,4,8,2, 7,5};
    initMatrix(X, data);
    dumpMatrix(X);

    printf("random 4 * 4 matrix and its inv: \n");
    // NTL::mat_GF2 Y;
    Mat<GF2> Y;
    // Y.SetDims(4,4);
    // initMatrixFromBit(Y, data);
    random(Y,4 ,4);
    dumpMatrix(Y);
    NTL::GF2 det = determinant(Y);
    while (det==0) {
        random(Y,4 ,4);
        det = determinant(Y);

        
    }
    Mat<GF2> Z;
    Z.SetDims(4,4);
    Z = inv(Y);
    dumpMatrix(Z);
    
    printf("------------------------------\n");
    printf("random affine augment matrix:\n");
    
    genRandomAffineAugmentedMatrix(X, Y, 32);
    // dumpMatrix(X);
    // dumpMatrix(Y);

    NTL::vec_GF2 vg;
    random(vg, 33);
    vg.put(32, 1);

    dumpVector(vg);
    printf("%u\n", get32FromVec(vg));

    NTL::vec_GF2 vr;
    vr = X*vg;
    dumpVector(vr);
    printf("%u\n", get32FromVec(vr)); 
    
    vr = Y*vr;
    dumpVector(vr);
    printf("%u\n", get32FromVec(vr));      
    
    printf("------------------------------\n");
    printf("random 2 affine augmented matrix xor before and after:\n");

    NTL::vec_GF2 va, vb, vc, v_result;
    random(va, 33);
    random(vb, 33);
    printf("vector a: %u\n", get32FromVec(va));
    printf("vector b: %u\n", get32FromVec(vb));
    
    printf("before affine:\n");
    v_result = va + vb;
    dumpVector(v_result);

    printf("after affine and deaffine:\n");
    v_result = X*va + X*vb;
    v_result = Y*v_result;
    dumpVector(v_result);

    printf("single xor between u32:\n");
    getVecFrom32(v_result, get32FromVec(va)^get32FromVec(vb));
    dumpVector(v_result);

    printf("------------------------------\n");
    printf("random affine matrix:\n");
    
    NTL::vec_GF2 V;
    genRandomAffineMatrix(X, Y, V,32);
    // dumpMatrix(X);
    // dumpMatrix(Y);
    random(vg, 32);

    dumpVector(vg);
    printf("%u\n", get32FromVec(vg));

    vr = X*vg+V;
    dumpVector(vr);
    printf("%u\n", get32FromVec(vr)); 
    
    vr = Y*(vr+V);
    dumpVector(vr);
    printf("%u\n", get32FromVec(vr));   

    printf("------------------------------\n");
    printf("random 3 affine matrix xor before and after:\n");

    random(va, 32);
    random(vb, 32);
    random(vc, 32);
    v_result = va + vb + vc;
    printf("before affine:\n");
    printf("vector a: %u\n", get32FromVec(va));
    printf("vector b: %u\n", get32FromVec(vb));
    printf("vector c: %u\n", get32FromVec(vc));
    printf("vector xor: %u\n", get32FromVec(v_result));
    printf("\n");
    
    printf("after affine:\n");
    va = X*va+V;
    vb = X*vb+V;
    vc = X*vc+V;
    printf("vector a: %u\n", get32FromVec(va));
    printf("vector b: %u\n", get32FromVec(vb));
    printf("vector c: %u\n", get32FromVec(vc));
    // dumpVector(v_result);
    printf("\n");

    printf("after deaffine:\n");
    v_result = va+vb+vc;
    v_result = Y*(v_result+V);
    printf("vector xor: %u\n", get32FromVec(v_result));

    // dumpVector(v_result);

    // dumpVector(v_result);

    printf("------------------------------\n");
    printf("random matrix diag test:\n");

    NTL::mat_GF2 ma, mb, mc;
    random(ma, 2, 2);
    random(mb, 3, 3);
    combineDiagMat(mc, ma, mb);

    printf("matrix A:\n");
    dumpMatrix(ma);
    printf("matrix B:\n");
    dumpMatrix(mb);
    printf("matrix After combined:\n");
    dumpMatrix(mc);

    printf("------------------------------\n");
    printf("CSL matrix test:\n");

    printf("x + (x<<<2) + (x<<<10) + (x<<<18) + (x<<<24) : \n");

    random(va, 32);
    long r = get32FromVec(va);
    long csl_r = (r) ^  CSL32((r),  2) ^				\
                        CSL32((r), 10) ^				\
                        CSL32((r), 18) ^				\
                        CSL32((r), 24);
    printf("random int:\t%X\n", r);
    dumpVector(va);
    printf("after CSL:\t%X\n", csl_r);
    X.SetDims(1,32);
    initMatrixFromBit(X, &csl_r);
    dumpMatrix(X);
    printf("\n");
    long csl_m =    (1<<31) ^ \
                    (1<<29) ^ \
                    (1<<21) ^ \
                    (1<<13) ^ \
                    (1<<7);
    long csl_ms[32];
    csl_ms[0] = csl_m;
    for (i=1; i<32; i++) {
        long t = csl_ms[i-1];
        csl_ms[i] = (((t>> 1) & 0x7fffffff) | (t << (31)));
    } 
    // for (i=0; i<32; i++) {
    //     printf("0x%4X, ", csl_ms[i]);
    // }
    // printf("\n");
                    
    ma.SetDims(32,32);
    initMatrixFromBit(ma, csl_ms);
    // dumpMatrix(ma);

    va = ma*va;
    r = get32FromVec(va);
    printf("random int after CSL matrix:\t%X\n", r);
    dumpVector(va);
    
    printf("CSL matrix:\n");
    dumpMatrix(ma);

    printf("------------------------------\n");
    printf("Diag matrix xor test:\n");
    printf("random 3 affine diag matrix xor before and after:\n");

    random(va, 32);
    random(vb, 32);
    random(vc, 32);
    v_result = va + vb + vc;
    printf("before affine:\n");
    printf("vector a: %u\n", get32FromVec(va));
    printf("vector b: %u\n", get32FromVec(vb));
    printf("vector c: %u\n", get32FromVec(vc));
    printf("vector xor: %u\n", get32FromVec(v_result));
    printf("\n");
    
    mat_GF2 Xi[4];
    for (i=0; i<4; i++) {
        genRandomInvMatrix(Xi[i], X, 8);
    }
    mat_GF2 s, t; 
    combineDiagMat(s, Xi[0], Xi[1]);
    combineDiagMat(t, s, Xi[2]);
    combineDiagMat(X, t, Xi[3]);
    // dumpMatrix(X);
    printf("after affine:\n");
    long a,b,c;
    a = get32FromVec(va);
    b = get32FromVec(vb);
    c = get32FromVec(vc);
    vec_GF2 vt;
    initVecFromBit(vt, a>>24, 8);
    va = Xi[0]*vt;
    initVecFromBit(vt, (a>>16 & 0xff), 8);
    va.append( Xi[1]*vt); 
    initVecFromBit(vt, (a>>8 & 0xff), 8);
    va.append( Xi[2]*vt); 
    initVecFromBit(vt, (a    & 0xff), 8);
    va.append( Xi[3]*vt); 

    initVecFromBit(vt, b >> 24, 8);
    vb = Xi[0]*vt;
    initVecFromBit(vt, (b >> 16 & 0xff), 8);
    vb.append( Xi[1]*vt); 
    initVecFromBit(vt, (b >> 8 & 0xff), 8);
    vb.append( Xi[2]*vt); 
    initVecFromBit(vt, (b    & 0xff), 8);
    vb.append( Xi[3]*vt); 

    initVecFromBit(vt, c >>24, 8);
    vc = Xi[0]*vt;
    initVecFromBit(vt, (c >>16 & 0xff), 8);
    vc.append( Xi[1]*vt); 
    initVecFromBit(vt, (c >>8 & 0xff), 8);
    vc.append( Xi[2]*vt); 
    initVecFromBit(vt, (c    & 0xff), 8);
    vc.append( Xi[3]*vt); 
    
    printf("vector a: %u\n", get32FromVec(va));
    printf("vector b: %u\n", get32FromVec(vb));
    printf("vector c: %u\n", get32FromVec(vc));
    // dumpVector(v_result);
    printf("\n");

    printf("after deaffine:\n");
    v_result = va+vb+vc;
    v_result = inv(X)*(v_result);
    printf("vector xor: %u\n", get32FromVec(v_result));


    



    return 0;
}

