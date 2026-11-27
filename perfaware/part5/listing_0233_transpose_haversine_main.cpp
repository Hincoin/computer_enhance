/* ========================================================================

   (C) Copyright 2026 by Molly Rocket, Inc., All Rights Reserved.
   
   This software is provided 'as-is', without any express or implied
   warranty. In no event will the authors be held liable for any damages
   arising from the use of this software.
   
   Please see https://computerenhance.com for more information
   
   ======================================================================== */

/* ========================================================================
   LISTING 233
   ======================================================================== */

#include "listing_0208_combined_haversine_test.cpp"
#include "listing_0218_wide_haversine_math.cpp"

#include <assert.h>

static __attribute__((noinline)) f64 CompilerLoadedHaversine(haversine_setup Setup)
{
    assert((Setup.PairCount % 4) == 0);
    
    u64 PairCount = Setup.PairCount;
    haversine_pair *Pairs = Setup.Pairs;
    
    __m256d Sum = _mm256_setzero_pd();
    
    __m256d SumCoef = _mm256_set1_pd((2.0*QUESTIONABLE_EARTH_RADIUS) / (f64)PairCount);
    for(u64 PairIndex = 0; PairIndex < PairCount; PairIndex += 4)
    {
        asm volatile("nop\nnop\nnop"::: "memory");
        haversine_pair Pair0 = Pairs[PairIndex + 0];
        haversine_pair Pair1 = Pairs[PairIndex + 1];
        haversine_pair Pair2 = Pairs[PairIndex + 2];
        haversine_pair Pair3 = Pairs[PairIndex + 3];
        
        __m256d lat1 = _mm256_setr_pd(Pair0.Y0, Pair1.Y0, Pair2.Y0, Pair3.Y0);
        __m256d lat2 = _mm256_setr_pd(Pair0.Y1, Pair1.Y1, Pair2.Y1, Pair3.Y1);
        __m256d lon1 = _mm256_setr_pd(Pair0.X0, Pair1.X0, Pair2.X0, Pair3.X0);
        __m256d lon2 = _mm256_setr_pd(Pair0.X1, Pair1.X1, Pair2.X1, Pair3.X1);
        asm volatile("nop\nnop\nnop"::: "memory");

        __m256d NegRadC = _mm256_set1_pd(-0.01745329251994329577);
        __m256d NegPi = _mm256_set1_pd(-3.14159265358979323846);
        __m256d HalfRadC = _mm256_set1_pd(0.00872664625997164788);
        __m256d HalfPi = _mm256_set1_pd(1.57079632679489661923);
        __m256d Half = _mm256_set1_pd(0.5);
        __m256d One = _mm256_set1_pd(1.0);
        __m256d Deg180 = _mm256_set1_pd(180.0);
        
        __m256d AbsDLat = FAbs4x(_mm256_sub_pd(lat2, lat1));
        __m256d AbsDLon = FAbs4x(_mm256_sub_pd(lon2, lon1));
        __m256d AbsLat1 = FAbs4x(lat1);
        __m256d AbsLat2 = FAbs4x(lat2);
        
        __m256d ALC0Pred = _mm256_cmp_pd(AbsDLat, Deg180, _CMP_LT_OS);
        __m256d ALC0 = _mm256_andnot_pd(ALC0Pred, NegPi);
        
        __m256d ALC3Pred = _mm256_cmp_pd(AbsDLon, Deg180, _CMP_LT_OS);
        __m256d ALC3 = _mm256_andnot_pd(ALC3Pred, NegPi);
        
        __m256d S1 = SineCoreWithPrefix4x(NegRadC, AbsLat1, HalfPi);
        __m256d S2 = SineCoreWithPrefix4x(NegRadC, AbsLat2, HalfPi);
        __m256d S0 = SineCoreWithPrefix4x(HalfRadC, AbsDLat, ALC0);
        __m256d S3 = SineCoreWithPrefix4x(HalfRadC, AbsDLon, ALC3);
        
        __m256d S12 = _mm256_mul_pd(S1, S2);
        __m256d S33 = _mm256_mul_pd(S3, S3);
        __m256d S1233 = _mm256_mul_pd(S12, S33);
        __m256d a = _mm256_fmadd_pd(S0, S0, S1233);
        
        __m256d NeedsTransform = _mm256_cmp_pd(Half, a, _CMP_LT_OS);
        __m256d RangeA = _mm256_blendv_pd(a, _mm256_sub_pd(One, a), NeedsTransform);
        __m256d R = ArcsineCoreFromSquared4x(RangeA);
        __m256d RangeR = _mm256_blendv_pd(R, _mm256_sub_pd(HalfPi, R), NeedsTransform);
        
        Sum = _mm256_fmadd_pd(SumCoef, RangeR, Sum);
    }
    
    __m256d SumSwap = _mm256_permute4x64_pd(Sum, 2 + (3 << 2));
    __m256d SumHalf = _mm256_add_pd(Sum, SumSwap);
    __m256d SumHalfSwap = _mm256_permute_pd(SumHalf, 1);
    __m256d SumFull = _mm256_add_pd(SumHalf, SumHalfSwap);
    
    return _mm256_cvtsd_f64(SumFull);
}

static __attribute__((noinline)) f64 GatherHaversine(haversine_setup Setup)
{
    assert((Setup.PairCount % 4) == 0);
    
    u64 PairCount = Setup.PairCount;
    haversine_pair *Pairs = Setup.Pairs;
    
    __m256d Sum = _mm256_setzero_pd();
    
    __m256d SumCoef = _mm256_set1_pd((2.0*QUESTIONABLE_EARTH_RADIUS) / (f64)PairCount);
    for(u64 PairIndex = 0; PairIndex < PairCount; PairIndex += 4)
    {
        asm volatile("nop\nnop\nnop":::);
        f64 *Values = (f64 *)&Pairs[PairIndex];
        
        __m128i Indexes = _mm_setr_epi32(0, 4, 8, 12);
        __m256d lon1 = _mm256_i32gather_pd(Values + 0, Indexes, 8);
        __m256d lat1 = _mm256_i32gather_pd(Values + 1, Indexes, 8);
        __m256d lon2 = _mm256_i32gather_pd(Values + 2, Indexes, 8);
        __m256d lat2 = _mm256_i32gather_pd(Values + 3, Indexes, 8);
        asm volatile("nop\nnop\nnop":::);
        __m256d NegRadC = _mm256_set1_pd(-0.01745329251994329577);
        __m256d NegPi = _mm256_set1_pd(-3.14159265358979323846);
        __m256d HalfRadC = _mm256_set1_pd(0.00872664625997164788);
        __m256d HalfPi = _mm256_set1_pd(1.57079632679489661923);
        __m256d Half = _mm256_set1_pd(0.5);
        __m256d One = _mm256_set1_pd(1.0);
        __m256d Deg180 = _mm256_set1_pd(180.0);
        
        __m256d AbsDLat = FAbs4x(_mm256_sub_pd(lat2, lat1));
        __m256d AbsDLon = FAbs4x(_mm256_sub_pd(lon2, lon1));
        __m256d AbsLat1 = FAbs4x(lat1);
        __m256d AbsLat2 = FAbs4x(lat2);
        
        __m256d ALC0Pred = _mm256_cmp_pd(AbsDLat, Deg180, _CMP_LT_OS);
        __m256d ALC0 = _mm256_andnot_pd(ALC0Pred, NegPi);
        
        __m256d ALC3Pred = _mm256_cmp_pd(AbsDLon, Deg180, _CMP_LT_OS);
        __m256d ALC3 = _mm256_andnot_pd(ALC3Pred, NegPi);
        
        __m256d S1 = SineCoreWithPrefix4x(NegRadC, AbsLat1, HalfPi);
        __m256d S2 = SineCoreWithPrefix4x(NegRadC, AbsLat2, HalfPi);
        __m256d S0 = SineCoreWithPrefix4x(HalfRadC, AbsDLat, ALC0);
        __m256d S3 = SineCoreWithPrefix4x(HalfRadC, AbsDLon, ALC3);
        
        __m256d S12 = _mm256_mul_pd(S1, S2);
        __m256d S33 = _mm256_mul_pd(S3, S3);
        __m256d S1233 = _mm256_mul_pd(S12, S33);
        __m256d a = _mm256_fmadd_pd(S0, S0, S1233);
        
        __m256d NeedsTransform = _mm256_cmp_pd(Half, a, _CMP_LT_OS);
        __m256d RangeA = _mm256_blendv_pd(a, _mm256_sub_pd(One, a), NeedsTransform);
        __m256d R = ArcsineCoreFromSquared4x(RangeA);
        __m256d RangeR = _mm256_blendv_pd(R, _mm256_sub_pd(HalfPi, R), NeedsTransform);
        
        Sum = _mm256_fmadd_pd(SumCoef, RangeR, Sum);
    }
    
    __m256d SumSwap = _mm256_permute4x64_pd(Sum, 2 + (3 << 2));
    __m256d SumHalf = _mm256_add_pd(Sum, SumSwap);
    __m256d SumHalfSwap = _mm256_permute_pd(SumHalf, 1);
    __m256d SumFull = _mm256_add_pd(SumHalf, SumHalfSwap);
    
    return _mm256_cvtsd_f64(SumFull);
}

static __attribute__((noinline)) f64 TranposeHaversine(haversine_setup Setup)
{
    assert((Setup.PairCount % 4) == 0);
    
    u64 PairCount = Setup.PairCount;
    haversine_pair *Pairs = Setup.Pairs;
    
    __m256d Sum = _mm256_setzero_pd();
    
    __m256d SumCoef = _mm256_set1_pd((2.0*QUESTIONABLE_EARTH_RADIUS) / (f64)PairCount);
    for(u64 PairIndex = 0; PairIndex < PairCount; PairIndex += 4)
    {
        f64 *Values = (f64 *)&Pairs[PairIndex];
        asm volatile("nop\nnop\nnop":::);
        __m256d Row0 = _mm256_loadu_pd(Values);
        __m256d Row1 = _mm256_loadu_pd(Values + 4);
        __m256d Row2 = _mm256_loadu_pd(Values + 8);
        __m256d Row3 = _mm256_loadu_pd(Values + 12);
        
        __m256d Lo01 = _mm256_unpacklo_pd(Row0, Row1);
        __m256d Hi01 = _mm256_unpackhi_pd(Row0, Row1);
        
        __m256d Lo23 = _mm256_unpacklo_pd(Row2, Row3);
        __m256d Hi23 = _mm256_unpackhi_pd(Row2, Row3);
        
        __m256d lon1 = _mm256_permute2f128_pd(Lo01, Lo23, 0x20);
        __m256d lat1 = _mm256_permute2f128_pd(Hi01, Hi23, 0x20);
        __m256d lon2 = _mm256_permute2f128_pd(Lo01, Lo23, 0x31);
        __m256d lat2 = _mm256_permute2f128_pd(Hi01, Hi23, 0x31);
        asm volatile("nop\nnop\nnop":::);
        __m256d NegRadC = _mm256_set1_pd(-0.01745329251994329577);
        __m256d NegPi = _mm256_set1_pd(-3.14159265358979323846);
        __m256d HalfRadC = _mm256_set1_pd(0.00872664625997164788);
        __m256d HalfPi = _mm256_set1_pd(1.57079632679489661923);
        __m256d Half = _mm256_set1_pd(0.5);
        __m256d One = _mm256_set1_pd(1.0);
        __m256d Deg180 = _mm256_set1_pd(180.0);
        
        __m256d AbsDLat = FAbs4x(_mm256_sub_pd(lat2, lat1));
        __m256d AbsDLon = FAbs4x(_mm256_sub_pd(lon2, lon1));
        __m256d AbsLat1 = FAbs4x(lat1);
        __m256d AbsLat2 = FAbs4x(lat2);
        
        __m256d ALC0Pred = _mm256_cmp_pd(AbsDLat, Deg180, _CMP_LT_OS);
        __m256d ALC0 = _mm256_andnot_pd(ALC0Pred, NegPi);
        
        __m256d ALC3Pred = _mm256_cmp_pd(AbsDLon, Deg180, _CMP_LT_OS);
        __m256d ALC3 = _mm256_andnot_pd(ALC3Pred, NegPi);
        
        __m256d S1 = SineCoreWithPrefix4x(NegRadC, AbsLat1, HalfPi);
        __m256d S2 = SineCoreWithPrefix4x(NegRadC, AbsLat2, HalfPi);
        __m256d S0 = SineCoreWithPrefix4x(HalfRadC, AbsDLat, ALC0);
        __m256d S3 = SineCoreWithPrefix4x(HalfRadC, AbsDLon, ALC3);
        
        __m256d S12 = _mm256_mul_pd(S1, S2);
        __m256d S33 = _mm256_mul_pd(S3, S3);
        __m256d S1233 = _mm256_mul_pd(S12, S33);
        __m256d a = _mm256_fmadd_pd(S0, S0, S1233);
        
        __m256d NeedsTransform = _mm256_cmp_pd(Half, a, _CMP_LT_OS);
        __m256d RangeA = _mm256_blendv_pd(a, _mm256_sub_pd(One, a), NeedsTransform);
        __m256d R = ArcsineCoreFromSquared4x(RangeA);
        __m256d RangeR = _mm256_blendv_pd(R, _mm256_sub_pd(HalfPi, R), NeedsTransform);
        
        Sum = _mm256_fmadd_pd(SumCoef, RangeR, Sum);
    }
    
    __m256d SumSwap = _mm256_permute4x64_pd(Sum, 2 + (3 << 2));
    __m256d SumHalf = _mm256_add_pd(Sum, SumSwap);
    __m256d SumHalfSwap = _mm256_permute_pd(SumHalf, 1);
    __m256d SumFull = _mm256_add_pd(SumHalf, SumHalfSwap);
    
    return _mm256_cvtsd_f64(SumFull);
}


static __attribute__((noinline)) f64 CompilerCustomLoadedHaversine(haversine_setup Setup)
{
    assert((Setup.PairCount % 4) == 0);
    
    u64 PairCount = Setup.PairCount;
    haversine_pair *Pairs = Setup.Pairs;
    
    __m256d Sum = _mm256_setzero_pd();
    
    __m256d SumCoef = _mm256_set1_pd((2.0*QUESTIONABLE_EARTH_RADIUS) / (f64)PairCount);
    for(u64 PairIndex = 0; PairIndex < PairCount; PairIndex += 4)
    {
        asm volatile("nop\nnop\nnop"::: "memory");
        haversine_pair Pair0 = Pairs[PairIndex + 0];
        haversine_pair Pair1 = Pairs[PairIndex + 1];
        haversine_pair Pair2 = Pairs[PairIndex + 2];
        haversine_pair Pair3 = Pairs[PairIndex + 3];

        __m128d xmm3 = _mm_load_sd(&Pair2.Y0);
        __m128d xmm4 = _mm_load_sd(&Pair0.X0);
        __m128d xmm6 = _mm_load_sd(&Pair0.Y0);
        __m128d xmm1 = _mm_load_sd(&Pair0.X1);

        xmm3 = _mm_loadh_pd(xmm3, &Pair3.Y0); // xmm3 = [Y0_2,  Y0_3]
        __m128d xmm5 = _mm_load_sd(&Pair0.Y1);
        xmm6 = _mm_loadh_pd(xmm6, &Pair1.Y0); // [Y0_0, Y0_1]
        __m128d xmm2 = _mm_load_sd(&Pair2.Y1);
        xmm2 = _mm_loadh_pd(xmm2, &Pair3.Y1); // [Y1_2, Y1_3]
        __m256d ymm6 = _mm256_setr_m128d(xmm6, xmm3); // Y0_0, Y0_1, Y0_2, Y0_3
        xmm3 = _mm_loadh_pd(xmm5, &Pair1.Y1); //[Y1_0,  Y1_1]
        xmm5 = _mm_load_sd(&Pair2.X0);
        xmm5 = _mm_loadh_pd(xmm5, &Pair3.X0); // [X0_2, X0_3]
        __m256d ymm2 = _mm256_setr_m128d(xmm3, xmm2); //Y1_0, Y1_1, Y1_2, Y1_3

        xmm3 = _mm_loadh_pd(xmm4, &Pair1.X0); //[X0_0,  X0_1]
        xmm4 = _mm_load_sd(&Pair2.X1);
        xmm4 = _mm_loadh_pd(xmm4, &Pair3.X1); //[X1_2,  X1_3]
        __m256d ymm3 = _mm256_setr_m128d(xmm3, xmm5); //X0_0, X0_1, X0_2, X0_3
        xmm1 = _mm_loadh_pd(xmm1, &Pair1.X1); //[X1_0,  X1_1]
        __m256d ymm1 = _mm256_setr_m128d(xmm1, xmm4); // X1_0, X1_1, X1_2, X1_3



        

        
        __m256d lat1 = ymm6; //_mm256_setr_pd(Pair0.Y0, Pair1.Y0, Pair2.Y0, Pair3.Y0);
        __m256d lat2 = ymm2; //_mm256_setr_pd(Pair0.Y1, Pair1.Y1, Pair2.Y1, Pair3.Y1);
        __m256d lon1 = ymm3; // _mm256_setr_pd(Pair0.X0, Pair1.X0, Pair2.X0, Pair3.X0);
        __m256d lon2 = ymm1; // _mm256_setr_pd(Pair0.X1, Pair1.X1, Pair2.X1, Pair3.X1);
        asm volatile("nop\nnop\nnop"::: "memory");

        __m256d NegRadC = _mm256_set1_pd(-0.01745329251994329577);
        __m256d NegPi = _mm256_set1_pd(-3.14159265358979323846);
        __m256d HalfRadC = _mm256_set1_pd(0.00872664625997164788);
        __m256d HalfPi = _mm256_set1_pd(1.57079632679489661923);
        __m256d Half = _mm256_set1_pd(0.5);
        __m256d One = _mm256_set1_pd(1.0);
        __m256d Deg180 = _mm256_set1_pd(180.0);
        
        __m256d AbsDLat = FAbs4x(_mm256_sub_pd(lat2, lat1));
        __m256d AbsDLon = FAbs4x(_mm256_sub_pd(lon2, lon1));
        __m256d AbsLat1 = FAbs4x(lat1);
        __m256d AbsLat2 = FAbs4x(lat2);
        
        __m256d ALC0Pred = _mm256_cmp_pd(AbsDLat, Deg180, _CMP_LT_OS);
        __m256d ALC0 = _mm256_andnot_pd(ALC0Pred, NegPi);
        
        __m256d ALC3Pred = _mm256_cmp_pd(AbsDLon, Deg180, _CMP_LT_OS);
        __m256d ALC3 = _mm256_andnot_pd(ALC3Pred, NegPi);
        
        __m256d S1 = SineCoreWithPrefix4x(NegRadC, AbsLat1, HalfPi);
        __m256d S2 = SineCoreWithPrefix4x(NegRadC, AbsLat2, HalfPi);
        __m256d S0 = SineCoreWithPrefix4x(HalfRadC, AbsDLat, ALC0);
        __m256d S3 = SineCoreWithPrefix4x(HalfRadC, AbsDLon, ALC3);
        
        __m256d S12 = _mm256_mul_pd(S1, S2);
        __m256d S33 = _mm256_mul_pd(S3, S3);
        __m256d S1233 = _mm256_mul_pd(S12, S33);
        __m256d a = _mm256_fmadd_pd(S0, S0, S1233);
        
        __m256d NeedsTransform = _mm256_cmp_pd(Half, a, _CMP_LT_OS);
        __m256d RangeA = _mm256_blendv_pd(a, _mm256_sub_pd(One, a), NeedsTransform);
        __m256d R = ArcsineCoreFromSquared4x(RangeA);
        __m256d RangeR = _mm256_blendv_pd(R, _mm256_sub_pd(HalfPi, R), NeedsTransform);
        
        Sum = _mm256_fmadd_pd(SumCoef, RangeR, Sum);
    }
    
    __m256d SumSwap = _mm256_permute4x64_pd(Sum, 2 + (3 << 2));
    __m256d SumHalf = _mm256_add_pd(Sum, SumSwap);
    __m256d SumHalfSwap = _mm256_permute_pd(SumHalf, 1);
    __m256d SumFull = _mm256_add_pd(SumHalf, SumHalfSwap);
    
    return _mm256_cvtsd_f64(SumFull);
}



static test_function TestFunctions[] =
{
  //  {"ReferenceHaversine", ReferenceSumHaversine},
    {"CompilerLoadedHaversine", CompilerLoadedHaversine},
   // {"GatherHaversine", GatherHaversine},
   // {"TranposeHaversine", TranposeHaversine},
    {"CompilerCUSTOMHaversine", CompilerCustomLoadedHaversine},
};

int main(int ArgCount, char **Args)
{
    CombinedHaversineTest(ArgCount, Args, ArrayCount(TestFunctions), TestFunctions);
    return 0;
}


/*
    344c:       c5 fb 10 59 d0          vmovsd -0x30(%rcx),%xmm3
    3451:       c5 fb 10 61 88          vmovsd -0x78(%rcx),%xmm4
    3456:       c5 e1 16 59 f0          vmovhpd -0x10(%rcx),%xmm3,%xmm3
    345b:       c5 fb 10 71 90          vmovsd -0x70(%rcx),%xmm6
    3460:       c5 fb 10 49 98          vmovsd -0x68(%rcx),%xmm1
    3465:       c5 fb 10 69 a0          vmovsd -0x60(%rcx),%xmm5
    346a:       c5 c9 16 71 b0          vmovhpd -0x50(%rcx),%xmm6,%xmm6
    346f:       c4 e3 4d 18 f3 01       vinsertf128 $0x1,%xmm3,%ymm6,%ymm6
    3475:       c5 fb 10 59 c8          vmovsd -0x38(%rcx),%xmm3
    347a:       c5 fb 10 51 e0          vmovsd -0x20(%rcx),%xmm2
    347f:       c5 e9 16 11             vmovhpd (%rcx),%xmm2,%xmm2
    3483:       c5 d1 16 69 c0          vmovhpd -0x40(%rcx),%xmm5,%xmm5
    3488:       c4 e3 55 18 d2 01       vinsertf128 $0x1,%xmm2,%ymm5,%ymm2
    348e:       c5 d9 16 61 a8          vmovhpd -0x58(%rcx),%xmm4,%xmm4
    3493:       c5 fb 10 69 d8          vmovsd -0x28(%rcx),%xmm5
    3498:       c5 d1 16 69 f8          vmovhpd -0x8(%rcx),%xmm5,%xmm5
    349d:       c5 f1 16 49 b8          vmovhpd -0x48(%rcx),%xmm1,%xmm1
    34a2:       c5 f1 5c cc             vsubpd %xmm4,%xmm1,%xmm1
    34a6:       c5 d1 5c db             vsubpd %xmm3,%xmm5,%xmm3
    34aa:       c4 e3 75 18 c3 01       vinsertf128 $0x1,%xmm3,%ymm1,%ymm0*/