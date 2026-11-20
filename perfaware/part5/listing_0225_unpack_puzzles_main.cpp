/* ========================================================================

   (C) Copyright 2026 by Molly Rocket, Inc., All Rights Reserved.
   
   This software is provided 'as-is', without any express or implied
   warranty. In no event will the authors be held liable for any damages
   arising from the use of this software.
   
   Please see https://computerenhance.com for more information
   
   ======================================================================== */

/* ========================================================================
   LISTING 225
   ======================================================================== */

#include "listing_0220_puzzle_utils.cpp"
 
#define SHUFIMM_2x4(E0, E1) ((E0) | ((E1) << 4))
#define SHUFIMM_4x2(E0, E1, E2, E3) ((E0) | ((E1) << 2) | ((E2) << 4) | ((E3) << 6))
#define SHUFIMM_4x1(E0, E1, E2, E3) ((E0) | ((E1) << 1) | ((E2) << 2) | ((E3) << 3))

static m256 Puzzle0(m256 A, m256 B)
{
    // "abcd" "efgh" "ijkl" "mnop" "qrst" "uvwx" "1234" "5678";
    // "ABCD" "EFGH" "IJKL" "MNOP" "QRST" "UVWX" ".,;-" "=[]|";

    // ijIJ" "klKL" "mnMN" "opOP" "12.," "34;-" "56=[" "78]|"
    __m256 Result = _mm256_unpackhi_epi16(A.I, B.I);
    
    return M256From(Result);
}

static m256 Puzzle1(m256 A, m256 B)
{
    // "abcd" "efgh" "ijkl" "mnop" | "qrst" "uvwx" "1234" "5678";
    // "ABCD" "EFGH" "IJKL" "MNOP" | "QRST" "UVWX" ".,;-" "=[]|";
    // "aAiI" "bBjJ" "cCkK" "dDlL" "qQ1." "rR2," "sS3;" "tT4-"

    // aibj ckdl <> <> q1r2 s3t4
    // AIBJ CKDL <> <> Q.R, S;T-
    

    // X = unpackLo8(A,B): aAbB cCdD eEfF gGhH | qQrR sStT uUvV wWxX
    // Y = unpackHi8(A,B): iIjJ kKlL mMnN oOpP | 1.2, 3;4- 5=6[ 7]8|

    // unpackLo16(X,Y): aAiI bBjJ cCkK dDlL | qQ1. rR2, sS3; tT4- 

    // Is there a better way than epi8_shuffle?.....
    // UPDATE: Oct 22 2026: YES! There is a better way!
    
    /* Old way, Oct 20: 
    __m256 indices = _mm256_setr_epi8(
        0, 8, 1, 9, 2, 10, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1,
        0, 8, 1, 9, 2, 10, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1
    );
    __m256i intermediateA = _mm256_shuffle_epi8(
        A.I, indices
    );
    __m256i intermediateB = _mm256_shuffle_epi8(
        B.I, indices
    );
    */

    // New Way:
    __m256 X = _mm256_unpacklo_epi8(A.I, B.I);
    __m256 Y = _mm256_unpackhi_epi8(A.I, B.I);

    __m256 Result = _mm256_unpacklo_epi16(X, Y);
    
    return M256From(Result);
}

static m256 Puzzle2(m256 A, m256 B)
{
    // "abcd" "efgh" "ijkl" "mnop" | "qrst" "uvwx" "1234" "5678";
    // "ABCD" "EFGH" "IJKL" "MNOP" | "QRST" "UVWX" ".,;-" "=[]|";
    
    // "1i2j" "3k4l" "5m6n" "7o8p" | ".1,2" ";3-4" "=5[6" "]7|8"

    // "abcd" "efgh" "ijkl" "mnop" | "qrst" "uvwx" "1234" "5678";
    // "qrst" "uvwx" "1234" "5678" | "QRST" "UVWX" ".,;-" "=[]|";
    
    __m256 merged = _mm256_permute2x128_si256(A.I, B.I, 0b110001);

    

    // 1234 5678

    __m256 Result = _mm256_unpackhi_epi8(merged, A.I);
    
    return M256From(Result);
}

int main(void)
{
    char ValuesA[SIMD_BYTE_COUNT + 1] = "abcd" "efgh" "ijkl" "mnop" "qrst" "uvwx" "1234" "5678";
    char ValuesB[SIMD_BYTE_COUNT + 1] = "ABCD" "EFGH" "IJKL" "MNOP" "QRST" "UVWX" ".,;-" "=[]|";
    
    m256 A, B;
    A.D = _mm256_loadu_pd((f64 *)ValuesA);
    B.D = _mm256_loadu_pd((f64 *)ValuesB);
    
    //
    // NOTE(casey): Unpack puzzles
    //
    
    Match("Puzzle0", Puzzle0(A, B), "ijIJ" "klKL" "mnMN" "opOP" "12.," "34;-" "56=[" "78]|");
    Match("Puzzle1", Puzzle1(A, B), "aAiI" "bBjJ" "cCkK" "dDlL" "qQ1." "rR2," "sS3;" "tT4-");
    Match("Puzzle2", Puzzle2(A, B), "1i2j" "3k4l" "5m6n" "7o8p" ".1,2" ";3-4" "=5[6" "]7|8");
    
    return 0;
}