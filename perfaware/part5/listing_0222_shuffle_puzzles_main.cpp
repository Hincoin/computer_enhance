/* ========================================================================

   (C) Copyright 2026 by Molly Rocket, Inc., All Rights Reserved.

   This software is provided 'as-is', without any express or implied
   warranty. In no event will the authors be held liable for any damages
   arising from the use of this software.

   Please see https://computerenhance.com for more information

   ======================================================================== */

/* ========================================================================
   LISTING 222
   ======================================================================== */

#include "listing_0220_puzzle_utils.cpp"

#define SHUFIMM_2x4(E0, E1) ((E0) | ((E1) << 4))
#define SHUFIMM_4x2(E0, E1, E2, E3) ((E0) | ((E1) << 2) | ((E2) << 4) | ((E3) << 6))
#define SHUFIMM_4x1(E0, E1, E2, E3) ((E0) | ((E1) << 1) | ((E2) << 2) | ((E3) << 3))

static m256 Puzzle0(m256 A)
{
    //"abcd" "efgh" "ijkl" "mnop" "qrst" "uvwx" "1234" "5678";
    //"ijkl" "mnop" "abcd" "efgh" "1234" "5678" "qrst" "uvwx");
    //  2 3 0 1
    __m256 Result = _mm256_shuffle_epi32(A.I, 0b01'00'11'10);

    return M256From(Result);
}

static m256 Puzzle1(m256 A)
{
    //"abcd" "efgh" "ijkl" "mnop" "qrst" "uvwx" "1234" "5678";
    // "5678" "qrst" "ijkl" "uvwx" "1234" "abcd" "efgh" "mnop"
    // (7, 4, 2, 5, 6, 0, 1, 3)
    __m256i indices = _mm256_setr_epi32(7, 4, 2, 5, 6, 0, 1, 3);
    __m256 Result = _mm256_permutevar8x32_ps(A.P, indices);

    return M256From(Result);
}

static m256 Puzzle2(m256 A)
{
    // "abcd" "efgh" "ijkl" "mnop" "qrst" "uvwx" "1234" "5678"
    // "keia" "dbcm" "fglj" "nohp" "sr5t" "uq8w" "x123" "4v67"
    __m256i indices = _mm256_setr_epi8(
        10, 4, 8, 0, 3, 1, 2, 12, 5, 6, 11, 9, 13, 14, 7, 15,
        2, 1, 12, 3, 4, 0, 15, 6, 7, 8, 9, 10, 11, 5, 13, 14);
    __m256 Result = _mm256_shuffle_epi8(A.I, indices);

    return M256From(Result);
}

static m256 Puzzle3(m256 A)
{
    // "abcd" "efgh" "ijkl" "mnop" "qrst" "uvwx" "1234" "5678"
    // "3fch" "1267" "ab8d" "45eg"       "uxlq" "vmin" "opwj" "ktrs"
    // "1234" "5678" "abcd" "efgh"    "ijkl" "uvwx" "mnop" "qrst"
    __m256i indices1 = _mm256_setr_epi32(6, 7, 0, 1, 2, 5, 3, 4);
    __m256i shuffle1 = _mm256_permutevar8x32_ps(A.P, indices1);
    // return M256From(shuffle1);
    __m256i indices2 = _mm256_setr_epi8(2, 13, 10, 15, 0, 1, 5, 6, 8, 9, 7, 11, 3, 4, 12, 14, 4, 7, 3, 12, 5, 8, 0, 9, 10, 11, 6, 1, 2, 15, 13, 14);
    __m256 Result = _mm256_shuffle_epi8(shuffle1, indices2);
    return M256From(Result);

    // return M256From(Result);
}

static m256 Puzzle4(m256 A, m256 B)
{
    // A = "abcd" "efgh" "ijkl" "mnop" "qrst" "uvwx" "1234" "5678"
    // B = "ABCD" "EFGH" "IJKL" "MNOP" "QRST" "UVWX" ".,;-" "=[]|"

    // "ijkl" "mnop" "ABCD" "EFGH" "qrst" "uvwx" ".,;-" "=[]|"

   //64-bit shuffle_pd
   // 1, 0, 0, 1

    __m256 Result = _mm256_shuffle_pd(A.D, B.D, 0b1001);

    return M256From(Result);
}

static m256 Puzzle5(m256 A, m256 B)
{
    // A = "abcd" "efgh" "ijkl" "mnop" | "qrst" "uvwx" "1234" "5678"
    // B = "ABCD" "EFGH" "IJKL" "MNOP" | "QRST" "UVWX" ".,;-" "=[]|"

    
    // "ijkl" "abcd" "MNOP" "EFGH"  |   "1234" "qrst" "=[]|" "UVWX"

    // AA, BB 
    // 2, 0, 3, 1
    __m256 Result =  _mm256_shuffle_ps(A.P, B.P, 0b01'11'00'10);
    return M256From(Result);
}

static m256 Puzzle6(m256 A, m256 B)
{
    // A = "abcd" "efgh" "ijkl" "mnop" | "qrst" "uvwx" "1234" "5678"
    // B = "ABCD" "EFGH" "IJKL" "MNOP" | "QRST" "UVWX" ".,;-" "=[]|"

    __m256 lowerHalves = _mm256_permute2x128_si256(A.I, B.I, 0b0001'0011);
    // QRSTUVWX.,;-=[]|        qrstuvwx12345678 

    __m256i indices = _mm256_setr_epi8(
        10, 5, 12, -1, 7, 8, 15, 6, 3, 9, 0, 14, 4, 1, 2, 11,
        4, 14, 7, 15, -1, 8, 9, 10, 1, 2, 0, 5, -1, 11, 6, 3
    );
    //";V=0" "X.|W" "T,Q]" "URS-"  "u7x8" "0123" "rsqv" "04wt"
    __m256 Result = _mm256_shuffle_epi8(lowerHalves, indices);

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
    // NOTE(casey): One-input puzzles
    //
    
    Match("Puzzle0", Puzzle0(A), "ijkl" "mnop" "abcd" "efgh" "1234" "5678" "qrst" "uvwx");
    Match("Puzzle1", Puzzle1(A), "5678" "qrst" "ijkl" "uvwx" "1234" "abcd" "efgh" "mnop");
    Match("Puzzle2", Puzzle2(A), "keia" "dbcm" "fglj" "nohp" "sr5t" "uq8w" "x123" "4v67");
    Match("Puzzle3", Puzzle3(A), "3fch" "1267" "ab8d" "45eg" "uxlq" "vmin" "opwj" "ktrs");
    
    //
    // NOTE(casey): Two-input puzzles
    //
    
    Match("Puzzle4", Puzzle4(A, B), "ijkl" "mnop" "ABCD" "EFGH" "qrst" "uvwx" ".,;-" "=[]|");
    Match("Puzzle5", Puzzle5(A, B), "ijkl" "abcd" "MNOP" "EFGH" "1234" "qrst" "=[]|" "UVWX");
    Match("Puzzle6", Puzzle6(A, B), ";V=0" "X.|W" "T,Q]" "URS-" "u7x8" "0123" "rsqv" "04wt");
    
    return 0;
}