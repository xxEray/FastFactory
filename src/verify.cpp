#include <cstdio>
#include <cstdint>
#include <string>
#include <array>
#include <algorithm>

#include "bigint.h"

using ull = unsigned long long;
using u128 = __uint128_t;

constexpr int STUDENT_ID = 5;

constexpr int pow10[] = {
    1,
    10,
    100,
    1000,
    10000,
    100000,
    1000000,
    10000000,
    100000000
};

constexpr int LIMB_DIGITS = 14;
constexpr ull BASE = 100000000000000ULL; // 10^14

// 多个彼此独立的验证模数。
// 它们都大于本测试中的 n，因此不会出现 n! mod p = 0 的情况。
constexpr ull VERIFY_PRIMES[] = {
    1000000007ULL,
    1000000009ULL,
    998244353ULL,
    100005001ULL,
    100005023ULL
};

constexpr int PRIME_COUNT =
    sizeof(VERIFY_PRIMES) / sizeof(VERIFY_PRIMES[0]);


// ------------------------------------------------------------
// Product Tree
// ------------------------------------------------------------

BigInt get_prod(int l, int r) {
    if (r - l + 1 <= LEAF_PIVOT) {
        BigInt ret(l);
        for (int i = l + 1; i <= r; i++)
            ret *= static_cast<unsigned>(i);
        return ret;
    }

    int mid = (l + r) / 2;

    BigInt val = get_prod(l, mid);
    mul_eq(val, get_prod(mid + 1, r));

    return val;
}


// ------------------------------------------------------------
// 从 BigInt 独立计算 x mod p
//
// BigInt:
//     x = w[0] + w[1] * BASE + w[2] * BASE^2 + ...
//
// 注意：这里只使用 BigInt 的公开数据 w，
// 不使用 NTT / Montgomery / CRT 的任何中间结果。
// ------------------------------------------------------------

ull bigint_mod(const BigInt &x, ull p) {
    ull r = 0;

    for (auto it = x.w.rbegin(); it != x.w.rend(); ++it) {
        r = static_cast<ull>(
            (static_cast<u128>(r) * (BASE % p) + (*it % p)) % p
        );
    }

    return r;
}


// ------------------------------------------------------------
// 独立计算 n! mod p
//
// 故意不使用 BigInt、Product Tree、NTT。
// 就是最简单的：
//     1 * 2 * 3 * ... * n mod p
//
// 这样可以作为最终大整数结果的独立验证路径。
// ------------------------------------------------------------

ull factorial_mod(int n, ull p) {
    ull ret = 1;

    for (int i = 2; i <= n; ++i) {
        ret = static_cast<ull>(
            static_cast<u128>(ret) * static_cast<ull>(i) % p
        );
    }

    return ret;
}


// ------------------------------------------------------------
// Legendre:
//     v_p(n!) = floor(n/p) + floor(n/p^2) + ...
// ------------------------------------------------------------

long long legendre(int n, int p) {
    long long result = 0;

    while (n > 0) {
        n /= p;
        result += n;
    }

    return result;
}


// ------------------------------------------------------------
// 从最终 BigInt 独立计算十进制末尾 0
//
// 因为 BASE = 10^14，
// 每个完整的零 limb 对应 14 个十进制 0。
// ------------------------------------------------------------

long long bigint_trailing_zeros(const BigInt &x) {
    const auto &w = x.w;

    if (w.empty())
        return 0;

    size_t j = 0;
    long long zeros = 0;

    // 完整的 10^14 limb
    while (j < w.size() && w[j] == 0) {
        zeros += LIMB_DIGITS;
        ++j;
    }

    // 最低的非零 limb 中继续统计 10 的因子
    if (j < w.size()) {
        ull x0 = w[j];

        while (x0 != 0 && x0 % 10 == 0) {
            ++zeros;
            x0 /= 10;
        }
    }

    return zeros;
}


// ------------------------------------------------------------
// 主验证
// ------------------------------------------------------------

int main(int argc, char *argv[]) {

    if (argc != 2 ||
        argv[1][0] < '4' ||
        argv[1][0] > '8' ||
        argv[1][1] != 0) {

        fprintf(stderr, "Usage: %s [k]\n", argv[0]);
        return 1;
    }

    int k = argv[1][0] - '0';
    int n = pow10[k] + STUDENT_ID * 1000;

    printf("========================================\n");
    printf("Verification\n");
    printf("========================================\n");
    printf("k = %d\n", k);
    printf("n = %d\n", n);
    printf("\n");

    // --------------------------------------------------------
    // 1. 计算 n!
    // --------------------------------------------------------

    printf("[1] Computing n! ...\n");

    BigInt ans = get_prod(1, n);

    printf("    computation finished.\n");
    printf("\n");


    // --------------------------------------------------------
    // 2. 验证末尾 0
    // --------------------------------------------------------

    printf("[2] Verifying trailing zeros\n");

    // 从最终 BigInt 直接统计
    long long zeros_bigint = bigint_trailing_zeros(ans);

    // 数学公式：
    // v_10(n!) = min(v_2(n!), v_5(n!))
    // 对 n! 来说 v_2 > v_5，因此：
    // v_10(n!) = v_5(n!)
    long long zeros_legendre = legendre(n, 5);

    printf("    BigInt trailing zeros : %lld\n", zeros_bigint);
    printf("    Legendre v5(n!)       : %lld\n", zeros_legendre);

    bool zero_ok = (zeros_bigint == zeros_legendre);

    if (zero_ok) {
        printf("    RESULT                 : PASS\n");
    } else {
        printf("    RESULT                 : FAIL\n");
    }

    printf("\n");


    // --------------------------------------------------------
    // 3. 多模数验证
    // --------------------------------------------------------

    printf("[3] Multi-mod verification\n");
    printf("\n");

    bool mod_ok = true;

    for (int i = 0; i < PRIME_COUNT; ++i) {

        ull p = VERIFY_PRIMES[i];

        printf("    p = %llu\n",
               static_cast<unsigned long long>(p));

        // 从最终 BigInt 计算：
        ull result_from_bigint = bigint_mod(ans, p);

        // 完全独立地：
        // 1 * 2 * ... * n mod p
        ull result_independent = factorial_mod(n, p);

        printf("        BigInt mod p : %llu\n",
               static_cast<unsigned long long>(result_from_bigint));

        printf("        n! mod p     : %llu\n",
               static_cast<unsigned long long>(result_independent));

        bool ok = (result_from_bigint == result_independent);

        if (ok) {
            printf("        RESULT        : PASS\n");
        } else {
            printf("        RESULT        : FAIL\n");
            mod_ok = false;
        }

        printf("\n");
    }


    // --------------------------------------------------------
    // 4. 最终结果
    // --------------------------------------------------------

    bool all_ok = zero_ok && mod_ok;

    printf("========================================\n");

    if (all_ok) {
        printf("VERIFY PASS\n");
    } else {
        printf("VERIFY FAIL\n");
    }

    printf("========================================\n");

    return all_ok ? 0 : 1;
}