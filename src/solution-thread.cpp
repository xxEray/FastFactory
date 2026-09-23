#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <future>
#include <thread>
#include <string>
#include <gmpxx.h>

static mpz_class prod_seq(unsigned long lo, unsigned long hi) {
    if (lo > hi)
        return mpz_class(1);
    if (hi - lo < 64) {
        mpz_class r = 1;
        for (unsigned long i = lo; i <= hi; ++i)
            r *= i;
        return r;
    }
    unsigned long mid = lo + (hi - lo) / 2;
    return prod_seq(lo, mid) * prod_seq(mid + 1, hi);
}

static mpz_class prod_par(unsigned long lo, unsigned long hi, int depth) {
    if (lo > hi)
        return mpz_class(1);
    if (depth <= 0 || hi - lo < 64)
        return prod_seq(lo, hi);
    unsigned long mid = lo + (hi - lo) / 2;
    std::future<mpz_class> left =
        std::async(std::launch::async, prod_par, lo, mid, depth - 1);
    mpz_class right = prod_par(mid + 1, hi, depth - 1);
    return left.get() * right;
}

int main(int argc, char **argv) {
    if (argc < 2) {
        std::fprintf(stderr, "usage: %s N\n", argv[0]);
        return 1;
    }

    unsigned long n = std::strtoul(argv[1], nullptr, 10);

    unsigned threads = std::thread::hardware_concurrency();
    if (threads == 0)
        threads = 1;
    int depth = 0;
    while (depth < 20 && (1u << (depth + 1)) <= threads)
        depth++;

    mpz_class f = (n < 2) ? mpz_class(1) : prod_par(1, n, depth);

    std::string str = f.get_str();
    const char *s = str.c_str();
    std::size_t len = str.size();

    std::printf("DIGITS %zu\n", len);

    std::size_t tz = 0;
    while (tz < len && s[len - 1 - tz] == '0')
        tz++;

    std::printf("ZEROS\n%zu\n", tz);

    unsigned long long sum = 0;
    for (std::size_t i = 0; i < len; i++)
        sum += static_cast<unsigned long long>(s[i] - '0');
    std::printf("DIGITSUM %llu\n", sum);

    std::printf("HEAD50 ");
    for (std::size_t i = 0; i < 50 && i < len; i++)
        std::putchar(s[i]);
    std::putchar('\n');

    std::size_t core = len - tz;
    std::printf("TAIL50 ");
    if (core <= 50)
        std::fwrite(s, 1, core, stdout);
    else
        std::fwrite(s + core - 50, 1, 50, stdout);
    std::putchar('\n');

    return 0;
}
