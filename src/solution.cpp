// Powered by AI
// 获取正确答案 & 得到运行时间参考

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <gmp.h>

int main(int argc, char **argv) {
    if (argc < 2) {
        fprintf(stderr, "usage: %s N\n", argv[0]);
        return 1;
    }
    freopen("out.txt", "w", stdout);

    unsigned long n = std::strtoul(argv[1], nullptr, 10);

    printf("n %lu\n", n);

    mpz_t f;
    mpz_init(f);
    mpz_fac_ui(f, n);

    char *s = mpz_get_str(nullptr, 10, f);
    std::size_t len = std::strlen(s);

    printf("DIGITS   %zu\n", len);

    std::size_t tz = 0;
    while (tz < len && s[len - 1 - tz] == '0')
        tz++;

    printf("ZEROS    %zu\n", tz);

    unsigned long long sum = 0;
    for (std::size_t i = 0; i < len; i++)
        sum += static_cast<unsigned long long>(s[i] - '0');
    printf("DIGITSUM %llu\n", sum);

    printf("HEAD50   ");
    for (std::size_t i = 0; i < 50 && i < len; i++)
        std::putchar(s[i]);
    std::putchar('\n');

    std::size_t core = len - tz;
    printf("TAIL50   ");
    if (core <= 50) {
        std::fwrite(s, 1, core, stdout);
    } else {
        std::fwrite(s + core - 50, 1, 50, stdout);
    }
    std::putchar('\n');

    std::free(s);
    mpz_clear(f);
    return 0;
}