#include <string>
#include <omp.h>
#include "bigint.h"
#include "primeswing.h"

constexpr int LIMB_DIGITS = 14;

std::string pad14(ull x) {
	std::string t = std::to_string(x);
	return std::string(LIMB_DIGITS - t.size(), '0') + t;
}

int main(int argc, char *argv[]) {
	int n;
	if(argc != 2 || sscanf(argv[1], "%d", &n) != 1) {
		fprintf(stderr, "Usage: %s [n]", argv[0]);
		return 1;
	}
	SPEED_TICK(tp0);
	sieve(n);
	auto ans = factorial(n);
	SPEED_TICK(tp1);
	SPEED_ADD_TOTAL(prod, tp0, tp1);
	const auto &w = ans.w;

	size_t digits = (w.size() - 1) * LIMB_DIGITS + std::to_string(w.back()).size();

	int zeros = 0;
	size_t j = 0;
	while(j < w.size() && w[j] == 0) zeros += LIMB_DIGITS, j++;
	int low_zeros = 0;
	for(ull x = w[j]; x % 10 == 0; x /= 10) low_zeros++;
	zeros += low_zeros;

	long long digitsum = 0;
	for(ull x : w)
		while(x) digitsum += x % 10, x /= 10;

	std::string head;
	for(int i = (int)w.size() - 1; i >= 0 && head.size() < 50; i--)
		head += (i == (int)w.size() - 1) ? std::to_string(w[i]) : pad14(w[i]);
	head.resize(std::min<size_t>(head.size(), 50));

	int m = low_zeros + 50;
	std::string mstr;
	for(int k = 4; k >= 0; k--)
		mstr += pad14(j + k < w.size() ? w[j + k] : 0);
	std::string tail = mstr.substr(mstr.size() - m, m);
	tail.resize(m - low_zeros);
	size_t trimmed = digits - zeros;
	if(trimmed <= 50) {
		size_t p = tail.find_first_not_of('0');
		tail = p == std::string::npos ? std::string("0") : tail.substr(p);
	}

	printf("%-9s%zu\n", "DIGITS", digits);
	printf("%-9s%d\n", "ZEROS", zeros);
	printf("%-9s%lld\n", "DIGITSUM", digitsum);
	printf("%-9s%s\n", "HEAD50", head.c_str());
	printf("%-9s%s\n", "TAIL50", tail.c_str());
	SPEED_TICK(tp2);
	SPEED_ADD_TOTAL(output, tp1, tp2);
	return 0;
}