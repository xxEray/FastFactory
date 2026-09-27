#include <string>
#include <omp.h>
#include "bigint.h"

constexpr int STUDENT_ID = 5;

constexpr int pow10[] = {1, 10, 100, 1000, 10000, 100000, 1000000, 10000000, 100000000};

constexpr int LIMB_DIGITS = 14;

std::string pad14(ull x) {
	std::string t = std::to_string(x);
	return std::string(LIMB_DIGITS - t.size(), '0') + t;
}

BigInt get_prod(int l, int r) {
	if(r - l + 1 <= LEAF_PIVOT) {
		BigInt ret(l);
		for(int i = l + 1; i <= r; i++) ret *= static_cast<unsigned>(i);
		return ret;
	}
	int mid = (l + r) / 2;
	BigInt val = get_prod(l, mid);
	mul_eq(val, get_prod(mid + 1, r));
	return val;
}

int main(int argc, char *argv[]) {
	if(argc != 2 || argv[1][0] < '4' || argv[1][0] > '8' || argv[1][1] != 0) {
		fprintf(stderr, "Usage: %s [k]", argv[0]);
		return 1;
	}
	int n = pow10[argv[1][0] - '0'] + STUDENT_ID * 1000;
	SPEED_TICK(tp0);
	auto ans = get_prod(1, n);
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