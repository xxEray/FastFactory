#include <string>
#include "bigint.h"

constexpr int STUDENT_ID = 5;

constexpr int pow10[] = {1, 10, 100, 1000, 10000, 100000, 1000000, 10000000, 100000000};

BigInt get_prod(int l, int r) {
	if(l == r) return BigInt(l);
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
	auto ans = get_prod(1, n);
	std::string s = std::to_string(ans.w.back());
	for(int i = (int)ans.w.size() - 2; i >= 0; i--) {
		std::string t = std::to_string(ans.w[i]);
		s += std::string(14 - t.size(), '0') + t;
	}
	int zeros = 0;
	for(int i = (int)s.size() - 1; i >= 0 && s[i] == '0'; i--) zeros++;
	long long digitsum = 0;
	for(char c : s) digitsum += c - '0';
	std::string trimmed = s.substr(0, s.size() - zeros);
	std::string head = s.substr(0, 50);
	std::string tail = trimmed.substr(trimmed.size() > 50 ? trimmed.size() - 50 : 0);
	printf("%-9s%zu\n", "DIGITS", s.size());
	printf("%-9s%d\n", "ZEROS", zeros);
	printf("%-9s%lld\n", "DIGITSUM", digitsum);
	printf("%-9s%s\n", "HEAD50", head.c_str());
	printf("%-9s%s\n", "TAIL50", tail.c_str());
	return 0;
}