#include "biginteger.h"

int main() {
  BigInteger bi1, bi2;
  std::cin >> bi1 >> bi2;
  Rational num1(bi1), denom1(bi2);
  if (bi1 < bi2) {
    std::cout << "LESS";
  } else if (bi1 == bi2) {
    std::cout << "EQ";
  } else if (bi1 > bi2) {
    std::cout << "Higher";
  }
  return 0;
}