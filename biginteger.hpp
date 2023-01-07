#pragma once
#include <cctype>
#include <cmath>
#include <compare>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

enum Ordering { Greater = 1, Less = -1, Equal = 0 };

class BigInteger {
 public:
  BigInteger() = default;
  BigInteger(int n) : is_positive_(n >= 0), digits_(1, std::abs(n)) {}

  explicit operator bool() const { return !isZero(); }

  [[nodiscard]] bool isPositive() const { return is_positive_; }
  void setSign(bool is_positive) { is_positive_ = isZero() || is_positive; }
  BigInteger operator-() const;

  BigInteger& operator+=(const BigInteger& added);
  BigInteger& operator-=(const BigInteger& subtrahend);
  BigInteger& operator*=(const BigInteger& multiplier);
  BigInteger& operator/=(const BigInteger& divider);
  BigInteger& operator%=(const BigInteger& mod);

  BigInteger& operator--();
  BigInteger operator--(int);
  BigInteger& operator++();
  BigInteger operator++(int);

  std::strong_ordering operator<=>(const BigInteger& bigint) const;

  [[nodiscard]] std::vector<int>& digits() { return digits_; }

  [[nodiscard]] std::string toString() const;

  [[nodiscard]] static int radix() { return kRadix; }
  [[nodiscard]] static int base() { return kBase; }

  static void swap(BigInteger& l_bi, BigInteger& r_bi);

  ~BigInteger() = default;

 private:
  static const int kRadix{1000000000};
  static const int kBase{10};
  static const int kPow{9};
  bool is_positive_{true};
  std::vector<int> digits_;

  [[nodiscard]] const std::vector<int>& digits() const { return digits_; }
  [[nodiscard]] size_t len() const { return digits_.size(); }
  [[nodiscard]] bool isZero() const {
    return digits_.empty() || (len() == 1 && digits_[0] == 0);
  }

  void absAdd(const BigInteger& bigint);
  void absSub(const BigInteger& subtrahend);
  void absMul(const BigInteger& multiplier);
  void absDiv(const BigInteger& divider);

  static int findEnoughMultiplier(const BigInteger& incomplete_divisible,
                                  const BigInteger& divider);
  static int absCmp(const BigInteger& l_bi, const BigInteger& r_bi);

  void removeLeadingZeros();
};

BigInteger operator ""_bi(unsigned long long bigint) {
  BigInteger num;
  for (long long radix(BigInteger::radix()); bigint > 0; bigint /= radix) {
    num.digits().emplace_back(bigint % radix);
  }
  return num;
}

void BigInteger::swap(BigInteger& l_bi, BigInteger& r_bi) {
  std::swap(l_bi.digits_, r_bi.digits_);
  std::swap(l_bi.is_positive_, r_bi.is_positive_);
}

BigInteger BigInteger::operator-() const {
  BigInteger negative_bi(*this);
  negative_bi.setSign(!isPositive());
  return negative_bi;
}

std::istream& operator>>(std::istream& bi_in, BigInteger& bigint) {
  bigint = 0_bi;
  char digit;
  bi_in.get(digit);
  while (!bi_in.eof() && std::isspace(digit) != 0) {
    bi_in.get(digit);
  }
  char is_positive('+');
  if (digit == '-') {
    is_positive = digit;
    bi_in.get(digit);
  }
  while (!bi_in.eof() && std::isspace(digit) == 0) {
    bigint *= BigInteger::base();
    bigint += static_cast<int>(digit - '0');
    bi_in.get(digit);
  }
  bigint.setSign(is_positive == '+');
  return bi_in;
}
std::ostream& operator<<(std::ostream& bi_out, const BigInteger& bigint) {
  return bi_out << bigint.toString();
}

BigInteger& BigInteger::operator--() {
  *this -= 1_bi;
  return *this;
}
BigInteger BigInteger::operator--(int) {
  BigInteger copy(*this);
  operator--();
  return copy;
}
BigInteger& BigInteger::operator++() {
  *this += 1_bi;
  return *this;
}
BigInteger BigInteger::operator++(int) {
  BigInteger copy(*this);
  operator++();
  return copy;
}

std::strong_ordering BigInteger::operator<=>(const BigInteger& bigint) const {
  switch (absCmp(*this, bigint)) {
    case Ordering::Greater:
      return isPositive() ? std::strong_ordering::greater
                          : std::strong_ordering::less;
    case Ordering::Less:
      if (isPositive() == bigint.isPositive()) {
        return isPositive() ? std::strong_ordering::less
                            : std::strong_ordering::greater;
      }
      return isPositive() ? std::strong_ordering::greater
                          : std::strong_ordering::less;
  }
  if (isPositive() != bigint.isPositive()) {
    return isPositive() ? std::strong_ordering::greater
                        : std::strong_ordering::less;
  }
  return std::strong_ordering::equal;
}
bool operator==(const BigInteger& l_bi, const BigInteger& r_bi) {
  return std::is_eq(l_bi <=> r_bi);
}

BigInteger operator+(const BigInteger& l_bi, const BigInteger& r_bi) {
  BigInteger res(l_bi);
  res += r_bi;
  return res;
}
BigInteger operator-(const BigInteger& l_bi, const BigInteger& r_bi) {
  BigInteger res(l_bi);
  res -= r_bi;
  return res;
}
BigInteger operator*(const BigInteger& l_bi, const BigInteger& r_bi) {
  BigInteger res(l_bi);
  res *= r_bi;
  return res;
}
BigInteger operator/(const BigInteger& l_bi, const BigInteger& r_bi) {
  BigInteger res(l_bi);
  res /= r_bi;
  return res;
}
BigInteger operator%(const BigInteger& l_bi, const BigInteger& r_bi) {
  BigInteger res(l_bi);
  res %= r_bi;
  return res;
}

BigInteger& BigInteger::operator+=(const BigInteger& added) {
  if (isPositive() == added.isPositive()) {
    absAdd(added);
  } else {
    switch (absCmp(*this, added)) {
      case Ordering::Equal:
        *this = 0_bi;
        break;
      case Ordering::Greater:
        absSub(added);
        break;
      case Ordering::Less:
        BigInteger tmp(added);
        std::swap(*this, tmp);
        absSub(tmp);
        break;
    }
  }
  removeLeadingZeros();
  return *this;
}
BigInteger& BigInteger::operator-=(const BigInteger& subtrahend) {
  setSign(!isPositive());
  *this += subtrahend;
  setSign(!isPositive());
  return *this;
}
BigInteger& BigInteger::operator*=(const BigInteger& multiplier) {
  bool new_is_positive(isPositive() == multiplier.isPositive());
  absMul(multiplier);
  setSign(new_is_positive);
  removeLeadingZeros();
  return *this;
}
BigInteger& BigInteger::operator/=(const BigInteger& divider) {
  bool new_is_positive(isPositive() == divider.isPositive());
  absDiv(divider);
  setSign(new_is_positive);
  removeLeadingZeros();
  return *this;
}
BigInteger& BigInteger::operator%=(const BigInteger& mod) {
  *this -= (*this / mod) * mod;
  return *this;
}

void BigInteger::absAdd(const BigInteger& bigint) {
  int overflow(0);
  for (size_t i(0); i < std::max(len(), bigint.len()) || overflow > 0; ++i) {
    if (i == digits_.size()) {
      digits_.emplace_back(0);
    }
    digits_[i] += (i < bigint.len() ? bigint.digits()[i] : 0) + overflow;
    overflow = static_cast<int>(digits_[i] >= radix());
    digits_[i] -= overflow > 0 ? radix() : 0;
  }
}
void BigInteger::absSub(const BigInteger& subtrahend) {
  int underflow = 0;
  for (size_t i(0); i < subtrahend.len() || underflow > 0; ++i) {
    digits_[i] -=
        underflow + (i < subtrahend.len() ? subtrahend.digits_[i] : 0);
    underflow = static_cast<int>(digits_[i] < 0);
    digits_[i] += underflow > 0 ? radix() : 0;
  }
}
void BigInteger::absMul(const BigInteger& multiplier) {
  BigInteger composition;
  composition.digits_.resize(len() + multiplier.len() - 1);
  long long overflow = 0;
  for (size_t rank(0); rank < composition.len(); ++rank) {
    long long column = overflow % radix();
    overflow /= radix();
    for (size_t i(std::max(static_cast<int>(rank - multiplier.len()) + 1, 0));
         i < std::min(rank + 1, len()); ++i) {
      column +=
          static_cast<long long>(digits_[i]) * multiplier.digits_[rank - i];
      overflow += column / radix();
      column %= radix();
    }
    composition.digits_[rank] = static_cast<int>(column);
  }
  while (overflow > 0) {
    composition.digits_.emplace_back(overflow % radix());
    overflow /= radix();
  }
  swap(*this, composition);
}
void BigInteger::absDiv(const BigInteger& divider) {
  BigInteger quotient;
  BigInteger incomplete_divisible;
  for (long long i(static_cast<long long>(len()) - 1LL); i >= 0; --i) {
    incomplete_divisible *= radix();
    quotient *= radix();
    incomplete_divisible += digits_[i];
    if (absCmp(incomplete_divisible, divider) < 0) {
      continue;
    }
    int enough_mul = findEnoughMultiplier(incomplete_divisible, divider);
    quotient += enough_mul;
    incomplete_divisible -=
        (divider.isPositive() ? divider : -divider) * enough_mul;
  }
  swap(*this, quotient);
}

std::string BigInteger::toString() const {
  if (isZero()) {
    return "0";
  }
  std::stringstream bi_stream;
  if (!isPositive()) {
    bi_stream << '-';
  }
  bi_stream << digits_.back();
  char fill_null = bi_stream.fill('0');
  for (long long i(static_cast<long long>(len()) - 2LL); i >= 0; --i) {
    bi_stream << std::setw(kPow) << digits_[i];
  }
  bi_stream.fill(fill_null);
  return bi_stream.str();
}
void BigInteger::removeLeadingZeros() {
  while (len() > 1 && digits_.back() == 0) {
    digits_.pop_back();
  }
  if (len() == 1 && digits_[0] == 0) {
    is_positive_ = true;
  }
}
int BigInteger::absCmp(const BigInteger& l_bi, const BigInteger& r_bi) {
  if (l_bi.isZero() && r_bi.isZero()) {
    return Ordering::Equal;
  }
  if (l_bi.len() < r_bi.len()) {
    return Ordering::Less;
  }
  if (l_bi.len() > r_bi.len()) {
    return Ordering::Greater;
  }
  for (size_t i(l_bi.len() - 1); i != size_t(-1); --i) {
    if (l_bi.digits_[i] < r_bi.digits_[i]) {
      return Ordering::Less;
    }
    if (l_bi.digits_[i] > r_bi.digits_[i]) {
      return Ordering::Greater;
    }
  }
  return Ordering::Equal;
}
int BigInteger::findEnoughMultiplier(const BigInteger& incomplete_divisible,
                                     const BigInteger& divider) {
  int lowest = 1;
  int biggest = radix();
  while (biggest - lowest > 1) {
    int median = (biggest + lowest) / 2;
    BigInteger res = divider * static_cast<BigInteger>(median);
    if (absCmp(res, incomplete_divisible) <= Ordering::Equal) {
      lowest = median;
    } else {
      biggest = median;
    }
  }
  return lowest;
}

BigInteger GCD(BigInteger l_bi, BigInteger r_bi) {
  while (r_bi) {
    l_bi %= r_bi;
    BigInteger::swap(l_bi, r_bi);
  }
  return l_bi;
}

class Rational {
 public:
  Rational() = default;
  explicit Rational(int num) : numerator_(num) {}
  Rational(const BigInteger& bigint) : numerator_(bigint) {}

  explicit operator double() const { return std::stod(toString()); }

  [[nodiscard]] bool isPositive() const { return numerator_.isPositive(); }
  Rational operator-() const;

  [[nodiscard]] std::pair<BigInteger, BigInteger> fraction() const {
    return {numerator_, denominator_};
  }

  std::strong_ordering operator<=>(const Rational& rat) const;

  Rational& operator+=(const Rational& fraction);
  Rational& operator-=(const Rational& subtrahend);
  Rational& operator*=(const Rational& multiplier);
  Rational& operator/=(const Rational& divider);

  [[nodiscard]] std::string asDecimal(size_t precision) const;
  [[nodiscard]] std::string toString() const;

  ~Rational() = default;

 private:
  BigInteger numerator_{0};
  BigInteger denominator_{1};

  void shortenFraction();
  void add(const Rational& fraction);

  void setSign(bool is_positive) { numerator_.setSign(is_positive); }
};

Rational Rational::operator-() const {
  Rational negative_fraction(*this);
  negative_fraction.setSign(!isPositive());
  return negative_fraction;
}
void Rational::shortenFraction() {
  setSign(numerator_.isPositive() == denominator_.isPositive());
  denominator_.setSign(true);
  BigInteger gcd = GCD(numerator_, denominator_);
  gcd.setSign(true);
  numerator_ /= gcd;
  denominator_ /= gcd;
}

std::strong_ordering Rational::operator<=>(const Rational& rat) const {
  BigInteger l_rat_num = numerator_ * rat.denominator_;
  BigInteger r_rat_num = rat.numerator_ * denominator_;
  return l_rat_num < r_rat_num ? std::strong_ordering::less
                               : std::strong_ordering::greater;
}
bool operator==(const Rational& l_rat, const Rational& r_rat) {
  return l_rat.fraction() == r_rat.fraction();
}

Rational& Rational::operator+=(const Rational& fraction) {
  add(fraction);
  shortenFraction();
  return *this;
}
Rational& Rational::operator-=(const Rational& subtrahend) {
  setSign(!isPositive());
  *this += subtrahend;
  setSign(!isPositive());
  return *this;
}
Rational& Rational::operator*=(const Rational& multiplier) {
  numerator_ *= multiplier.numerator_;
  denominator_ *= multiplier.denominator_;
  shortenFraction();
  return *this;
}
Rational& Rational::operator/=(const Rational& divider) {
  numerator_ *= divider.denominator_;
  denominator_ *= divider.numerator_;
  shortenFraction();
  return *this;
}

Rational operator+(const Rational& l_rat, const Rational& r_rat) {
  Rational res(l_rat);
  res += r_rat;
  return res;
}
Rational operator-(const Rational& l_rat, const Rational& r_rat) {
  Rational res(l_rat);
  res -= r_rat;
  return res;
}
Rational operator*(const Rational& l_rat, const Rational& r_rat) {
  Rational res(l_rat);
  res *= r_rat;
  return res;
}
Rational operator/(const Rational& l_rat, const Rational& r_rat) {
  Rational res(l_rat);
  res /= r_rat;
  return res;
}

std::string Rational::asDecimal(size_t precision) const {
  size_t boost(precision);
  BigInteger num(isPositive() ? numerator_ : -numerator_);
  BigInteger precision_boost(1);
  while (boost-- > 0) {
    precision_boost *= BigInteger::base();
  }
  num *= precision_boost;
  num /= denominator_;
  BigInteger leftover(num % precision_boost);
  num /= precision_boost;
  if (precision == 0) {
    return num.toString();
  }
  std::string dec_str(!isPositive() ? 1 : 0, '-');
  dec_str += num.toString();
  dec_str.push_back('.');
  dec_str += std::string(precision - leftover.toString().length(), '0');
  dec_str += leftover.toString();
  return dec_str;
}
std::string Rational::toString() const {
  std::string fraction(numerator_.toString());
  if (numerator_ % denominator_ != 0_bi) {
    fraction.push_back('/');
    fraction += denominator_.toString();
  }
  return fraction;
}

void Rational::add(const Rational& fraction) {
  numerator_ *= fraction.denominator_;
  numerator_ += fraction.numerator_ * denominator_;
  denominator_ *= fraction.denominator_;
}
