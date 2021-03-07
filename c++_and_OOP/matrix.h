#include<iostream>
#include<vector>
#include<string>
#include <math.h>
#include <utility>
#include <cassert>

//I prefer using snake case
#define residue Residue
#define get_inverse getInverse
#define get_primitive_root getPrimitiveRoot
#define to_string toString 
#define big_integer BigInteger
#define rational Rational 
#define as_decimal asDecimal 

#define matrix Matrix
#define get_row getRow
#define get_column getColumn
#define square_matrix SquareMatrix
 
class big_integer {
private:
    static int const BASE = 10;
    static int const HALF_BASE = 5;
    static int const mod = 998244353;
    static const int root = 15311432;
    static const int root_inv = 469870224;
    static const int root_pw = 1 << 23;       
    std::vector<int> num = std::vector<int>(1, 0);
    bool neg = 0;
    void leading_zeros();
    void shift();
    void format();
    static size_t min(size_t a, size_t b) { return a < b ? a : b; }
    static size_t max(size_t a, size_t b) { return a > b ? a : b; }
    static int mod_mlt(int a, int b) { return (a * 1ll * b) % mod; }
    static int mod_inv(int x) { return mod_bin_pow(x, mod - 2); }
    static int mod_sum(int a, int b);
    static int mod_sub(int a, int b);
    static int mod_bin_pow(int a, int n); 
    static void ntt(std::vector<int>& a, bool invert);
    static void multiply(std::vector<int>& a, std::vector<int> b);
    void addition(const big_integer& x);
    void subtraction(const big_integer& x);
    void summation(const big_integer& x, bool sub);
    void division(const big_integer& x, bool mod);
    void ban_negative_zero();
    void divide_by_two();
public:
    big_integer() = default;
    big_integer(const std::string& str) { 
        if (!str.size())
            return;
        neg = (str[0] == '-');
        num.resize(str.size() - neg);
        auto num_it = num.end() - 1; 
        for (auto str_it = str.begin() + neg; str_it != str.end(); ++str_it, --num_it) 
            if (*str_it == '-') {
                ++num_it;
                neg = 1;
            } else {
                *num_it = (*str_it - '0');
            }
        ban_negative_zero();
    }
    big_integer(int x) {
        if (x == 0) 
            return;
        if (x < 0) {
            neg = 1;
            x = -x;
        }
        num.clear();
        while (x) {
            num.push_back(x % BASE);
            x /= BASE;
        }
    }
    big_integer(const char* s) : big_integer(std::string(s)) {}
    big_integer(const big_integer& x) : num(std::vector<int>(x.num.begin(), x.num.end())), neg(x.neg) {}
    ~big_integer() = default;
    big_integer& operator+=(const big_integer& x) &;
    big_integer& operator-=(const big_integer& x) &;
    big_integer& operator*=(const big_integer& x) &;
    big_integer& operator/=(const big_integer& x) &;
    big_integer& operator%=(const big_integer& x) &;
    big_integer& operator=(big_integer x) &;
    big_integer& operator++() &;
    big_integer operator++(int) &;
    big_integer& operator--() &;
    big_integer operator--(int) &;
    big_integer operator-() const;
    bool negative() const;
    bool is_zero() const;
    bool is_one() const;
    void change_sign();
    bool less_abs(const big_integer& x) const;
    int digit(size_t index) const;
    int& digit(size_t index);
    size_t length() const;
    std::string to_string() const;
    void swap(big_integer& b) &;
    explicit operator bool() const { return !is_zero(); }      
    big_integer multiply_by_power_of_ten(size_t ppw) const;
    explicit operator int() const {
        int ans = 0;
        for (size_t pos = length(); pos;) {
            ans = ans * BASE + num[--pos];
        }
        return ans * (neg ? -1 : 1);
    }
};

int big_integer::mod_sum(int a, int b) {
    a += b;
    return a < mod ? a : a - mod;
}

int big_integer::mod_sub(int a, int b) {
        a -= b;
    return a < 0 ? a + mod : a;
}

int big_integer::mod_bin_pow(int a, int n) {
    int res = 1;
    while (n) {
        if(n & 1) {
            res = mod_mlt(a, res), n ^= 1;
        } else {
            a = mod_mlt(a, a), n >>= 1;
        }
    }
    return res;    
}

//Number theoretic transform (I thought that's well-known name)
//time complexity O(N * log(N))  
void big_integer::ntt(std::vector<int>& a, bool invert) {
    size_t n = a.size();
    for (size_t i = 1, j = 0; i < n; ++i) {
        size_t bit = n >> 1;
        for (; j >= bit; bit >>= 1) {
            j -= bit;
        }
        j += bit;
        if (i < j)
            std::swap(a[i], a[j]);
    } 
    for (size_t hlen = 1, len = 2; len <= n; hlen = len, len <<= 1) {
        int wlen = invert ? root_inv : root;
        for (int i = len; i < root_pw; i <<= 1) {
            wlen = mod_mlt(wlen, wlen);
        }
        auto j = a.begin() + hlen;
        for (auto i = a.begin(); i != a.end(); i += hlen, j += hlen) {
            auto ni = i + len;
            int w = 1;
            for (; j != ni; ++i, ++j) {
                int j_val = mod_mlt(*j, w);
                *j = mod_sub(*i, j_val);
                *i = mod_sum(*i, j_val);
                w = mod_mlt(w, wlen);    
            }
        }
    }
    if (invert) {
        int n_inv = mod_inv(n);
        for (size_t i = 0; i < n; ++i) {
            a[i] = mod_mlt(a[i], n_inv);
        }
    }
}

void big_integer::multiply(std::vector<int>& a, std::vector<int> b){
    size_t n = 1;
    while (n < std::max(a.size(), b.size())) {
        n <<= 1;
    }
    n <<= 1;
    a.resize(n, 0);
    b.resize(n, 0);
    ntt(a, 0);
    ntt(b, 0);
    for (size_t i = 0; i < n; ++i) {
        a[i] = mod_mlt(a[i], b[i]);
    }
    ntt(a, 1);    
}

big_integer big_integer::operator-() const {
    big_integer cp(*this);
    cp.change_sign();
    return cp;
}

bool big_integer::is_zero() const{
    return num.back() == 0;
}

bool big_integer::is_one() const {
    return length() == 1 && num[0] == 1;    
}

void big_integer::change_sign(){    
    if (!is_zero())
        neg = !neg;    
}

big_integer abs(big_integer q) {
    if (q.negative())
        q.change_sign();
    return q;
}

big_integer operator""_bi(const char* s) {
    return big_integer(s);
}

bool big_integer::less_abs(const big_integer& x) const {
    if (length() != x.length())
        return length() < x.length();
    for (size_t pos = length(); pos;) {
        --pos;
        if (num[pos] != x.num[pos])
            return num[pos] < x.num[pos];
    }
    return 0;        
}

big_integer& big_integer::operator+=(const big_integer& x) & {
    summation(x, negative() != x.negative());
    return *this;
}

big_integer& big_integer::operator-=(const big_integer& x) & {
    summation(x, negative() == x.negative());
    return *this;    
}

big_integer& big_integer::operator*=(const big_integer& x) & {
    neg ^= x.neg;
    multiply(num, x.num);
    format();
    return *this;
}

big_integer& big_integer::operator/=(const big_integer& x) & {
    division(x, 0);
    return *this;
}

big_integer& big_integer::operator%=(const big_integer& x) & {
    division(x, 1);
    return *this;    
}

big_integer& big_integer::operator=(big_integer x) & {
    swap(x);
    return *this;
}
    
void big_integer::swap(big_integer& x) & {
    std::swap(neg, x.neg);
    std::swap(num, x.num);
}

big_integer& big_integer::operator++() & {
    *this += 1_bi;
    return *this;
}

big_integer& big_integer::operator--() & {
    *this -= 1_bi;
    return *this;
}

big_integer big_integer::operator++(int) & {
    big_integer res(*this);
    *this += 1_bi;
    return res;
}

big_integer big_integer::operator--(int) & {
    big_integer res(*this);
    *this -= 1_bi;
    return res;
}

big_integer operator+(big_integer lhs, const big_integer& rhs) {
    return lhs += rhs;    
}

big_integer operator-(big_integer lhs, const big_integer& rhs) {
    return lhs -= rhs;    
}

big_integer operator*(big_integer lhs, const big_integer& rhs) {
    return lhs *= rhs;
}

big_integer operator/(big_integer lhs, const big_integer& rhs) {
    return lhs /= rhs;
}

big_integer operator%(big_integer lhs, const big_integer& rhs) {
    return lhs %= rhs;
}

bool operator<(const big_integer& lhs, const big_integer& rhs) {
    if (lhs.negative() != rhs.negative())
        return lhs.negative();
    if (lhs.less_abs(rhs))
        return !lhs.negative();
    if (rhs.less_abs(lhs))
        return lhs.negative();
    return 0;     
}

bool operator>(const big_integer& lhs, const big_integer& rhs) {
    return rhs < lhs;
}

bool operator<=(const big_integer& lhs, const big_integer& rhs) {
    return !(rhs < lhs);
}

bool operator>=(const big_integer& lhs, const big_integer& rhs) {
    return !(lhs < rhs);
}

bool operator==(const big_integer& lhs, const big_integer& rhs) {
    return !(lhs < rhs || rhs < lhs);
}

bool operator!=(const big_integer& lhs, const big_integer& rhs) {
    return lhs < rhs || rhs < lhs;
}    

std::istream& operator>>(std::istream& in, big_integer& x) {
    std::string input;
    in >> input;
    x = big_integer(input);
    return in;
}

std::ostream& operator<<(std::ostream& out, const big_integer& x) {
    out << x.to_string();
    return out;
}

void big_integer::leading_zeros() {
    while (num.size() > 1 && !num.back()) {
        num.pop_back();    
    }
}

void big_integer::shift() {
    int bonus = 0;
    for (size_t pos = 0; pos < num.size(); ++pos) {
        num[pos] += bonus;
        bonus = num[pos] / BASE;
        num[pos] %= BASE;
    }
}

void big_integer::format() {
    shift();
    leading_zeros();
    ban_negative_zero(); 
}

void big_integer::ban_negative_zero() {
    if (!num.back())
        neg = 0;    
}

void big_integer::addition(const big_integer& x) {
    int bonus = 0;
    size_t nsz = max(length(), x.length());
    num.resize(nsz, 0);
    for (size_t pos = 0; pos < nsz; ++pos) {
        num[pos] = num[pos] + (pos < x.length() ? x.num[pos] : 0) + bonus;
        bonus = 0;
        while (num[pos] >= BASE) {
            num[pos] -= BASE;
            ++bonus;
        }        
    }    
    if (bonus) 
        num.push_back(bonus);
}

//less_abs works linear only if numbers have same length - so it's O(changed digits)
void big_integer::subtraction(const big_integer& x) {
    if (less_abs(x)) {  
        bool new_neg = !neg;
        big_integer cp(*this);
        *this = big_integer(x);
        subtraction(cp);
        neg = new_neg;
        return;                
    }
    int bonus = 0;
    for (size_t pos = 0; bonus || pos < x.length(); ++pos) {
        num[pos] = num[pos] - bonus - (pos < x.length() ? x.num[pos] : 0);
        bonus = 0;
        while (num[pos] < 0) {
            num[pos] += BASE;
            ++bonus;
        }
    }
    leading_zeros();
    ban_negative_zero();        
}

void big_integer::summation(const big_integer& x, bool sub) {
    if (sub)
        subtraction(x);
    else
        addition(x);
}

void big_integer::divide_by_two() {
    bool bonus = 0;
    for (size_t pos = length(); pos;) {
        --pos;
        bool nbonus = (num[pos] & 1);
        num[pos] >>= 1;
        if (bonus)
            num[pos] += HALF_BASE;
        bonus = nbonus;            
    }
    leading_zeros();
}

/*
restoring qoutient's binnary representation
O(N) times try to add a power of two to the answer
every time it takes O(N) time to recalculate all stuff
total coplexity: O(N^2)
*/
    
void big_integer::division(const big_integer& x, bool need_rem){
    big_integer quot(0); 
    bool rem_neg = neg;
    bool quot_neg = neg ^ x.neg;
    big_integer rem = abs(*this);
    size_t p = 0;
    big_integer ppw(1);  
    big_integer x_ppw = abs(x);
    while (x_ppw <= rem) {
        x_ppw += x_ppw;
        ppw += ppw;
        ++p;
    }
    while (p) {
        --p;
        ppw.divide_by_two();
        x_ppw.divide_by_two();
        if (x_ppw <= rem) {
            rem -= x_ppw;
            quot += ppw;
        }
    }
    if (quot_neg) 
        quot.change_sign();
    if (rem_neg) 
        rem.change_sign();
    swap(need_rem ? rem : quot);        
}

bool big_integer::negative() const {
    return neg;
}

int big_integer::digit(size_t index) const {
    return num[index];
}

int& big_integer::digit(size_t index) {
    return num[index];
}    

size_t big_integer::length() const {
    return num.size();
}

std::string big_integer::to_string() const {
    std::string res;
    res.resize(length() + neg);
    if (neg)
        res[0] = '-';    
    for (size_t i = 0, j = length() - 1; i < length(); ++i, --j) {
        res[j + neg] = '0' + num[i];
    }
    return res;    
}

big_integer gcd(big_integer a, big_integer b) {
    return b ? gcd(b, a%b) : a;
}

big_integer big_integer::multiply_by_power_of_ten(size_t ppw) const {
    if (is_zero())
        return *this;
    big_integer res;
    res.neg = neg;
    res.num.resize(length() + ppw, 0);
    for (size_t i = 0; i < length(); ++i) {
        res.num[i + ppw] = num[i];
    }
    return res;
}

class rational {
private:
    big_integer num;
    big_integer den = big_integer(1);
    void format();
    void invert();
public:
    rational() : num(0_bi), den(1_bi) {}
    rational(const big_integer& num, const big_integer& den) : num(num), den(den) { format(); }
    rational(const big_integer& x) : num(x), den(1_bi) {}     
    rational(int x) : rational(big_integer(x)) {}
    rational(const rational& q) : num(q.num), den(q.den) {}
    ~rational() = default;
    rational& operator+=(const rational& q) &;
    rational& operator-=(const rational& q) &;
    rational& operator*=(const rational& q) &;
    rational& operator/=(const rational& q) &;
    rational& operator=(rational q) &;
    rational operator-() const;
    void change_sign();
    bool negative() const;
    std::string to_string() const;
    std::string as_decimal(size_t precision = 0) const;
    big_integer numerator() const;
    big_integer denominator() const;
    big_integer& numerator();
    big_integer& denominator();
    rational get_inverse() const;
    void swap(rational& q) &;
    size_t priority() const { return num.is_zero() ? 0 : num.length(); }
    explicit operator double() const {                                              
        double res = 0;
        bool was_point = 0;
        bool neg = 0;
        double ppw = 1;    
        std::string dec = as_decimal(100);
        for (char c : dec) {
            if (c == '-') {
                neg = 1;
                continue;        
            }    
            if (c == '.') {
                was_point = 1;
                continue;
            }
            double digit = static_cast<int>(c - '0');
            if (was_point) {
                ppw *= 0.1;
                res += ppw * digit;
            } else {
                res = res * 10 + digit;
            }
        }
        return res * (neg ? -1 : 1);
    }
    explicit operator bool() const { return !num.is_zero(); }
};

big_integer rational::numerator() const {
    return num;
}

big_integer rational::denominator() const {
    return den;
}

big_integer& rational::numerator() {
    return num;
}

big_integer& rational::denominator() {
    return den;
}
    
bool rational::negative() const {
    return num.negative();
}

void rational::change_sign() {
    num.change_sign();
}

rational rational::get_inverse() const {
    return rational(den, num);
}

void rational::swap(rational& q) & {
    std::swap(num, q.num);
    std::swap(den, q.den);
}

void rational::invert() {
    num.swap(den);
    if (den.negative()) {
        num.change_sign();
        den.change_sign();
    }
}

void rational::format() {
    big_integer g = gcd(abs(num), abs(den));
    if (g != 1_bi) {
        num /= g;
        den /= g;
    }
    if (den.negative()) {
        num.change_sign();
        den.change_sign();
    }
}

rational& rational::operator=(rational q) & {
    swap(q);
    return *this;
}

rational& rational::operator*=(const rational& q) & {
    num *= q.num;
    den *= q.den;
    format();
    return *this;    
}

rational& rational::operator/=(const rational& q) & {
    num *= q.den;
    den *= q.num;
    format();
    return *this;
}

rational& rational::operator+=(const rational& q) & {
    num *= q.den;
    num += den * q.num;
    den *= q.den;
    format();
    return *this;
}

rational& rational::operator-=(const rational& q) & {
    num *= q.den;
    num -= den * q.num;
    den *= q.den;
    format();
    return *this;    
}

rational rational::operator-() const {
    rational cp(*this);
    cp.change_sign();
    return cp;
}

rational operator+(rational lhs, const rational& rhs) {
    return lhs += rhs;
}

rational operator-(rational lhs, const rational& rhs) {
    return lhs -= rhs;     
} 

rational operator*(rational lhs, const rational& rhs) {
    return lhs *= rhs;
}

rational operator/(rational lhs, const rational& rhs) {
    return lhs /= rhs;    
}

bool operator<(const rational& lhs, const rational& rhs) {
    return lhs.numerator() * rhs.denominator() < lhs.denominator() * rhs.numerator();
}

bool operator>(const rational& lhs, const rational& rhs) {
    return rhs < lhs;
}

bool operator<=(const rational& lhs, const rational& rhs) {
    return !(rhs < lhs);
}

bool operator>=(const rational& lhs, const rational& rhs) {
    return !(lhs < rhs);
}

bool operator==(const rational& lhs, const rational& rhs) {
    return !(lhs < rhs || rhs < lhs);
}

bool operator!=(const rational& lhs, const rational& rhs) {
    return lhs < rhs || rhs < lhs;
}

std::string rational::to_string() const {
    return num.to_string() + (den.is_one() ? "" : "/" + den.to_string());
}

std::string rational::as_decimal(size_t precision) const {
    bool neg = negative();
    big_integer integer = abs(num);
    big_integer mantissa = (integer % den).multiply_by_power_of_ten(precision + 1) / den;
    integer /= den;
    if (mantissa.is_zero()) 
        return integer.to_string() + (precision ? "+" + std::string(precision, '0') : "");
    int need = 10 - mantissa.digit(0);
    if (need < 6) {
        mantissa += need;
        if (mantissa.length() > precision + 1)
            ++integer;    
    }
    std::string int_part = (neg ? "-" : "") + integer.to_string();
    if (!precision) 
        return int_part;
    std::string mant_part = mantissa.to_string();
    if (mant_part.size() <= precision)
        mant_part = std::string(precision + 1 - mant_part.size(), '0') + mant_part;
    return int_part + "." + mant_part.substr(0, precision);
}

rational abs(rational q) {
    if (q.negative())
        q.change_sign();
    return q;
}

//input big_integer and convert to rational
std::istream& operator>>(std::istream& in, rational& q) {
    return in >> q.numerator();
}

std::ostream& operator<<(std::ostream& out, const rational& q) {
    out << q.to_string();
    return out;
}

namespace details {

template <bool arg>
struct my_static_assert {
    static void check() {
        int* a = new int[arg ? 1 : -1];
        delete a;
    }
};

template <unsigned int number, unsigned int divisor>
struct factor_helper {
    static const unsigned int value = std::min(number % divisor == 0 ? divisor : number, 
                                               factor_helper<number, divisor - 1>::value);
};

template <unsigned int number>
struct factor_helper<number, 1> {
    static const unsigned int value = number;
};

template <unsigned int number>
struct square_root {
    static const unsigned int value = sqrt(number);        
};

//find smallest prime factor of number or 1 if there isn't any 
template <unsigned int number>
struct prime_factor {    
    static const unsigned int value = factor_helper<number, square_root<number>::value>::value;
};

template <unsigned int number, unsigned int p>
struct power_of_odd_prime_helper {
    static const bool value = (number % p == 0 && power_of_odd_prime_helper<number / p, p>::value);
};

template <unsigned int p> 
struct power_of_odd_prime_helper<1, p> {
    static const bool value = (p > 2);
};

template <unsigned int p> 
struct power_of_odd_prime_helper<0, p> {
    static const bool value = 0;
};

template <unsigned int number>
struct is_power_of_odd_prime {
    static const bool value = power_of_odd_prime_helper<number, prime_factor<number>::value>::value;    
};

template <unsigned int number>    
struct has_primitive_root {
    //p_power (checkk if number is power of negative prime * (1 or 2))
    static const bool p_power = is_power_of_odd_prime<(number % 2 == 0 ? number / 2 : number)>::value;
    static const bool value = (number == 2 || number == 4 || p_power);
};

unsigned int gcd(unsigned int a, unsigned int b) {
    return b ? gcd(b, a % b) : a;
}

unsigned int lcm(unsigned int a, unsigned int b) {
    return a / gcd(a, b) * b;
}

}; //namespace details

template <unsigned int number>
bool const has_primitive_root_v = details::has_primitive_root<number>::value;

template <unsigned int number>
bool const is_prime_v = (details::prime_factor<number>::value == number && number > 1);

template <unsigned int modulus>
class residue {
private:
    //factorization of Carmichael's function of modulus in (p, a) form
    //static std::vector<std::pair<unsigned int, unsigned int>> fac;
    //calculating Carmichael's function of modulus and factorizing it
    static unsigned int precalc(std::vector<std::pair<unsigned int, unsigned int>>& fac);
    unsigned int value = 0;
    static int const intmod = modulus;
public:    
    residue() : value(0) {}
    explicit residue(unsigned int x) : value(x % intmod) {}
    residue(int x) {
        x %= intmod;
        if (x < 0)
            x += intmod;
        value = x;
    }
    residue& operator+=(const residue& x) &;
    residue& operator-=(const residue& x) &;
    residue& operator*=(const residue& x) &;
    residue& operator/=(const residue& x) &;
    residue pow(long long x) const;
    //operator of exponentiation
    residue& operator^=(long long x) &;
    residue operator-() const;
    void change_sign() { *this = -(*this); } 
    bool operator==(const residue& x) const { return value == x.value; }
    bool operator!=(const residue& x) const { return value != x.value; }
    void swap(residue& x) &;    
    residue get_inverse() const;
    unsigned int order() const;
    static residue get_primitive_root();
    //priority for gaussian algorithm
    size_t priority() const { return value != 0; }
    explicit operator int() const { return value; }
    explicit operator bool() const { return value != 0; }
};     

template <unsigned int modulus>
void residue<modulus>::swap(residue<modulus>& x) & {
    std::swap(value, x.value);
}    

template <unsigned int modulus>
std::istream& operator>>(std::istream& in, residue<modulus>& x) { 
    int input;
    in >> input;
    x = residue<modulus>(input);
    return in; 
}

template <unsigned int modulus>
std::ostream& operator<<(std::ostream& out, const residue<modulus>& x) { 
    out << static_cast<int>(x);
    return out; 
}

template <unsigned int modulus> 
residue<modulus>& residue<modulus>::operator+=(const residue<modulus>& x) & {
    value += x.value;
    if (value >= modulus)
        value -= modulus;
    return *this;
}

template <unsigned int modulus>
residue<modulus> residue<modulus>::operator-() const {
    return residue(value ? modulus - value : 0);    
}

template <unsigned int modulus> 
residue<modulus>& residue<modulus>::operator-=(const residue<modulus>& x) & {
    return *this += -x;
}

template <unsigned int modulus> 
residue<modulus>& residue<modulus>::operator*=(const residue<modulus>& x) & {
    value = (value * 1ll * x.value) % modulus;
    return *this;                                 
}

template <unsigned int modulus>
residue<modulus> residue<modulus>::pow(long long x) const {
    residue res(1);
    residue a(*this);
    while (x) {
        if (x & 1)
            res *= a, x ^= 1;
        else
            a *= a, x >>= 1;            
    }
    return res;
}

template <unsigned int modulus> 
residue<modulus> residue<modulus>::get_inverse() const {
    details::my_static_assert<is_prime_v<modulus>>::check();
    return pow(modulus - 2);
}

template <unsigned int modulus> 
residue<modulus>& residue<modulus>::operator^=(long long x) & {
    *this = (x >= 0 ? pow(x) : get_inverse().pow(-x));
    return *this;
}

template <unsigned int modulus> 
residue<modulus>& residue<modulus>::operator/=(const residue<modulus>& x) & {
    return *this *= x.get_inverse();
}

template <unsigned int modulus>
residue<modulus> operator+(const residue<modulus>& a, const residue<modulus>& b) {
    residue<modulus> res(a);
    res += b;
    return res; 
}

template <unsigned int modulus>
residue<modulus> operator-(const residue<modulus>& a, const residue<modulus>& b) {
    residue<modulus> res(a);
    res -= b;
    return res; 
}

template <unsigned int modulus>
residue<modulus> operator/(const residue<modulus>& a, const residue<modulus>& b) {
    residue<modulus> res(a);
    res /= b;
    return res; 
}

template <unsigned int modulus>
residue<modulus> operator*(const residue<modulus>& a, const residue<modulus>& b) {
    residue<modulus> res(a);
    res *= b;
    return res; 
}

template <unsigned int modulus>
residue<modulus> operator^(const residue<modulus>& a, long long b) {
    residue<modulus> res(a);
    res ^= b;
    return res; 
}

//find Carmichael function of modulus using formula
//lambda = lcm({lambda(p1^a1), ...lambda(pk^ak))
//lambda(p^a) = phi(p^a) if (p is odd) or (a < 2) and phi(p^a) / 2 otherwise
template <unsigned int modulus>
unsigned int residue<modulus>::precalc(std::vector<std::pair<unsigned int, unsigned int>>& fac) {
    static std::vector<std::pair<unsigned int, unsigned int>> rfac;
    static unsigned int ans = 0;
    if (fac.size())
        return ans;
    if (ans) {
        fac = rfac;
        return ans;
    }
    ans = 1;
    unsigned int left = modulus;
    for (unsigned int p = 2, p_sq = 4; p_sq <= left; p_sq += p + p + 1, ++p) {
        if (left % p)
            continue;
        unsigned int p_a = 1, a = 0;
        while (left % p == 0) {
            p_a *= p;
            ++a;
            left /= p;
        }
        unsigned int phi = p_a - p_a / p;
        unsigned int lambda = (p > 2 || a < 3 ? phi : phi / 2);
        ans = details::lcm(ans, lambda);    
    }
    if (left > 1) 
        ans = details::lcm(ans, left - 1);
    left = ans;
    fac.reserve(30);    
    for (unsigned int p = 2, p_sq = 4; p_sq <= left; p_sq += p + p + 1, ++p) {
        unsigned int a = 0;
        while (left % p == 0) {
            ++a;
            left /= p;                    
        }
        if (a > 0)
            fac.push_back(std::make_pair(p, a));
    }
    if (left > 1) 
        fac.push_back(std::make_pair(left, 1));
    rfac = fac;    
    return ans;
}

//Iliminating prime factors of Carmichael function one by one
//Time complexity : O(sqrt(modulus) + Q * log^2(modulus)), where Q - number of calls of function
template <unsigned int modulus> 
unsigned int residue<modulus>::order() const {
    static residue const one(1);
    static std::vector<std::pair<unsigned int, unsigned int>> fac;
    if (details::gcd(value, modulus) != 1 || value == 0)
        return 0;
    if (value == 1)
        return 1;
    unsigned int ans = precalc(fac);
    for (auto [p, a] : fac) {
        for (unsigned int j = 1; j <= a; ++j) {
            unsigned int to_check = ans / p;
            if (pow(to_check) != one)
                break;
            ans = to_check;
        }
    }
    return ans;
}

//Brute force and check using factorization of Carmichael's function
//time complexity O(ans * log^2(modulus)) 
//~ O(log^6(modulus)) under the assumption of correctness of Riemann's hypothesis
template <unsigned int modulus> 
residue<modulus> residue<modulus>::get_primitive_root() {
    details::my_static_assert<has_primitive_root_v<modulus>>::check();
    static bool was = 0;
    static residue ans(1);
    static residue const one(1);
    static std::vector<std::pair<unsigned int, unsigned int>> fac;
    static unsigned int lambda = precalc(fac);
    if (modulus == 2)
        return one;
    if (was)
        return ans;    
    for (bool good = 0; !good; ){
        ans += one;
        if (details::gcd(ans.value, modulus) != 1)
            continue;
        good = 1;
        for (auto [p, a] : fac) {
            if (ans.pow(lambda / p) == one) {
                good = 0;
                break;
            }
        }
    }
    was = 1;
    return ans;
};

template <typename field>
using vv = std::vector<std::vector<field>>;

namespace details{

//divide(merge) matrix a into(from) four equal-sized blocks
//a - square matrix with side equal to some power of two     
template <typename field> 
void blocks(vv<field>& a, vv<field>& a_11, vv<field>& a_12, vv<field>& a_21, vv<field>& a_22, bool divide = 1) {
    size_t n = (divide ? a.size() / 2 : a_11.size());
    a_11.resize(n);
    a_12.resize(n);
    a_21.resize(n);
    a_22.resize(n);    
    a.resize(2 * n);
    for (int brow = 0; brow < 2; ++brow) {
        for (size_t i = 0; i < n; ++i) {
            a_11[i].resize(n);
            a_12[i].resize(n);
            a_21[i].resize(n);
            a_22[i].resize(n);    
            a[i + n * brow].resize(2 * n);
            auto a_it = a[i + n * brow].begin();
            auto it = (brow ? a_21[i].begin() : a_11[i].begin());
            auto end = (brow ? a_21[i].end() : a_11[i].end());
            while (it != end) {
                if (divide)    
                    *it = *a_it;
                else
                    *a_it = *it;
                ++a_it;
                ++it;
            }
            it = (brow ? a_22[i].begin() : a_12[i].begin());
            end = (brow ? a_22[i].end() : a_12[i].end());
            while (it != end) {
                if (divide)    
                    *it = *a_it;
                else
                    *a_it = *it;
                ++a_it;
                ++it;
            }
        }
    }                    
}

//caluclate += or -= for two matrices
template <typename field> 
void summation(vv<field>& a, const vv<field>& b, bool sub = 0) {
    auto begin_it_a = a.begin();
    auto begin_it_b = b.begin();
    while (begin_it_a != a.end()) {
        auto it_a = begin_it_a->begin();
        auto it_b = begin_it_b->begin();
        auto a_end = begin_it_a->end();
        while (it_a != a_end) {
            if (sub)
                *it_a -= *it_b;
            else
                *it_a += *it_b;
            ++it_a;
            ++it_b;
        }
        ++begin_it_a;
        ++begin_it_b;
    }
}


//calculate + or - of two matrices
template <typename field> 
vv<field> summ(const vv<field>& a, const vv<field>& b, bool sub = 0) {
    vv<field> res(a);
    summation(res, b, sub);
    return res;
}

template <typename field> 
vv<field> operator+(const vv<field>& a, const vv<field>& b) {
    return summ(a, b);
}

template <typename field> 
vv<field> operator-(const vv<field>& a, const vv<field>& b) {
    return summ(a, b, 1);
}

template <typename field> 
vv<field>& operator+=(vv<field>& a, const vv<field>& b) {
    summation(a, b);
    return a;
}

template <typename field> 
vv<field>& operator-=(vv<field>& a, const vv<field>& b) {
    summation(a, b, 1);
    return a;
}

template <typename field> 
vv<field> multiplication(vv<field>& a, vv<field>& b);

template <typename field> 
vv<field> operator*(vv<field>& a, vv<field>& b) {
    return multiplication(a, b);
}

//calculate * of two matrices using Vinograd-Strassen algorithm

//references are not const but a, b are not changeble
//that's for not_using const_cast in bad way or copypasting in merge/divide
template <typename field> 
vv<field> multiplication(vv<field>& a, vv<field>& b) {     
    if (a.size() == 1) 
        return vv<field>(1, std::vector<field>(1, a[0][0] * b[0][0]));
    using matrix_ = vv<field>;
    if (a.size() <= 64) {
        size_t n = a.size();
        vv<field> res(n, std::vector<field>(n));
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < n; ++j) {
                for (size_t k = 0; k < n; ++k) {
                    res[i][j] += a[i][k] * b[k][j];
                }
            }
        }
        return res;
    }         
    matrix_ a_11, a_12, a_21, a_22;
    matrix_ b_11, b_12, b_21, b_22;
    blocks<field>(a, a_11, a_12, a_21, a_22);
    blocks<field>(b, b_11, b_12, b_21, b_22);

    matrix_ s_1 = a_21 + a_22;
    matrix_ s_2 = s_1 - a_11;
    matrix_ s_3 = a_11 - a_21;
    matrix_ s_4 = a_12 - s_2;
    matrix_ s_5 = b_12 - b_11;
    matrix_ s_6 = b_22 - s_5;
    matrix_ s_7 = b_22 - b_12;
    matrix_ s_8 = s_6 - b_21;

    matrix_ p_1 = s_2 * s_6;
    matrix_ p_2 = a_11 * b_11;
    matrix_ p_3 = a_12 * b_21;
    matrix_ p_4 = s_3 * s_7;
    matrix_ p_5 = s_1 * s_5;
    matrix_ p_6 = s_4 * b_22;
    matrix_ p_7 = a_22 * s_8;

    matrix_ t_1 = p_1 + p_2;
    matrix_ t_2 = t_1 + p_4;

    matrix_ c_11 = p_2 + p_3;
    matrix_ c_12 = t_1 + p_5 + p_6;
    matrix_ c_21 = t_2 - p_7;
    matrix_ c_22 = t_2 + p_5;

    matrix_ c;
    blocks<field>(c, c_11, c_12, c_21, c_22, 0);
    return c;
}

template <typename field>
void concatenate(vv<field>& lhs, const vv<field>& rhs) {
    assert(lhs.size() == rhs.size());
    if (!lhs.size())
        return;
    size_t pref_col = lhs[0].size();
    size_t new_col = pref_col + rhs[0].size();
    for (size_t i = 0; i < lhs.size(); ++i) {
        lhs[i].resize(new_col);
        std::copy(rhs[i].begin(), rhs[i].end(), lhs[i].begin() + pref_col);
    }        
}

//make submatrix with rows [ln, rn) and columns [lm, rm)
template <typename field>
vv<field> submatrix(const vv<field>& a, size_t ln, size_t rn, size_t lm, size_t rm) {
    assert(ln <= rn && lm <= rm && rn <= a.size());
    vv<field> res(rn - ln);
    for (size_t i = ln; i < rn; ++i) {
        assert(rm <= a[i].size());
        res[i - ln] = std::vector<field>(a[i].begin() + lm, a[i].begin() + rm);    
    }
    return res;            
}


//there is two usage-types
// rank_type - calculate rank of (n x m) matrix
// !rank_type - take matix with size (n x n) of (n x 2n) 
// and calculate detrminant of left (n x n)
// by modifying it to unit form
// returns pair but only one value is guaranteed to be correct
template <typename field> 
std::pair<field, size_t> gauss(vv<field>& a, bool rank_type = 0) {                                     
    static const field zero(0);
    if (!a.size() || !a[0].size()) 
        return std::make_pair(zero, 0);
    size_t m = (rank_type ? a[0].size() : a.size());
    size_t n = a.size();
    field det(1);
    size_t rank = 0;
    for (size_t col = 0; col < m; ++col) {
        if (rank == n)
            return std::make_pair(det, rank);
        size_t pivot = rank;
        for (size_t row = rank + 1; row < n; ++ row) {
            if (a[row][col].priority() > a[pivot][col].priority()) 
                pivot = row;        
        }
        if (a[pivot][col].priority() == 0) {
            if (rank_type)
                continue;
            return std::make_pair(zero, 0);
        }
        if (pivot != rank) {
            std::swap(a[rank], a[pivot]);
            det.change_sign();
        }
        if (!rank_type)
            det *= a[rank][col];
        field inv = a[rank][col].get_inverse();
        for (auto it = a[rank].begin() + col; it != a[rank].end(); ++it) {
            if (*it)
                *it *= inv;    
        }
        for (size_t row = (rank_type ? rank + 1 : 0); row < n; ++row) {
            if (row == rank || !a[row][col])
                continue;
            const field& mlt = a[row][col];
            auto pivot_it = a[rank].begin() + col + 1;
                   for (auto it = a[row].begin() + col + 1; it != a[row].end(); ++it, ++pivot_it) {
                       if (*pivot_it) 
                           *it -= (*pivot_it) * mlt;    
            }
            a[row][col] = zero;    
        }
        ++rank;
    }
    return std::make_pair(det, rank);
}

}; // namespace details


template <size_t n, size_t m, typename field = rational> 
class matrix {
private:
    vv<field> a = vv<field>(n, std::vector<field>(m));
public:
    matrix() = default;
    explicit matrix(const field& x) : a(vv<field>(n, std::vector<field>(m, x))) {} 
    explicit matrix(const vv<field>& a) : a(a) {
        bool good_shape = (n == a.size());
        for (size_t i = 0; i < n && good_shape; ++i) {
            good_shape &= (a[i].size() == m);
        }
        assert(good_shape);
    }     
    matrix(std::initializer_list<std::vector<field>> a) : a(a) {} 
    matrix& operator+=(const matrix& q) &;
    matrix& operator-=(const matrix& q) &;
    matrix& operator*= (const field& q) &;
    matrix operator-() const;
    std::vector<field> get_row(size_t i) const { return a[i]; };
    std::vector<field> get_column(size_t j) const;
    const std::vector<field>& operator[](size_t i) const { return a[i]; }
    std::vector<field>& operator[](size_t i) { return a[i]; }
    matrix<m, n, field> transposed() const;    
    void swap(matrix& q) & { a.swap(q.a); }
    const vv<field>& data() const { return a; }
    size_t rank() const;
    //returns ellements if matrix in sz * sz format
    vv<field> squared_data(size_t sz) const;
    matrix& operator*=(const matrix<m, m, field>& q) &; 
    ///only for square matrices
    void transpose();
    field trace() const;         
    field det() const;      
    //return submatrix with rows [ln, ln + n1) and columns [lm, lm + m1)
    //operator of exponentiation
    matrix& operator^=(long long q) &;  
    void invert();                  
    matrix inverted() const;        
};
        
template <size_t n, typename field = rational>
using square_matrix = matrix<n, n, field>;

//unit matrix
template <size_t n, typename field = rational>
square_matrix<n, field> E() {
    square_matrix<n, field> res;
    for (size_t i = 0; i < n; ++i) {
        res[i][i] = static_cast<field>(1);
    }
    return res;
}

template <size_t n, size_t m, typename field> 
std::vector<field> matrix<n, m, field>::get_column(size_t j) const {
    std::vector<field> res(n);
    for (size_t i = 0; i < n; ++i) {
        res[i] = a[i][j];
    }
    return res;
}

//input in simple format , just spaces and '\n'
template <size_t n, size_t m, typename field> 
std::istream& operator>>(std::istream& in, matrix<n, m, field>& a) {
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < m; ++j) {
            in >> a[i][j];
        }
    }
    return in;
}

//output in (, , ) format
template <size_t n, size_t m, typename field> 
std::ostream& operator<<(std::ostream& out, matrix<n, m, field>& a) {
    for (size_t i = 0; i < n; ++i) {
        out << "(";
        for (size_t j = 0; j < m; ++j) {
            out << a[i][j] << (j + 1 == m ? ")\n" : ", ");
        }
    }
    return out;
}

template <size_t n, size_t m, typename field> 
bool operator==(const matrix<n, m, field>& a, const matrix<n, m, field>& b) {    
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < m; ++j) {
            if (a[i][j] != b[i][j])
                return 0;
        }
    }    
    return 1;
}

template <size_t n, size_t m, typename field> 
bool operator!=(const matrix<n, m, field>& a, const matrix<n, m, field>& b) {
    return !(a == b);
}

template <size_t n, size_t m, typename field> 
field matrix<n, m, field>::trace() const {
    details::my_static_assert<n == m>::check();
    field res(a[0][0]);
    for (size_t i = 1; i < n; ++i) {
        res += a[i][i];
    }
    return res;
}

template <size_t n, size_t m, typename field>
vv<field> matrix<n, m, field>::squared_data(size_t sz) const {
    vv<field> res(sz, std::vector<field>(sz));
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < m; ++j) {
            res[i][j] = a[i][j];
        }                          
    }
    return res;    
}

template <size_t n, size_t m, typename field> 
matrix<n, m, field> operator^(matrix<n, m, field> a, long long q) {
    details::my_static_assert<n == m>::check();
    matrix<n, m, field> res = E<n, field>();
    while (q) {
        if (q & 1) 
            res *= a, q ^= 1;
        else
            a *= a, q >>= 1; 
    }
    return res;
}
    
template <size_t n, size_t m, typename field> 
matrix<n, m, field>& matrix<n, m, field>::operator^=(long long q) & {
    details::my_static_assert<n == m>::check();
    *this = *this ^ q;
    return *this;
} 

template <size_t n, size_t m, typename field> 
matrix<m, n, field> matrix<n, m, field>::transposed() const {
    matrix<m, n, field> res;
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < m; ++j) {
            res[j][i] = a[i][j];
        }
    }
    return res;
}

template <size_t n, size_t m, typename field> 
void matrix<n, m, field>::transpose() {
    details::my_static_assert<n == m>::check();
    *this = transposed();
}

template <size_t n, size_t m, typename field> 
matrix<n, m, field>& matrix<n, m, field>::operator+=(const matrix<n, m, field>& q) & {
    details::operator+=(a, q.a);
    return *this;
}

template <size_t n, size_t m, typename field> 
matrix<n, m, field>& matrix<n, m, field>::operator-=(const matrix<n, m, field>& q) & {
    details::operator-=(a, q.a);
    return *this;
}

template <size_t n, size_t m, typename field> 
matrix<n, m, field> matrix<n, m, field>::operator-() const {
    matrix res;
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < m; ++j) {
            res[i][j] = -res[i][j];
        }
    }
    return res;
}

template <size_t n, size_t m, typename field> 
matrix<n, m, field>& matrix<n, m, field>::operator*=(const field& q) & {
    for (auto beg_it = a.begin(); beg_it != a.end(); ++beg_it) {
        for (auto it = beg_it->begin(); it != beg_it->end(); ++it) {
            *it *= q;
        }
    } 
    return *this;
}

template <size_t n, size_t m, typename field> 
matrix<n, m, field> operator*(const matrix<n, m, field>& a, const field& q) {
    matrix<n, m, field> res(a);
    res *= q;
    return res;
}

template <size_t n, size_t m, typename field> 
matrix<n, m, field> operator*(const field& q, const matrix<n, m, field>& a) {
    matrix<n, m, field> res(a);
    res *= q;
    return res;
}

template <size_t n, size_t m, typename field> 
matrix<n, m, field> operator+(const matrix<n, m, field>& a, const matrix<n, m, field>& b) {
    matrix<n, m, field> res(a);
    res += b;
    return res;    
} 

template <size_t n, size_t m, typename field> 
matrix<n, m, field> operator-(const matrix<n, m, field>& a, const matrix<n, m, field>& b) {
    matrix<n, m, field> res(a);
    res -= b;
    return res;    
} 

template <size_t n, size_t m1, size_t m2, typename field>
matrix<n, m1 + m2, field> operator|(const matrix<n, m1, field>& lhs, const matrix<n, m2, field>& rhs) {
    vv<field> res(lhs.data());
    concatenate(res, rhs.data());
    return matrix<n, m1 + m2, field>(res);
} 

template <size_t n, size_t m, typename field>                      
size_t matrix<n, m, field>::rank() const {
    vv<field> cp(a);                                
    return details::gauss(cp, 1).second;
}

template <size_t n, size_t m, typename field>                      
field matrix<n, m, field>::det() const {
    details::my_static_assert<n == m>::check();
    vv<field> cp(a);
    return details::gauss(cp).first;
}

template <size_t n, size_t m, typename field> 
matrix<n, m, field> matrix<n, m, field>::inverted() const {                  
    details::my_static_assert<n == m>::check();
    vv<field> matr = data();
    details::concatenate(matr, E<n, field>().data());
    assert(details::gauss(matr).first);
    matrix<n, m, field> res(details::submatrix(matr, 0, n, n, n + n));
    return res;
}

template <size_t n, size_t m, typename field> 
void matrix<n, m, field>::invert() {  
    *this = inverted();
}                

template <size_t n, size_t m, size_t k, typename field>
matrix<n, k, field> operator*(const matrix<n, m, field>& a, const matrix<m, k, field>& b) {
    size_t sz = 1;
    size_t mx = std::max(n, std::max(m, k));
    while (sz < mx) {
        sz <<= 1;
    }
    vv<field> ma = a.squared_data(sz);
    vv<field> mb = b.squared_data(sz);
    vv<field> res = details::operator*(ma, mb);
    res.resize(n);
    for (size_t i = 0; i < n; ++i) {
        res[i].resize(k);
    }
    return matrix<n, k, field>(res);                                                                         
}  

template <size_t n, size_t m, typename field> 
matrix<n, m, field>& matrix<n, m, field>::operator*=(const matrix<m, m, field>& a) & {
    *this = *this * a;
    return *this;                                                                  
}

#undef to_string 
#undef big_integer
#undef rational
#undef as_decimal
#undef residue
#undef get_inverse
#undef get_primitive_root
#undef matrix
#undef get_row
#undef get_column
#undef square_matrix