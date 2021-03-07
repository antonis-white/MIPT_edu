#include<iostream>
#include<vector>
#include<string>

//I prefer using snake_case 
#define to_string toString 
#define big_integer BigInteger
#define rational Rational 
#define as_decimal asDecimal 
 
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
    big_integer den;
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
    rational inverse() const;
    void swap(rational& q) &;
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

rational rational::inverse() const {
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

std::ostream& operator<<(std::ostream& out, const rational& q) {
    out << q.to_string();
    return out;
}

#undef to_string 
#undef big_integer
#undef rational
#undef as_decimal
