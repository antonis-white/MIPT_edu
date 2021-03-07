#include <iostream>
#include <vector>
#include <math.h>
#include <utility>

#define residue Residue
#define get_inverse getInverse
#define get_primitive_root getPrimitiveRoot

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
    explicit residue(int x) {
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
    bool operator==(const residue& x) const { return value == x.value; }
    bool operator!=(const residue& x) const { return value != x.value; }
    void swap(residue& x) &;    
    residue get_inverse() const;
    unsigned int order() const;
    static residue get_primitive_root();
    explicit operator int() const { return value; }
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

#undef residue
#undef get_inverse
#undef get_primitive_root