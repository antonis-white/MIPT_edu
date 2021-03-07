#include<iostream>
#include<cstring>

size_t max(size_t a, size_t b) {
    return a > b ? a : b;
}

size_t min(size_t a, size_t b) {
    return a < b ? a : b;
}

class String {
  private:
    size_t len = 0;
    size_t buf_sz = 0;
    char* str = nullptr;
    static const size_t MIN_BUF_SZ = 4;
    static const int EXP_COEF = 2;
    void check_buffer(size_t clen);      
  public:
    String() = default;
    String(const char& ch) { push_back(ch); }
    String(size_t len, char ch) { resize(len, ch); }
    String(const String& par) {
        resize(par.len);
        replace(par, 0, 0, par.len);            
    }
    String(const char* c_str) {
        for (size_t i = 0; c_str[i]; ++i) {
            push_back(c_str[i]);
        }    
    }       
    ~String() { delete[] str; }
    char operator[](size_t index) const;
    char& operator[](size_t index);
    String& operator=(const String& s);
    void resize(size_t nlen, char ch = '\0');        
    void assign(size_t start, size_t count, char ch); 
    void assign(char ch);                         
    void replace(const String& from, size_t to_start, size_t from_start, size_t count);
    void reverse(size_t from, size_t to);        
    void reverse();           
    size_t length() const;      
    char& front();          
    char& back();
    char front() const;
    char back() const;                
    void push_back(char ch);   
    void pop_back();               
    size_t find(const String& pat) const;
    size_t rfind(const String& pat) const;
    String substr(size_t start, size_t count) const;
    bool empty() const;            
    void clear();         
    void swap(String& s);
};

size_t max_Pi_fun(const String& s, size_t& where) {
    size_t* pi = new size_t[s.length()];
    pi[0] = 0;
    where = 0;
    size_t ans = 0;
    for (size_t i = 1; i < s.length(); ++i) {
        size_t cur_pi = pi[i - 1];
        while (cur_pi > 0 && s[i] != s[cur_pi]) {
            cur_pi = pi[cur_pi - 1];
        }
        if (s[i] == s[cur_pi])
            ++cur_pi;
        if (cur_pi > ans)
            ans = cur_pi, where = i;          
        pi[i] = cur_pi;
    }    
    delete[] pi;
    return ans;
}

char String::operator[](size_t index) const {
    return str[index];
}

char& String::operator[](size_t index) {
    return str[index];
}               

String& String::operator=(const String& s) {
    String cp(s);
    swap(cp);
    return *this;
}

std::istream& operator>>(std::istream &in, String& s) {         
    s.clear();
    char inp;
    while (!in.fail()) {
        inp = in.get();
        if (!isspace(inp)) {
            s.push_back(inp);
            break;
        }
    }
    while (!in.fail()) {
        inp = in.get();
        if (inp == EOF || isspace(inp))
            break;
        s.push_back(inp);
    }
    return in;
}

std::ostream& operator<<(std::ostream &out, const String& s) {
    for (size_t i = 0; i < s.length(); ++i) {
        out.put(s[i]);
    }
    return out;
}
    
bool operator<(const String& a, const String& b) {
    size_t a_it = 0;
    size_t b_it = 0;
    for (; a_it < a.length() && b_it < b.length(); ++a_it, ++b_it) {
        if (a[a_it] < b[b_it])
            return 1;
        if (a[a_it] > b[b_it])
            return 0;        
    }
    return a_it == a.length() && b_it < b.length();
}

bool operator>(const String& a, const String& b) {
    return b < a;
}

bool operator<=(const String& a, const String& b) {
    return !(a > b);
}

bool operator>=(const String& a, const String& b) {
    return !(a < b);
}

bool operator==(const String& a, const String& b) {
    return !(a < b) && !(b < a);
}

bool operator!=(const String a, const String b) {
    return !(a == b);
}
                                     
String& operator+=(String& a, const String& b) {
    size_t plen = a.length();
    a.resize(a.length() + b.length());
    a.replace(b, plen, 0, b.length());
    return a;    
}

String operator+(const String& a, const String& b) {
    String res(a);
    res += b;
    return res;    
}
    
bool String::empty() const {
    return len == 0;
}

size_t String::length() const {
    return len;
}

void String::clear() {
    if (!str)
        return;
    len = buf_sz = 0;
    delete[] str;
    str = nullptr;    
}

void String::assign(char ch) {
    assign(0, len, ch);
}

void String::assign(size_t start, size_t count, char ch) {
    memset(str + start, ch, count);
}

void String::replace(const String& from, size_t to_start, size_t from_start, size_t count) {
    if (count)
        memcpy(str + to_start, from.str + from_start, count);    
}

char& String::front() {
    return str[0];
}

char& String::back() {
    return str[len - 1];            
}

char String::front() const {
    return str[0];
}

char String::back() const {
    return str[len - 1];
}

void String::push_back(char ch) {
    check_buffer(len + 1);
    str[len++] = ch;
}

void String::pop_back() {
    check_buffer(--len);
}

void String::resize(size_t nlen, char ch) {
    if (nlen == len)
        return;
    if (nlen > len) {
        check_buffer(nlen);
        size_t plen = len;
        len = nlen;
        assign(plen, nlen - plen, ch);
        return;    
    }        
    len = nlen;
    check_buffer(len);            
}

void String::check_buffer(size_t clen) {
    size_t new_buf_sz = buf_sz;
    if (clen > buf_sz || clen * EXP_COEF * EXP_COEF <= buf_sz)
        new_buf_sz = max(MIN_BUF_SZ, clen * EXP_COEF);
    if (buf_sz == new_buf_sz)
        return;
    char* new_buf = new char[new_buf_sz];
    size_t to_copy = min(len, new_buf_sz);     
    if (to_copy)
        memcpy(new_buf, str, to_copy);
    delete[] str;
    str = new_buf;
    buf_sz = new_buf_sz;
}

String String::substr(size_t start, size_t count) const {
    String res;
    res.resize(count);
    res.replace(*this, 0, start, count);
    return res;    
}

void String::reverse(size_t from, size_t to) {
    for (size_t l_it = from, r_it = to; l_it + 1 <= r_it; ++l_it, --r_it) {
        std::swap(str[l_it], str[r_it]);
    }
}

void String::reverse() {
    if (len > 1)
        reverse(0, len-1);    
}

void String::swap(String& s) {
    std::swap(len, s.len);
    std::swap(buf_sz, s.buf_sz);
    std::swap(str, s.str);
}

size_t String::find(const String& pat) const {
    size_t ind = len;
    size_t max_pi = max_Pi_fun(pat + '\0' + *this, ind);
    if (max_pi != pat.length())
        return len;
    return ind - 2*pat.length();                
}

size_t String::rfind(const String& pat) const {
    String cp(*this);
    String rpat(pat);
    cp.reverse();
    rpat.reverse();
    size_t ind = cp.find(rpat);
    if (ind == len)
        return len;
    return len - ind - pat.length();
}