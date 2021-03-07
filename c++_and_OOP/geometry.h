#include <iostream>
#include <math.h>
#include <algorithm>
#include <iomanip>
#include <vector>
#include <cassert>


/*-------------------------------DISCLAIMER---------------------------//
In this code snake_case is used instead of camelCase, so defines

"class std::vector" and "class vector" are significantly DIFFERENT 

Warning: Some bitwise compression is used in "class polygon"
(don't know if it's necessary, but why not)

//--------------------------------------------------------------------*/ 
#define point Point
#define line Line
#define shape Shape
#define polygon Polygon
#define ellipse Ellipse
#define circle Circle
#define rectangle Rectangle
#define square Square
#define triangle Triangle
#define get_vertices getVertices
#define is_convex isConvex
#define circumscribed_circle circumscribedCircle
#define inscribed_circle inscribedCircle
#define Euler_line EulerLine
#define nine_points_circle ninePointsCircle
#define is_congruent_to isCongruentTo
#define is_similar_to isSimilarTo
#define contains_point containsPoint

double const eps = 1e-6;
double const PI = acos(-1);

double to_radians(double alpha){
    return alpha / 180 * PI;
}

bool equal(double a, double b){
    return std::abs(a - b) < eps;
}

double less(double a, double b) {
    return a < b && !equal(a, b);
}

bool less_equal(double a, double b) {
    return less(a, b) || equal(a, b);    
}

struct vector;
class line;

struct point{
    double x = 0;
    double y = 0;
    point() = default;
    point(double x, double y) : x(x), y(y) {}
    point& operator+=(const vector& v) &;
    point& operator-=(const vector& v) &;
    void reflex(const line& axis);
    void reflex(const point& center);
    void rotate(const point& center, double angle);
    void scale(const point& center, double coefficient);
    static double distance(const point& a, const point& b);
    static point middle(const point& a, const point& b);
    ~point() = default;
};

std::istream& operator>>(std::istream& in, point& p) {
    in >> p.x >> p.y;
    return in;
}

std::ostream& operator<<(std::ostream& out, const point& p) {
    out << p.x << ' ' << p.y;
    return out;
}

bool operator==(const point& lhs, const point& rhs) {
    return equal(lhs.x, rhs.x) && equal(lhs.y, rhs.y);
}

bool operator!=(const point& lhs, const point& rhs) {
    return !(lhs == rhs);
}

struct vector{
    double x = 0;                                            
    double y = 0;
    vector() = default;
    vector(double x, double y) : x(x), y(y) {}
    explicit vector(const point& p) : x(p.x), y(p.y) {}
    vector(const point& a, const point& b) : x(b.x - a.x), y(b.y - a.y) {}
    //normal for vector
    vector normal() const;
    double length() const;
    //normalizeing vector for given length
    void normalize(double len = 1);
    vector normalized(double len = 1) const;
    vector& operator*=(double k) &;
    vector& operator/=(double k) &;
    vector& operator+=(const vector& v) &;
    vector& operator-=(const vector& v) &;
    vector& operator&=(double angle) &;  //counterclockwise rotation by alpha (in radians)
    vector operator-() const;
    void reflex(const vector& v);         
    explicit operator point() const { return point(x, y); }
    void swap(point& p) & { std::swap(x, p.x); std::swap(y, p.y); }
    ~vector() = default;
};

point& point::operator+=(const vector& v) & {
    x += v.x;
    y += v.y;
    return *this;
}

point& point::operator-=(const vector& v) & {
    return *this += -v;
} 

void vector::normalize(double len) {
    if (equal(length(), 0)) {
        x = len;
        y = 0;
        return;
    }
    double k = len / length();
    *this *= k;
    return;
}

vector vector::normalized(double len) const {
    vector v(*this);
    v.normalize(len);
    return v;    
}

std::istream& operator>>(std::istream& in, vector& p) {
    in >> p.x >> p.y;
    return in;
}

std::ostream& operator<<(std::ostream& out, const vector& p) {
    out << p.x << ' ' << p.y;
    return out;
}

vector vector::normal() const{
    return vector(y, -x);
}

double vector::length() const {    
    return sqrt(x * x + y * y);
}

double point::distance(const point& a, const point& b) { 
    return vector(a, b).length(); 
}

//scalar multiplication of vectors       
double operator*(const vector& a, const vector& b) { 
    return a.x * b.x + a.y * b.y;
}

//pseudo-vector multiplicationof vectors
double operator^(const vector& a, const vector& b) { 
    return a.x * b.y - a.y * b.x;
}
            
vector& vector::operator*=(double k) & {
    x *= k;
    y *= k;
    return *this;
}

vector& vector::operator/=(double k) & {
    x /= k;
    y /= k;
    return *this;
}

vector operator*(double k, vector v) {
    return v *= k;
}

vector operator*(vector v, double k) {
    return v *= k;
}

vector operator/(vector v, double k) {
    return v /= k;
}

vector& vector::operator+=(const vector& v) & {
    x += v.x;
    y += v.y;
    return *this;
}

vector& vector::operator-=(const vector& v) & {
    x -= v.x;
    y -= v.y;
    return *this;
}

vector operator+(vector a, const vector& b) {
    return a += b;
}

vector operator-(vector a, const vector& b) {
    return a -= b;
}

vector& vector::operator&=(double angle) & {
    double nx = x * cos(angle) - y * sin(angle);
    double ny = x * sin(angle) + y * cos(angle);
    x = nx;
    y = ny;
    return *this;
}

vector operator&(vector v, double angle) {
    return v &= angle;
}

bool operator==(const vector& lhs, const vector& rhs) {
    return equal(lhs.x, rhs.x) && equal(lhs.y, rhs.y);
}

bool operator!=(const vector& lhs, const vector& rhs) {
    return !(lhs == rhs);
}

//check if two vector are collinear
bool operator||(const vector& lhs, const vector& rhs) {     
    return equal(0, lhs ^ rhs);
}

bool operator|(const vector& lhs, const vector& rhs) {     
    return (lhs || rhs) && !less(lhs * rhs, 0);        
}

vector vector::operator-() const {
    return vector(-x, -y);
}

vector bisector(vector a, const vector& b) {
    a.normalize(b.length());
    return (a + b);    
}

vector proection(const vector& v, const vector& to) {   
    return to * (v * to) / to.length() / to.length();
}

void point::reflex(const point& center) {
    *this -= 2 * vector(center, *this);
}

void point::scale(const point& center, double coefficient) {
    vector v(center, *this);
    v *= coefficient;
    *this = center;
    *this += v;
}

void point::rotate(const point& center, double angle) {
    vector v(center, *this);
    v &= angle;
    *this = center;
    *this += v;
}

double angle_between(const vector& a, const vector& b) {
    if (equal(std::min(a.length(), b.length()), 0))
        return 0;
    double cos_ = (a * b) / a.length() / b.length();
    cos_ = std::min(cos_, (double)1.0);
    cos_ = std::max(cos_, (double)-1.0);
    return acos(cos_);
}

point point::middle(const point& a, const point& b) {
    return static_cast<point>((static_cast<vector>(a) + static_cast<vector>(b)) / 2);
}

//angle of counter-clockwise rotation that moves a to b
double ccw_angle_between(const vector& a, const vector& b) {         
    double ang = angle_between(a, b);
    return (a & ang) == b ? ang : 2 * PI - ang;
}

double cw_angle_between(const vector& a, const vector& b) {
    double ang = angle_between(a, b);
    return (b & ang) == a ? ang : 2 * PI - ang;
}

bool clockwise(const vector& a, const vector& b) {
    return less_equal(a ^ b, 0);
}

class line{
private:
    point p;
    vector dir;
    void calc_canonic(double& a, double& b, double& c) const;
public:
    line() = default;
    line(const point& p, const vector& dir) : p(p), dir(dir.normalized()) {}
    line(const point& a, const point& b) : p(a), dir(vector(a, b)) {}
    line(const point& p, double k) : p(p), dir(vector(1, k)) {}
    line(double k, double b) : line(point(0, b), k) {}
    bool contains_point(const point& p) const;
    point random_point() const;
    vector direction() const;
    void rotate(const point& center, double angle);
    vector normal() const { return dir.normal(); } 
    void reflex(const line& axis);
    void reflex(const point& center);
    void scale(const point& center, double angle);
    bool intersect(const line& l, point& res) const;
    double distance(const point& from) const;
    ~line() = default;
};

void point::reflex(const line& axis) {
    *this -= 2 * proection(static_cast<vector>(*this) - static_cast<vector>(axis.random_point()), axis.normal());    
}

double angle_between(const line& a, const line& b) {       
    double ang = angle_between(a.direction(), b.direction());
    return 2 * ang > PI ? PI - ang : ang;
}

bool line::contains_point(const point& q) const {
    return vector(p, q) || dir;
}

point line::random_point() const {
    return p;
}
vector line::direction() const {
    return dir;
}

void line::calc_canonic(double& a, double& b, double& c) const {
    vector n = dir.normal();
    a = n.x;
    b = n.y;
    c = -(p.x * a + p.y * b);
}

void line::rotate(const point& center, double alpha) {
    p.rotate(center, alpha);
    dir &= alpha;
}

bool operator==(const line& lhs, const line& rhs) {
    return lhs.contains_point(rhs.random_point()) && (lhs.direction() || rhs.direction());    
}

//check if two lines are parallel
bool operator||(const line& lhs, const line& rhs) {  
    return lhs.direction() || rhs.direction();
}

bool line::intersect(const line& l, point& res) const {
    if (*this || l)
        return 0;
    double a, b, c;
    double la, lb, lc;
    calc_canonic(a, b, c);
    l.calc_canonic(la, lb, lc);
    double zn = normal() ^ l.normal();
    res.x = -(vector(c, b) ^ vector(lc, lb)) / zn;
    res.y = -(vector(a, c) ^ vector(la, lc)) / zn;
    return 1;
}

void line::reflex(const line& axis) {
    if (*this == axis) 
        return;
    if (*this || axis) {
        p.reflex(axis);
        return;             
    }
    point pint;
    intersect(axis, pint);
    if (p == pint) {
        p += dir;
    }
    p.reflex(axis);
    dir = vector(p, pint);    
}

void line::scale(const point& center, double coefficient) {
    point pp = p;
    pp += dir;
    p.scale(center, coefficient);
    pp.scale(center, coefficient);
    dir = vector(p, pp);
}

void line::reflex(const point& center) {
    point pp = p;
    pp += dir;
    vector v(p, center);
    p += 2 * v;
    pp += 2 * v;
    dir = vector(p, pp);
}

double line::distance(const point& from) const {     
    return proection(vector(p, from), normal()).length();
}

void vector::reflex(const vector& v) {
    point a;
    line l(a, v);
    point b = a;
    b += *this;
    b.reflex(l);
    *this = vector(a, b);
}
                                                                
line middle_perpendicular(const point& a, const point& b) {
    return line(point::middle(a, b), vector(a, b).normal());
}

class shape{
public:
    virtual double perimeter() const = 0;
    virtual long double area() const = 0;
    bool is_congruent_to(const shape& another) const;
    bool is_similar_to(const shape& another) const;   
    virtual bool contains_point(const point& p) const = 0;
    virtual void reflex(const point& center) = 0;     
    virtual void reflex(const line& axis) = 0;
    virtual void rotate(const point& center, double ang) = 0;
    virtual void scale(const point& p, double coefficient) = 0;
    virtual ~shape() = default;            
};

class ellipse: public shape{
protected:
    //focuses
    point f1;     
    point f2;
    //ellipse  
    double a_ = 0;  
    double b_ = 0;
public:
    ellipse() = default;
    ellipse(const point& f1, const point& f2, double sum_len) : f1(f1), f2(f2), a_(sum_len * 0.5) {
        double c_ = focus_distance();
        b_ = sqrt(a_ * a_ - c_ * c_);
    }
    ellipse(const ellipse& e) : f1(e.f1), f2(e.f2), a_(e.a_), b_(e.b_) {}    
    double perimeter() const final override;
    long double area() const final override;
    double a() const { return a_; }
    double b() const { return b_; }
    bool contains_point(const point& p) const final override;
    void reflex(const point& center) final override;
    void reflex(const line& axis) final override;
    void rotate(const point& center, double ang) final override;
    void scale(const point& p, double coefficient) final override;
    double focus_distance() const { return point::distance(f1, f2) / 2; }
    std::pair<point, point> focuses() const { return std::make_pair(f1, f2); }
    std::pair<line, line> directrices() const;
    double eccentricity() const { return focus_distance() / a_; }
    point center() const { return point::middle(f1, f2); }    
    ~ellipse() = default;
};

double ellipse::perimeter() const {
    return PI * ((a_ + b_) * 3 - sqrt((a_ * 3 + b_) * (a_ + b_ * 3)));
}

long double ellipse::area() const {
    return PI * a_ * b_;
}

std::pair<line, line> ellipse::directrices() const {
    vector v(f2, f1);
    v.normalize(a_ / eccentricity());
    point c = center();
    c += v;
    line d1(c, v.normal());
    c -= v + v;
    line d2(c, v.normal());
    return std::make_pair(d1, d2);        
}    

void ellipse::rotate(const point& center, double angle) {
    angle = to_radians(angle);
    f1.rotate(center, angle);
    f2.rotate(center, angle);    
}
                                            
void ellipse::reflex(const point& center) {
    f1.reflex(center);
    f2.reflex(center);
}
    
void ellipse::reflex(const line& axis) {
    f1.reflex(axis);
    f2.reflex(axis);
}

void ellipse::scale(const point& center, double coefficient) {
    f1.scale(center, coefficient);
    f2.scale(center, coefficient);
    a_ *= std::abs(coefficient);
    b_ *= std::abs(coefficient);    
}

bool ellipse::contains_point(const point& p) const {
    return less_equal(point::distance(p, f1) + point::distance(p, f2), 2 * a_); 
}

class polygon: public shape{
protected:
    std::vector<point> points;
    ///mask of calculated quantities
    //0-th bit - is perimeter calculated
    //1-th bit - is area calculated
    //2-th bit - is convexity calculated
    //3-th bit - is polygon convex
    mutable char mask = 0;
    mutable double perimeter_ = 0;
    mutable long double area_ = 0;
public:
    polygon() = default;
    polygon(const std::vector<point>& points) : points(points) {}
    polygon(const std::initializer_list<point>& lst) : points(lst) {} 
    double perimeter() const final override;
    long double area() const final override;
    bool is_convex() const;    
    bool contains_point(const point& p) const final override;
    size_t vertices_count() const { return points.size(); }
    std::vector<point> get_vertices() const { return points; }
    void reflex(const point& center) final override;
    void reflex(const line& axis) final override;
    void rotate(const point& center, double angle) final override;
    void scale(const point& center, double coefficient) final override;
    point operator[](size_t index) const { return points[index]; }
    bool check_coefficient(const std::vector<point>& v, size_t pivot, double k) const;
    ~polygon() = default;
};

bool polygon::contains_point(const point& p) const {                  
    double delta = 0;
    for (auto pr = points.begin(), it =  pr + 1; pr != points.end(); ++it, ++pr) {
        if (it == points.end())
            it = points.begin();
        if (p == *it)
            return 1;
        vector vpr(p, *pr);
        vector vit(p, *it);
        //lies on side of polygon
        if ((vpr || vit) && !(vpr | vit))     
            return 1;
        double ang = angle_between(vpr, vit) * (clockwise(vpr, vit) ? 1 : -1);
        delta += ang;
    }
    return abs(delta) > PI;
}

void polygon::reflex(const point& center) {
    for (point& p: points) {
        p.reflex(center);    
    }
}

double polygon::perimeter() const {
    if (mask & 1)
        return perimeter_;
    mask |= 1;
    perimeter_ = 0;
    for (auto it = points.begin(), nxt = it + 1; it != points.end(); ++it, ++nxt) {
        if (nxt == points.end())
            nxt = points.begin();
        perimeter_ += point::distance(*it, *nxt);
    }
    return perimeter_;    
}

long double polygon::area() const {
    if (mask & 2)
        return area_;
    mask |= 2;
    area_ = 0;
    auto pr = points.end() - 1;
    for (auto it = points.begin(); it != points.end(); ++it, ++pr) {
        if (pr == points.end())
            pr = points.begin();
        long double dx = pr->x - it->x;
        long double dy = pr->y + it->y;
        area_ += dx * dy;
    }
    area_ = abs(area_) * 0.5;
    return area_;
}

void polygon::reflex(const line& axis) {
    for (point& p: points) {
        p.reflex(axis);
    }
}

void polygon::rotate(const point& center, double angle) {
    angle = to_radians(angle);
    for (point& p : points) {
        p.rotate(center, angle);
    }
}

bool polygon::is_convex() const {              
    //magic of bitwise comprsession
    //it will apear in some functions :)
    if (mask & 4)           
        return mask & 8;
    mask |= 4;
    bool was[2] = {};
    for (auto fir = points.begin(), sec = fir + 1, thi = sec + 1; fir != points.end(); ++fir, ++sec, ++thi) {
        if (sec == points.end())
            sec = points.begin();
        if (thi == points.end())
            thi = points.begin();
        was[clockwise(vector(*fir, *sec), vector(*fir, *thi))] = 1;    
    }
    bool ans = was[0] ^ was[1];
    if (ans) 
        mask |= 8;    
    return ans;
}

void polygon::scale(const point& center, double coefficient) {
    for (point& p : points) {
        p.scale(center, coefficient);
    }
    if (mask & 1)
        perimeter_ *= abs(coefficient);
    if (mask & 2)
        area_ *= coefficient * coefficient;
}

class circle: public ellipse{
public:    
    circle() = default;
    circle(const point& center, double radius) {
        f1 = f2 = center;
        a_ = b_ = radius;
    }
    double radius() const { return a_; }
    ~circle() = default;
};

class rectangle: public polygon{
public:
    rectangle() = default;
    rectangle(const point& pa, const point& pb, double k) {                      
        points.resize(4);
        points[0] = pa;
        points[2] = pb;
        vector c(pa, pb);
        if (k < 1)
            k = 1 / k;
        vector a((k * c.y + c.x) / (k * k + 1), (c.y - k * c.x) / (k * k + 1));
        if (!clockwise(a, c)) 
            a.reflex(c);
        points[1] = pa;
        points[0] += a;
        points[3] = pa;
        points[3] += c - a;            
    }
    point center() const { return point::middle(points[0], points[2]); }
    std::pair<line, line> diagonals() const; 
    ~rectangle() = default;
};

std::pair<line, line> rectangle::diagonals() const { 
    return std::make_pair(line(points[0], points[2]), line(points[1], points[3])); 
}

class square: public rectangle{
public:
    square() = default;
    square(const point& pa, const point& pb) : rectangle(pa, pb, 1) {}
    circle circumscribed_circle() const { return circle(center(), point::distance(points[0], points[2]) / 2); }
    circle inscribed_circle() const { return circle(center(), point::distance(points[0], points[1]) / 2); }
    ~square() = default;
};                  

class triangle: public polygon{
public:
    triangle() = default;
    triangle(const point& a, const point& b, const point& c) : polygon({a, b, c}) {}
    circle inscribed_circle() const;    
    circle circumscribed_circle() const;  
    circle nine_points_circle() const;
    point orthocenter() const;
    point centroid() const;
    line Euler_line() const { return line(centroid(), orthocenter()); } 
    ~triangle() = default;
};

point triangle::centroid() const {
    vector va(points[0]);
    vector vb(points[1]);
    vector vc(points[2]); 
    return static_cast<point>((va + vb + vc) / 3); 
}
    

point triangle::orthocenter() const {
    line ha(points[0], vector(points[1], points[2]).normal());
    line hb(points[1], vector(points[0], points[2]).normal());
    point res;
    ha.intersect(hb, res);
    return res;
}

circle triangle::nine_points_circle() const {
    point ma = point::middle(points[0], points[1]);
    point mb = point::middle(points[1], points[2]);
    point mc = point::middle(points[0], points[2]);
    return triangle(ma, mb, mc).circumscribed_circle();
} 

circle triangle::inscribed_circle() const {
    line bis_a(points[0], bisector(vector(points[0], points[1]), vector(points[0], points[2])));
    line bis_b(points[1], bisector(vector(points[1], points[0]), vector(points[1], points[2])));
    point cen;
    bis_a.intersect(bis_b, cen);
    return circle(cen, line(points[0], points[1]).distance(cen));
}

circle triangle::circumscribed_circle() const {
    line mid_per_a = middle_perpendicular(points[1], points[2]);
    line mid_per_b = middle_perpendicular(points[0], points[2]);
    point cen;
    mid_per_a.intersect(mid_per_b, cen);
    return circle(cen, point::distance(points[0], cen));
}

int are_both_ellipses(const shape& lhs, const shape& rhs) {
    int ellipses = 0;
    try {
        dynamic_cast<const polygon&>(lhs);        
    } catch (const std::bad_cast& er) {
        ++ellipses;
    }
    try {
        dynamic_cast<const polygon&>(rhs);         
    } catch (const std::bad_cast& er) {
        ++ellipses;
    }
    if (ellipses == 1)
        return -1;
    return ellipses ? 1 : 0;    
}

bool congruent_ellipses(const ellipse& lhs, const ellipse& rhs) {                            
    return equal(lhs.a(), rhs.a()) && equal(lhs.b(), rhs.b());    
}

bool polygon::check_coefficient(const std::vector<point>& v, size_t pivot, double k) const {                  
    size_t n = v.size();
    size_t lpr = n - 1;
    size_t lcur = 0;
    size_t lnxt;
    size_t rpr = (pivot ? pivot : n) - 1;
    size_t rcur = pivot;
    size_t rnxt;
    bool good = 1;
    for(; lcur < n && good; lpr = lcur, ++lcur, rpr = rcur, rcur = rnxt) {
        lnxt = (lcur + 1 == n ? 0 : lcur + 1);
        rnxt = (rcur + 1 == n ? 0 : rcur + 1);
        vector la(points[lcur], points[lpr]);
        vector lb(points[lcur], points[lnxt]);
        vector ra(v[rcur], v[rpr]);
        vector rb(v[rcur], v[rnxt]);
        good &= equal(la.length() * k, ra.length());
        good &= equal(lb.length() * k, rb.length());
        good &= equal(ccw_angle_between(la, lb), ccw_angle_between(ra, rb));
    }    
    return good;
}
    
bool congruent_polygons(const polygon& lhs, const polygon& rhs) {                                     
    size_t n = lhs.vertices_count();
    //different number of points
    if (n != rhs.vertices_count())
        return 0;
    std::vector<point> rhscp = rhs.get_vertices();
    for (int rev = 0; rev < 2; ++rev) {
        if (rev) 
            reverse(rhscp.begin(), rhscp.end());
        for (size_t pivot = 0; pivot < n; ++pivot) {
            if (lhs.check_coefficient(rhscp, pivot, 1))
                return 1;
        }
    }
    return 0;
}

bool similar_polygons(const polygon& lhs, const polygon& rhs) {                             
    size_t n = lhs.vertices_count();
    //different number of points
    if (n != rhs.vertices_count())
        return 0;
    std::vector<point> rhscp = rhs.get_vertices();
    double llen = point::distance(lhs[1], lhs[0]);
    for (int rev = 0; rev < 2; ++rev) {
        if (rev) 
            reverse(rhscp.begin(), rhscp.end());
        for (size_t pivot = 0; pivot < n; ++pivot) {
            double rlen = point::distance(rhscp[pivot], rhscp[(pivot + 1 == n ? 0 : pivot + 1)]);
            if (lhs.check_coefficient(rhscp, pivot, rlen / llen))
                return 1;                               
        }
    }        
    return 0;
}

bool similar_ellipses(const ellipse& lhs, const ellipse& rhs) {                                            
    return equal(lhs.a() / rhs.a(), lhs.b() / rhs.b());
}

bool shape::is_congruent_to(const shape& another) const {                                   
    int tps = are_both_ellipses(*this, another);
    if (tps == -1)
        return 0;
    if (tps == 0) 
        return congruent_polygons(dynamic_cast<const polygon&>(*this), dynamic_cast<const polygon&>(another));
    return congruent_ellipses(dynamic_cast<const ellipse&>(*this), dynamic_cast<const ellipse&>(another));
}

bool shape::is_similar_to(const shape& another) const {                           
    int tps = are_both_ellipses(*this, another);
    if (tps == -1)
        return 0;
    if (tps == 0) 
        return similar_polygons(dynamic_cast<const polygon&>(*this), dynamic_cast<const polygon&>(another));
    return similar_ellipses(dynamic_cast<const ellipse&>(*this), dynamic_cast<const ellipse&>(another));    
}

bool operator==(const polygon& lhs, const polygon& rhs) {
    size_t n = lhs.vertices_count();
    //different number of points
    if (n != rhs.vertices_count())
        return 0;
    size_t pivot = n;
    point check = lhs[0];
    for (size_t i = 0; i < n; ++i) {
        if (check == rhs[i]) {
            pivot = i;
            break; 
        }
    }          
    //no matching point 
    if (pivot == n) 
        return 0;
    //moving in one direction
    bool good = 1;
    for (size_t li = 0, ri = pivot; li < n && good; ++li, ri = (ri + 1 == n ? 0 : ri + 1)) {
        good &= (lhs[li] == rhs[ri]);
    }
    if (good)
        return 1;
    //oposite direction
    good = 1;
    for (size_t li = 0, ri = pivot; li < n && good; ++li, ri = (ri == 0 ? n - 1 : ri - 1)) {
        good &= (lhs[li] == rhs[ri]);
    }
    return good;
}

bool operator==(const ellipse& a, const ellipse& b) {
    if (!a.is_similar_to(b))
        return 0;
    std::pair<point, point> fa = a.focuses();
    std::pair<point, point> fb = b.focuses();
    if (fa.first == fb.first && fa.second == fb.second)
        return 1;
    if (fa.first == fb.second && fa.second == fb.first)
        return 1;
    return 0;
}

bool operator==(const shape& lhs, const shape& rhs) {                           
    int tps = are_both_ellipses(lhs, rhs);
    if (tps == -1)
        return 0;
    if (tps == 0)
        return dynamic_cast<const polygon&>(lhs) == dynamic_cast<const polygon&>(rhs);    
    return dynamic_cast<const ellipse&>(lhs) == dynamic_cast<const ellipse&>(rhs); 
}    

bool operator!=(const shape& lhs, const shape& rhs) {
    return !(lhs == rhs);
}

#undef point
#undef line
#undef shape
#undef polygon
#undef ellipse
#undef circle
#undef rectangle
#undef square
#undef triangle
#undef get_vertices
#undef is_convex
#undef circumscribed_circle
#undef inscribed_circle
#undef Euler_line
#undef nine_points_circle
#undef is_congruent_to
#undef is_similar_to
#undef contains_point