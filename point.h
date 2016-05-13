#ifndef TEAM_NOTEBOOK_INNOPOLIS_U_POINT_H
#define TEAM_NOTEBOOK_INNOPOLIS_U_POINT_H


template <class T>
struct point {
    T x, y;
    point() : x(), y() {}
    point(T x, T y) : x(x), y(y) {}
    point operator + (const point &r) const {
        return point(x + r.x, y + r.y);
    }
    point operator - (const point &r) const {
        return point(x - r.x, y - r.y);
    }
    point operator * (const T &r) const {
        return point(x * r, y * r);
    }
    point rot(T co, T si) const {
        return point(x * co - y * si, x * si + y * co);
    }
    point rot(T ang) const {
        return rot(cos(ang), sin(ang));
    }
    T sqlen() const {
        return abs(x * x + y * y);
    }
    long double len() const {
        return sqrtl(sqlen());
    }
};

template <class T>
T dot(const point<T> &l, const point<T> &r) {
    return l.x * r.x + l.y * r.y;
}

template <class T>
T cross(const point<T> &l, const point<T> &r) {
    return l.x * r.y - l.y * r.x;
}

typedef point<int> pti;
typedef point<long long> ptl;
typedef point<double> ptd;
typedef point<long double> ptld;


#endif //TEAM_NOTEBOOK_INNOPOLIS_U_POINT_H
