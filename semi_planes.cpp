#include <vector>
#include <iostream>
#include <deque>
#include <algorithm>
#include <iomanip>
#include <cmath>
#include <cstdio>
//#include "point.h"

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

using namespace std;

const long double eps = 1e-10;

typedef ptld pt;

int sgn(const long double &a) {
    if (a > eps)
        return 1;
    else if (a < -eps)
        return -1;
    else
        return 0;
}

pt rotate(const pt &p) {
    return pt(-p.y, p.x);
}

struct halfPlane {
    pt a, dir;
    halfPlane() {}
    halfPlane(const pt &_a, const pt &_dir) : a(_a), dir(_dir) {}

    int type() const {
        if (sgn(dir.y) != 0)
            return sgn(dir.y);
        return sgn(dir.x);
    }

    bool operator<(const halfPlane &rhs) const {
        if (type() != rhs.type())
            return type() < rhs.type();
        long double cr = cross(dir, rhs.dir);
        if (sgn(cr) != 0) {
            return sgn(cr) < 0;
        }
        return !include(rhs.a);
    }

    bool operator==(const halfPlane &rhs) const {
        if (type() != rhs.type())
            return false;
        long double cr = cross(dir, rhs.dir);
        return sgn(cr) == 0;
    }

    bool include(const pt &p) const {
        return sgn(cross(dir, p - a)) >= 0;
    }
};

pt getIntersection(const halfPlane &lhs, const halfPlane &rhs) {
    long double u = cross(lhs.a - rhs.a, rhs.dir);
    long double v = cross(lhs.dir, rhs.dir);
    long double t = -u / v;
    return lhs.a + lhs.dir * t;
}

const int N = 111111;

halfPlane hp[N];

int n;

int bot, top;

inline void intersect() {
    top = bot = 0;

    for (int i = 0; i < n; i++) {
        while (bot - top >= 2 && !hp[i].include(getIntersection(hp[bot - 1], hp[bot - 2])))
            bot--;
        while (bot - top >= 2 && !hp[i].include(getIntersection(hp[top], hp[top + 1])))
            top++;
        hp[bot++] = hp[i];
    }
    while (bot - top >= 2 && !hp[top].include(getIntersection(hp[bot - 1], hp[bot - 2])))
        bot--;
    while (bot - top >= 2 && !hp[bot - 1].include(getIntersection(hp[top], hp[top + 1])))
        top++;
}

int main() {
    while (scanf("%d", &n) == 1) {
        int sz = 0;
        for (int i = 0; i < n; ++i) {
            //double a = nxt(), b = nxt(), c = nxt(), d = nxt();
            double a, b, c, d;
            scanf("%lf%lf%lf%lf", &a, &b, &c, &d);
            pt A(a, b);
            pt dir(c - a, d - b);
            dir = dir * (1.0 / dir.len());
            hp[sz].a = A;
            hp[sz++].dir = dir;
        }
        n = sz;

        hp[n++] = halfPlane(pt(0, 0), pt(10000, 0));
        hp[n++] = halfPlane(pt(10000, 0), pt(0, 10000));
        hp[n++] = halfPlane(pt(10000, 10000), pt(-10000, 0));
        hp[n++] = halfPlane(pt(0, 10000), pt(0, -10000));

        sort(hp, hp + n);
        n = (int) (unique(hp, hp + n) - hp);

        intersect();

        if (bot - top < 2) {
            cout << "0.0" << endl;
            return 0;
        }

        long double sq = 0;

        pt pts[bot - top];
        sz = 0;

        for (int i = top; i < bot; ++i) {
            int k = i + 1;
            if (k == bot) {
                k = top;
            }
            pts[sz++] = getIntersection(hp[i], hp[k]);
        }

        for (int i = 0; i < sz; ++i) {
            sq += cross(pts[i], pts[(i + 1) % sz]);
        }
        sq = abs(sq);
        sq /= 2;
        cout << setprecision(1) << fixed;
        printf("%.1f\n", (double)sq);
    }
    return 0;
}