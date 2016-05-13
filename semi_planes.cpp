#include <vector>
#include <iostream>
#include <deque>
#include <algorithm>
#include <iomanip>
#include <cmath>
#include <cstdio>
#include "point.h"

const double eps = 1e-9;
using namespace std;

typedef ptd pt;

int sgn(const double &a) {
    if (a > eps)
        return 1;
    else if (a < -eps)
        return -1;
    else
        return 0;
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
        double cr = cross(dir, rhs.dir);
        if (sgn(cr) != 0) {
            return sgn(cr) < 0;
        }
        return !include(rhs.a);
    }

    bool operator==(const halfPlane &rhs) const {
        if (type() != rhs.type())
            return false;
        double cr = cross(dir, rhs.dir);
        return sgn(cr) == 0;
    }

    bool include(const pt &p) const {
        return sgn(cross(dir, p - a)) >= 0;
    }
};

pt getIntersection(const halfPlane &lhs, const halfPlane &rhs) {
    double u = cross(lhs.a - rhs.a, rhs.dir);
    double v = cross(lhs.dir, rhs.dir);
    double t = -u / v;
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

        double sq = 0;

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
        printf("%.1f\n", sq);
    }
    return 0;
}