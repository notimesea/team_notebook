#include <bits/stdc++.h>
#include "point.h"

using namespace std;

typedef ptl pt;

vector<pt> convex_hull_strict(vector<pt> a) {
    sort(a.begin(), a.end(), [&](const pt &l, const pt &r){
        if (l.y != r.y) return l.y < r.y;
        return l.x < r.x;
    });
    if (a.size() < 2) return a;
    int n = (int)a.size();
    pt p1 = a[0];
    pt p2 = a[n - 1];
    vector<pt> up(1, p1);
    vector<pt> down(1, p1);
    pt p21 = p2 - p1;
    for (int i = 1; i < n; ++i) {
        if (i == n - 1 || cross(p21, a[i] - p1) > 0) {
            while (up.size() > 1 && cross(up[up.size() - 1] - up[up.size() - 2], a[i] - up[up.size() - 1]) >= 0) {
                up.pop_back();
            }
            up.push_back(a[i]);
        }
        if (i == n - 1 || cross(p21, a[i] - p1) < 0) {
            while (down.size() > 1 &&
                   cross(down[down.size() - 1] - down[down.size() - 2], a[i] - down[down.size() - 1]) <= 0) {
                down.pop_back();
            }
            down.push_back(a[i]);
        }
    }
    a.clear();
    for (int i = 0; i < (int)up.size(); ++i) {
        a.push_back(up[i]);
    }
    for (int i = (int)down.size() - 2; i > 0; --i) {
        a.push_back(down[i]);
    }
    return a;
}


vector<pt> convex_hull_not_strict(vector<pt> a) {
    sort(a.begin(), a.end(), [&](const pt &l, const pt &r){
        if (l.x != r.x) return l.x < r.x;
        return l.y < r.y;
    });
    if (a.size() < 2) return a;
    int n = (int)a.size();
    pt p1 = a[0];
    pt p2 = a[n - 1];
    pt p21 = p2 - p1;
    bool is_on_line = true;
    for (int i = 0; i < (int)a.size(); ++i) {
        if (cross(p21, a[i] - p1) != 0) {
            is_on_line = false;
        }
    }
    if (is_on_line) return a;
    vector<pt> up(1, p1);
    vector<pt> down(1, p1);

    for (int i = 1; i < n; ++i) {
        if (i == n - 1 || cross(p21, a[i] - p1) >= 0) {
            while (up.size() > 1 && cross(up[up.size() - 1] - up[up.size() - 2], a[i] - up[up.size() - 1]) > 0) {
                up.pop_back();
            }
            up.push_back(a[i]);
        }
        if (i == n - 1 || cross(p21, a[i] - p1) <= 0) {
            while (down.size() > 1 &&
                   cross(down[down.size() - 1] - down[down.size() - 2], a[i] - down[down.size() - 1]) < 0) {
                down.pop_back();
            }
            down.push_back(a[i]);
        }
    }
    a.clear();
    for (int i = 0; i < (int)up.size(); ++i) {
        a.push_back(up[i]);
    }
    for (int i = (int)down.size() - 2; i > 0; --i) {
        a.push_back(down[i]);
    }
    return a;
}

