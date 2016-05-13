#include <bits/stdc++.h>
using namespace std;
typedef long long ll;
typedef long double ld;

struct Line {
    ll m, b;

    ll eval(ll x) {
        return m * x + b;
    }
};

struct Hull {
    vector <Line> lines;
    size_t ptr = 0;
    static int comp(Line a, Line b, Line c) {
        ld val = (b.m - a.m) * 1.0L * (c.b - a.b) + (b.b - a.b) * 1.0L * (a.m - c.m);
        if (val < -1e18) return -1;
        if (val > 1e18) return 1;

        ll v = (b.m - a.m) * (c.b - a.b) + (b.b - a.b) * (a.m - c.m);
        if (v < 0) return -1;
        if (v > 0) return 1;

        return 0;
    }
    void add(const Line &line) {
        if (!lines.empty()) assert(line.m > lines.back().m);
        while (lines.size() >= 2 && comp(lines[lines.size() - 2], lines[lines.size() - 1], line) >= 0 /* <= 0 for minima */) {
            lines.pop_back();
        }
        lines.push_back(line);
    }
    ll eval(long long x) {
        if (lines.empty()) {
            return LLONG_MIN;
        }
        ptr = min(ptr, lines.size() - 1);
        while (ptr + 1 < lines.size() && lines[ptr].eval(x) < lines[ptr + 1].eval(x)) {
            ++ptr;
        }
        return lines[ptr].eval(x);
    }
};

void testStatic() {
    int it = 100;
    while (it--) {
        vector <Line> lines;
        int n = 1000;
        for (int i = 0; i < n; ++i) {
            lines.push_back(Line{i - 500, rand() % 9999 - 5000});
        }
        sort(lines.begin(), lines.end(), [&](const Line &l, const Line &r) {
            return l.m < r.m || (l.m == r.m && l.b < r.b);
        });


        Hull h;
        for (auto l : lines) {
            h.add(l);
        }

        int m = 10000;

        for (int i = -m; i < m; ++i) {
            ll best = LLONG_MIN;
            for (auto l : lines) {
                best = max(best, l.eval(i));
            }
            assert(best == h.eval(i));
        }
    }
}

#define Line LineD

const ll is_query = -(1LL<<62);
struct Line {
    ll m, b;
    mutable function<const Line*()> succ;
    bool operator<(const Line& rhs) const {
        if (rhs.b != is_query) return m < rhs.m;
        const Line* s = succ();
        if (!s) return 0;
        ll x = rhs.m;
        return b - s->b < (s->m - m) * x;
    }
};
struct HullDynamic : public multiset<Line> { // will maintain upper hull for maximum
    bool bad(iterator y) {
        auto z = next(y);
        if (y == begin()) {
            if (z == end()) return 0;
            return y->m == z->m && y->b <= z->b;
        }
        auto x = prev(y);
        if (z == end()) return y->m == x->m && y->b <= x->b;
        ll delta = (x->b - y->b)*(z->m - y->m) - (y->b - z->b)*(y->m - x->m);

        ld d = (x->b - y->b) * 1.0L * (z->m - y->m) - (y->b - z->b) * 1.0L * (y->m - x->m);
        if (d > 1e18) {
            return true;
        }
        if (d < -1e18) {
            return false;
        }
        return delta >= 0;
    }

    void insert_line(ll m, ll b) {
        auto y = insert({ m, b });
        y->succ = [=] { return next(y) == end() ? 0 : &*next(y); };
        if (bad(y)) { erase(y); return; }
        while (next(y) != end() && bad(next(y))) erase(next(y));
        while (y != begin() && bad(prev(y))) erase(prev(y));
    }
    ll eval(ll x) {
        if (empty()) return LLONG_MIN;
        auto l = *lower_bound((Line) { x, is_query });
        return l.m * x + l.b;
    }
};



int main() {
    testStatic();
    return 0;
}
