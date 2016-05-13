#include <cmath>
#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cassert>

#include <algorithm>
#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <functional>
#include <map>
#include <set>

using namespace std;

#define forn(i, n) for (int i = 0; i < (int)(n); i++)
#define forit(i, a) for (__typeof((a).begin()) i = (a).begin(); i != (a).end(); i++)
#define sz(a) (int)(a).size()
#define all(a) (a).begin(), (a).end()
#define zero(a) memset(a, 0, sizeof(a))
#define pb push_back
#define mp make_pair

typedef long long ll;
typedef vector <int> vi;
typedef pair <int, int> pii;

#include <cmath>
#include <string>
#include <iostream>

using namespace std;

typedef long double dbl;

const dbl eps = 1e-9;

template <class T> T inline sqr(T x) { return x * x; }
template <class T> T inline sgn(T x) { return x < 0 ? -1 : 1; }

struct pnt
{
    dbl x, y, z;

    pnt( dbl _x, dbl _y, dbl _z ) : x(_x), y(_y), z(_z) { }
    pnt() { }

    pnt operator + ( pnt p ) const { return pnt(x + p.x, y + p.y, z + p.z); }
    pnt operator - ( pnt p ) const { return pnt(x - p.x, y - p.y, z - p.z); }
    pnt operator * ( dbl k ) const { return pnt(x * k, y * k, z * k); }
    pnt operator / ( dbl k ) const { return pnt(x / k, y / k, z / k); }
    pnt operator - () const { return pnt(-x, -y, -z); }

    pnt operator * ( pnt p ) const
    {
        return pnt(y * p.z - z * p.y,
                   z * p.x - x * p.z,
                   x * p.y - y * p.x);
    }
    dbl operator ^ ( pnt p ) const { return x * p.x + y * p.y + z * p.z; }

    bool operator == ( pnt p ) const { return fabs(x - p.x) + fabs(y - p.y) + fabs(z - p.z) < eps; }
    bool operator != ( pnt p ) const { return fabs(x - p.x) + fabs(y - p.y) + fabs(z - p.z) > eps; }
    bool operator < ( pnt p ) const
    {
        if (fabs(x - p.x) > eps) return x < p.x;
        if (fabs(y - p.y) > eps) return y < p.y;
        return z < p.z;
    }

    void read() { cin >> x >> y >> z; }

    dbl d2() const { return x * x + y * y + z * z; }
    dbl d() const { return sqrt(d2()); }

    pnt norm() const { return *this / d(); }
};

const int maxn = (int)1e3 + 10;

int n;
pnt p[maxn];

pnt getV( int i, int j, int k )
{
    pnt v = p[j] - p[i];
    dbl t = ((p[k] - p[i]) ^ v) / (v ^ v);
    return p[k] - (p[i] + v * t);
}

set < pair<int, pii> > m;

void go( int i, int j, int k )
{
    int t[] = {i, j, k};
    int x = 0;
    forn(y, 3)
        if (t[y] < t[x])
            x = y;
    pair<int, pii> state = mp(t[x], mp(t[(x + 1) % 3], t[(x + 2) % 3]));

    if (m.count(state))
        return;
    m.insert(state);

    forn(t, 2)
    {
        pnt no = getV(i, j, k).norm();
        dbl opt = 2;
        int ml = -1;
        forn(l, n)
            if (l != i && l != j && l != k)
            {
                pnt v = getV(i, j, l);
                dbl val = v ^ no;
                val = sgn(val) * sqr(val) / v.d2();
                if (val < opt)
                    ml = l, opt = val;
            }
        assert(ml != -1);
        go(i, ml, j);
        int tmp = i; i = j; j = k; k = tmp;
    }
}

int main()
{
    assert(freopen("e1.in", "r", stdin));
    assert(freopen("e1.out", "w", stdout));

    int tn;
    scanf("%d", &tn);
    while (scanf("%d", &n) == 1 && tn--)
    {
        forn(i, n)
            p[i].read();

        int mi = 0;
        forn(i, n)
            if (p[i] < p[mi])
                mi = i;

        int mj = !mi;
#define F(j) ((p[j] - p[mi]).norm().x)
        forn(j, n)
            if (j != mi)
            if (F(j) < F(mj))
                mj = j;

        int mk = -1;
        forn(i, n)
            if (i != mi && i != mj)
            {
                pnt no = ((p[mi] - p[i]) * (p[mj] - p[i])).norm();
                dbl d = no ^ p[i];
                int bad = 0;
                forn(j, n)
                    if (j != mi && j != mj && j != i)
                    {
                        assert(fabs((no ^ p[j]) - d) > eps);
                        if ((no ^ p[j]) < d)
                            bad = 1, j = n;
                    }
                if (!bad)
                    mk = i, i = n;
            }
        fprintf(stderr, "%d %d %d\n", mi, mj, mk);
        assert(mk != -1);
        swap(mk, mj);

        m.clear();
        go(mi, mj, mk);
        printf("%d\n", sz(m));
        forit(it, m)
            printf("3 %d %d %d\n", it->first, it->second.first, it->second.second);
    }
    return 0;
}