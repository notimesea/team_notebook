#include <bits/stdc++.h>
#define forn(i, n) for(int i = 0; i < (int)(n); ++i)
#define pb push_back

using namespace std;
typedef vector<int> vi;
typedef long long i64;

i64 ext_gcd (i64 a, i64 b, i64 & x, i64 & y) {
    if (a == 0) { x = 0; y = 1; return b; }
    i64 x1, y1, d = ext_gcd(b%a, a, x1, y1);
    x = y1 - (b / a) * x1; y = x1;
    return d;
}

int inv(i64 x, int p) {
	int s = p - 2;
	i64 r = 1;
	while (s) {
		if (s & 1) 
			r = (r * x) % p;
		x = (x * x) % p;
		s >>= 1;
	}
	return r;
}

bool isPrime(int x) {
	if (x < 2) return 0;
	for (int i = 2; i * i <= x; i++) if (x % i == 0) return 0;
	return 1;
}

vector<int> kto_p;
vector<vi> kto_r;
void init_kto(int n = 100) {
    kto_p.clear();
    for(int x = 1e9; (int)kto_p.size() < n; ++x)
        if( isPrime(x) )
            kto_p.pb(x);
    kto_r.assign(n, vi(n));
    forn(i,n)forn(j,n)if(j>i)
        kto_r[i][j]=inv(kto_p[i], kto_p[j]);
}
int get_by_kto(vi x) {
    int r(0), m(1);
    forn(i, x.size()) {
        forn(j, i) {
            i64 cur = (i64)(x[i] - x[j]) * kto_r[j][i];
            cur %= kto_p[i]; if( cur < 0 ) cur += kto_p[i];
            x[i] = cur;                  
        }
        r = r + (m * x[i]);
        m = m * kto_p[i];
    } return r;
}

int powmod (int a, int b, int p) {
	int res = 1;
	while (b)
		if (b & 1)
			res = int (res * 1ll * a % p),  --b;
		else
			a = int (a * 1ll * a % p),  b >>= 1;
	return res;
}
 
int generator (int p) {
	vector<int> fact;
	int phi = p-1,  n = phi;
	for (int i=2; i*i<=n; ++i)
		if (n % i == 0) {
			fact.push_back (i);
			while (n % i == 0)
				n /= i;
		}
	if (n > 1)
		fact.push_back (n);
 
	for (int res=2; res<=p; ++res) {
		bool ok = true;
		for (size_t i=0; i<fact.size() && ok; ++i)
			ok &= powmod (res, phi / fact[i], p) != 1;
		if (ok)  return res;
	}
	return -1;
}

int gauss (vector < vector<double> > a, vector<double> & ans) {
	const int inf = 1e9; const double eps = 1e-8;
	int n = (int) a.size();
	int m = (int) a[0].size() - 1;
 
	vector<int> where (m, -1);
	for (int col=0, row=0; col<m && row<n; ++col) {
		int sel = row;
		for (int i=row; i<n; ++i)
			if (abs (a[i][col]) > abs (a[sel][col]))
				sel = i;
		if (abs (a[sel][col]) < eps)
			continue;
		for (int i=col; i<=m; ++i)
			swap (a[sel][i], a[row][i]);
		where[col] = row;
 
		for (int i=0; i<n; ++i)
			if (i != row) {
				double c = a[i][col] / a[row][col];
				for (int j=col; j<=m; ++j)
					a[i][j] -= a[row][j] * c;
			}
		++row;
	}
 
	ans.assign (m, 0);
	for (int i=0; i<m; ++i)
		if (where[i] != -1)
			ans[i] = a[where[i]][m] / a[where[i]][i];
	for (int i=0; i<n; ++i) {
		double sum = 0;
		for (int j=0; j<m; ++j)
			sum += ans[j] * a[i][j];
		if (abs (sum - a[i][m]) > eps)
			return 0;
	}
 
	for (int i=0; i<m; ++i)
		if (where[i] == -1)
			return inf;
	return 1;
}

int main() {
	return 0;
}
