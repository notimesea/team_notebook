#include <vector>
#include <cassert>
#include <algorithm>
#include <iostream>
#include <numeric>

using namespace std;

typedef long long T;

pair <T, vector <int> > hungarian(const vector <vector <T> > &b) {
	T INF = (sizeof(T) == 4) ? INT_MAX / 2 : LLONG_MAX / 2;

	int n = b.size();
	int m = b[0].size();

	assert(n <= m);

	vector <vector <T> > a(n + 1, vector <T>(m + 1));

	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < m; ++j) {
			a[i + 1][j + 1] = b[i][j];
		}
	}

	vector<T> u(n+1), v(m+1);
	vector <int> p(m+1), way(m+1);

	for (int i=1; i<=n; ++i) {
		p[0] = i;
		int j0 = 0;
		vector<T> minv (m+1, INF);
		vector<char> used (m+1, false);
		do {
			used[j0] = true;
			int i0 = p[j0],  j1;
			T delta = INF;
			for (int j=1; j<=m; ++j)
				if (!used[j]) {
					T cur = a[i0][j]-u[i0]-v[j];
					if (cur < minv[j])
						minv[j] = cur,  way[j] = j0;
					if (minv[j] < delta)
						delta = minv[j],  j1 = j;
				}
			for (int j=0; j<=m; ++j)
				if (used[j])
					u[p[j]] += delta,  v[j] -= delta;
				else
					minv[j] -= delta;
			j0 = j1;
		} while (p[j0] != 0);
		do {
			int j1 = way[j0];
			p[j0] = p[j1];
			j0 = j1;
		} while (j0);
	}

	T cost = -v[0];

	vector<int> ans(n);
	for (int j=1; j<=m; ++j)
		if (p[j]) ans[p[j] - 1] = j - 1;

	T checkSum = 0;
	for (int i = 0; i < n; ++i) {
		checkSum += b[i][ans[i]];
	}
	assert(checkSum == cost);
	return make_pair(cost, ans);

}

void test() {
	int m = rand() % 10 + 1;
	int n = rand() % m + 1;
	cerr << n << " " << m << "\n";
	vector <vector <T> > a(n, vector <T>(m));
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < m; ++j) {
			a[i][j] = rand();
		}
	}

	vector <int> order(m);

	iota(order.begin(), order.end(), 0);

	T ans = LLONG_MAX;
	cerr << "I" << endl;
	do {
		T res = 0;
		for (int i = 0; i < n; ++i) {
			res += a[i][order[i]];
		}
		if (res < ans) {
			ans = res;
		}
	} while (next_permutation(order.begin(), order.end()));
	cerr << "O" << endl;
	cerr << "CHECK" << endl;
	assert(ans == hungarian(a).first);
	cerr << "OK" << endl;
}

int main() {
	int it = 10000;
	while (it--) {
		test();
	}
	return 0;
}