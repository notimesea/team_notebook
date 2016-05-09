#include <bits/stdc++.h>
#define forn(i, n) for(int i = 0; i < (int)(n); ++i)
#define pb push_back
using namespace std;
typedef vector<int> vi;
typedef pair<int, int> pii;
typedef vector<pii> vpi;

const int maxn = 100500;
int tin[maxn], timer = 0, fup[maxn], ans[maxn], numComps = 0;
struct Edge{
	int id, to;
};
vi st;
vector<Edge> g[maxn];
bool used[maxn];
bool isBridge[maxn];

void dfs(int v, int p = -1) {
	used[v] = true;
	fup[v] = tin[v] = timer++;
	st.push_back(v);
	forn (i, g[v].size()) {
		if (g[v][i].id == p)
			continue;
		int to = g[v][i].to;
		if (used[to]) 
			fup[v] = min(fup[v], tin[to]);
		else {
			dfs(to, g[v][i].id);
			if (fup[to] > tin[v]) { 
				numComps++;
				isBridge[g[v][i].id] = true;
				while (st.size()) {
					int x = st.back(); st.pop_back();
					ans[x] = numComps;
					if (x == to)
						break;
				}
			}
			fup[v] = min(fup[v], fup[to]);
		}
	} 
}

void testQuality() {
	int n = 3 + rand() % 100;
	forn (i, n) g[i].clear(), used[i] = 0; numComps = 0;
	int m = 1 + rand() % (n * n / 3);
	forn (i, m) isBridge[i] = 0;
	vpi edges(m);
	forn (i, m) {
		int x = rand() % n, y = rand() % n;
		edges[i] = make_pair(x, y);
		g[x].pb(Edge{i, y}); g[y].pb(Edge{i, x});
	}
	int comps = 0;
	forn (i, n) 
		if (!used[i]) {
			dfs(i);
			++comps;
			if (st.size()) {
				numComps++;
				for (int x: st) ans[x] = numComps;
				st.clear();
			}
		}
	forn (b, m) {
		vi used(n);
		int cur = 0;
		forn (i, n)
			if (!used[i]) {
				cur++;
				queue<int> q;
				q.push(i); used[i] = 1;
				while (!q.empty()) {
					int v = q.front(); q.pop();
					for (auto e: g[v])
						if (!used[e.to] && e.id != b) {
							q.push(e.to);
							used[e.to] = 1;
						}
				}
			}
		assert((comps < cur) == isBridge[b]);
	}
}

int main() {
	forn (t, 1000) testQuality();
	return 0;
}
