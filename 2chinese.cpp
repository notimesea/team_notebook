#include <bits/stdc++.h>

#define all(v) (v).begin(), (v).end()
#define pb push_back
#define forn(i, n) for(int i = 0; i < (int)(n); ++i)
#define all(v) (v).begin(), (v).end()

using namespace std;
typedef vector<int> vi;

struct TEdge {
    int to, w, id;
 
    TEdge(int to = 0, int w = 0, int id = 0)
        : to(to)
        , w(w)
        , id(id)
    {
    }
};
 
typedef vector< vector<TEdge> > TGraph;
 
void dfs0(int v, TGraph &e, vi &vis, vi &tre) {
    vis[v] = 1;
    for (TEdge w: e[v]) {
        if (w.w || vis[w.to]) continue;
        tre.pb(w.id);
        dfs0(w.to, e, vis, tre);
    }
}
 
void dfs_ord(int v, TGraph &e, vi &vis, vi &ord) {
    if (vis[v]) return;
    vis[v] = 1;
    for (TEdge w: e[v]) dfs_ord(w.to, e, vis, ord);
    ord.pb(v);
}
 
void dfs_comp(int v, TGraph &e, vi &vis, vi &comp, int cc) {
    if (vis[v]) return;
    comp[v] = cc;
    vis[v] = 1;
    for (TEdge w: e[v]) dfs_comp(w.to, e, vis, comp, cc);
}
 
void dfs_cc0(int v, TGraph &e, vi &vis, vi &comp, vi &tre) {
    if (vis[v]) return;
    vis[v] = 1;
    for (TEdge w: e[v]) {
        if (vis[w.to] || comp[w.to] != comp[v]) continue;
        tre.pb(w.id);
        dfs_cc0(w.to, e, vis, comp, tre);
    }
}
 
TGraph rev(TGraph g) {
    int n = g.size();
    TGraph rg(n);
    forn(i, n) for (TEdge w: g[i]) rg[w.to].pb(TEdge(i, w.w, w.id));
    return rg;
}
 
vi condense(TGraph &e) {
    int n = e.size();
    vi vis(n), ord;
    forn(i, n) dfs_ord(i, e, vis, ord);
    reverse(all(ord));
    vis.assign(n, 0);
    TGraph re = rev(e);
    vi comp(n);
    int cc = 0;
    for (int v: ord) {
        if (vis[v]) continue;
        dfs_comp(v, re, vis, comp, cc++);
    }
    return comp;
}
 
vi mst(int v, TGraph e) {
    int n = e.size();
    vi pot(n, 1e9);
    forn(i, n) for (TEdge w: e[i]) pot[w.to] = min(pot[w.to], w.w);
    TGraph e0(n);
    forn(i, n) for (TEdge &w: e[i]) {
        w.w -= pot[w.to];
        if (!w.w) e0[i].pb(w);
    }
    vi vis(n);
    vi tre;
    dfs0(v, e0, vis, tre);
    if ((int)tre.size() == n - 1) return tre;
    vi comp = condense(e0);
    int cc = *max_element(all(comp)) + 1;
    TGraph ne(cc);
    forn(i, n) for (TEdge w: e[i]) {
        if (comp[i] != comp[w.to]) ne[comp[i]].pb(TEdge(comp[w.to], w.w, w.id));
    }
    vi cmst = mst(comp[v], ne);
    set<int> cid(all(cmst));
    vis.assign(n, 0);
    dfs_cc0(v, e0, vis, comp, cmst);
    forn(i, n) for (TEdge w: e[i]) {
        if (cid.count(w.id)) dfs_cc0(w.to, e0, vis, comp, cmst);
    }
    return cmst;
}

bool hasAnswer(int v, const TGraph& e) {
    vi used(e.size());
    queue<int> q; q.push(v); used[v] = 1;
    while (!q.empty()) {
		int v = q.front(); q.pop();
		for (auto ee: e[v])
			if (!used[ee.to]) {
				used[ee.to] = 1;
				q.push(ee.to);
			}	
	}
    if (accumulate(all(used), 0) != (int)e.size()) {
		return false;
	}
	return true;
}
// 11183 - Teen Girl Squad
// tested uva.onlinejudge.org/.../problem=2124

void testQuality() {
	forn (t, 100) {
		int n = 6 + rand() % 3;
		int m = n + rand() % 14;
		vi from(m), to(m), w(m);
		forn(i, n) from[i] = rand() % n, to[i] = rand() % n, w[i] = rand() % 10000;
		TGraph g(n); forn (i, m) g[from[i]].pb(TEdge(to[i], w[i], i));
		vi perm(m); forn (i, n - 1) perm[m - 1 - i] = 1;
		int ans = INT_MAX;
		do {
			int cw = 0; forn (i, m) cw += w[i] * perm[i];
			vector<vi> g(n); forn (i, m) if (perm[i]) g[from[i]].pb(to[i]);
			queue<int> q; vi used(n); used[0] = 1; q.push(0);
			while (!q.empty()) {
				int v = q.front(); q.pop();
				for (int to: g[v]) if (!used[to]) {used[to] = 1; q.push(to);}
			}
			if (accumulate(all(used), 0) == n) ans = min(ans, cw);
		} while (next_permutation(all(perm)));
		
		assert(hasAnswer(0, g) == (ans != INT_MAX));
		if (ans != INT_MAX) {
			vi edges = mst(0, g);
			int res = 0; for (int id: edges) res += w[id];
			assert(ans == res);
		}
	}
}

int main() {
	testQuality();
	return 0;
}
