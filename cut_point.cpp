#include <bits/stdc++.h>
#define forn(i, n) for(int i = 0; i < (int)(n); ++i)
#define pb push_back
using namespace std;
typedef vector<int> vi;
typedef pair<int, int> pii;
typedef vector<pii> vpi;

const int maxn = 100500;
int tin[maxn], timer = 0, fup[maxn];
struct Edge{
	int id, to;
};
vector<Edge> g[maxn];
bool used[maxn];
bool isCutPoint[maxn];
vpi st;
vector<vpi> components;
void dfs(int v, int p = -1){
    used[v] = true;
    fup[v] = tin[v] = timer++;
    int k = 0;
	forn(i, g[v].size()) {
	    if (g[v][i].to == p)
			continue;
		int to = g[v][i].to;
        if (used[to])
            fup[v] = min(fup[v], tin[to]);
        else
        {
            ++k;
            pii cedge(v, to);
            st.pb(cedge);
            dfs(to, v);
            vpi cur;
            if (fup[to] >= tin[v])
            {
                if (p != -1) 
                    isCutPoint[v] = true;
                while (st.back() != cedge)
                {
                    cur.pb(st.back());
                    st.pop_back();
                }
                st.pop_back();
                cur.pb(cedge);
                components.pb(cur);
            }

            fup[v] = min(fup[v], fup[to]);
        }
    }
    if (p == -1)
        isCutPoint[v] = (k > 1);
}

void testQuality() {
	int n = 3 + rand() % 20; components.clear();
	forn (i, n) g[i].clear(), used[i] = 0, isCutPoint[i] = 0; 
	int m = 1 + rand() % (n * n / 3);
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
				components.pb(st);
				st.clear();
			}
		}
	forn (b, n) {
		vi used(n);
		used[b] = 1;
		int cur = 0;
		forn (i, n)
			if (!used[i]) {
				cur++;
				queue<int> q;
				q.push(i); used[i] = 1;
				while (!q.empty()) {
					int v = q.front(); q.pop();
					for (auto e: g[v])
						if (!used[e.to]) {
							q.push(e.to);
							used[e.to] = 1;
						}
				}
			}
		assert((comps < cur) == isCutPoint[b]);
	}
}

int main() {
	forn (t, 100) testQuality();
	return 0;
}
