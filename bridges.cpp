#include <bits/stdc++.h>
using namespace std;
typedef vector<int> vi;
typedef pair<int, int> pii;
typedef vector<pii> vpi;

const int maxn = 100500;
int n, m, tin[maxn], timer = 0, fup[maxn], ans[maxn], numComps = 0;
struct Edge{
	int id, tar;
};
vi st;
vector<Edge> g[maxn];
bool used[maxn];
void dfs(int v, int p) {
	used[v] = true;
	timer++;
	fup[v] = tin[v] = timer;
	st.push_back(v);
	for (int i = 0; i < (int)g[v].size(); i++) {
		if (g[v][i].id == p)
			continue;
		int to = g[v][i].tar;
		if (used[to]) 
			fup[v] = min(fup[v], tin[to]);
		else {
			dfs(to, g[v][i].id);
			if (fup[to] > tin[v]) { 
				numComps++;
				//IS_BRIDGE(to, v);
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

int main() {
	return 0;
}
