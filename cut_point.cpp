#include <bits/stdc++.h>

#define pb push_back
using namespace std;
typedef vector<int> vi;
typedef pair<int, int> pii;
typedef vector<pii> vpi;

const int maxn = 100500;
int n, m, tin[maxn], timer = 0, fup[maxn], ans[maxn], numComps = 0;
struct Edge{
	int id, tar;
};
vector<Edge> g[maxn];
bool used[maxn];
bool isCutPoint[maxn];
vpi st;
vector<vpi> components;
void dfs(int v, int p = -1){
    used[v] = true;
    timer++;
    fup[v] = tin[v] = timer;
    int kol = 0;
	for (int i = 0; i < (int)g[v].size(); i++) {
	    if (g[v][i].id == p)
			continue;
		int to = g[v][i].tar;
		
        if (to == p)
            continue;
        if (used[to])
            fup[v] = min(fup[v], tin[to]);
        else
        {
            ++kol;
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
        isCutPoint[v] = (kol > 1);
}

int main() {
	
	return 0;
}
