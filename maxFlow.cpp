//
// Created by Sergey Kiyan on 01.05.16.
//

#include <bits/stdc++.h>
using namespace std;
#define pb push_back

typedef long long D;
const D INF = LLONG_MAX;

struct Edge {
    int fr, to;
    D c, f;
    Edge() {}
    Edge(int fr, int to, D c) : fr(fr), to(to), c(c), f(0) {}
};

const int N = 5000;
struct Dinic {
    vector<Edge> e;
    vector<int> g[N];
    int q[N];
    int d[N];
    int ptr[N];
    int n, s, t;

    void addEdge(int fr, int to, D c) {
        g[fr].pb(int(e.size()));
        e.pb(Edge(fr, to, c));
        g[to].pb(int(e.size()));
        e.pb(Edge(to, fr, 0));
    }

    bool bfs() {
        memset(d, -1, n * sizeof(int));
        int q1 = 0, q2 = 0;
        q[q2++] = s;
        d[s] = 0;
        while (q1 < q2 && d[t] == -1) {
            int v = q[q1++];
            for (int id : g[v]) {
                if (e[id].f == e[id].c) continue;
                int to = e[id].to;

                if (d[to] == -1) {
                    d[to] = d[v] + 1;
                    q[q2++] = to;
                }
            }
        }
        return d[t] != -1;
    }

    D dfs(int v, D flow) {
        if (v == t || !flow) return flow;
        for (;ptr[v] < (int)g[v].size(); ++ptr[v]) {
            int id = g[v][ptr[v]];
            if (e[id].f == e[id].c) continue;
            int to = e[id].to;
            if (d[to] != d[v] + 1) continue;

            D pushed = dfs(to, min(flow, e[id].c - e[id].f));
            if (pushed) {
                e[id].f += pushed;
                e[id ^ 1].f -= pushed;
                return pushed;
            }
        }
        return 0;
    }

    D dinic(int S, int T) {
        s = S, t = T;
        D flow = 0;
        while (bfs()) {
            memset(ptr, 0, n * sizeof(int));
            while (int pushed = dfs(s, INF)) {
                flow += pushed;
            }
        }
        return flow;
    }

    void clearFlow() {
        for (Edge &ee : e) {
            ee.f = 0;
        }
    }
};

vector<vector<D> > flows;
vector<int> parent;
void homoryHu(Dinic &g) {
    int n = g.n;
    parent.assign(n, 0);
    flows.assign(n, vector<D>(n, INF));
    for (int i = 1; i < n; ++i) {
        g.clearFlow();
        D F = g.dinic(i, parent[i]);
        for (int j = i + 1; j < n; ++j) {
            if (g.d[j] != -1 && parent[j] == parent[i]) {
                parent[j] = i;
            }
        }
        flows[i][parent[i]] = flows[parent[i]][i] = F;
        for (int j = 0;j < i; ++j)
            flows[i][j] = flows[j][i] = min(F, flows[parent[i]][j]);
    }
    g.clearFlow();
}


// TEST

void testHomory() {
    const int V = 100;
    const int E = 1000;
    const int CAP = 1000;
    const int IT = 10;
    int it = IT;
    while (it--) {
        Dinic g;
        g.n = V;
        for (int i = 0; i < E; ++i) {
            int a = rand() % V;
            int b = rand() % V;
            while (a == b) {
                a = rand() % V;
                b = rand() % V;
            }
            D cap = rand() % CAP;
            g.addEdge(a, b, cap);
            g.addEdge(b, a, cap);
        }
        homoryHu(g);

        vector<vector<D> > res2(V, vector<D>(V, INF));
        for (int i = 0; i < V; ++i) {
            for (int j = 0; j < V; ++j) {
                if (i == j) continue;
                g.clearFlow();
                res2[i][j] = g.dinic(i, j);
            }
        }
        assert(flows == res2);
    }
}


int main() {
    //testHomory();

    //http://www.spoj.com/problems/FASTFLOW/

    int n, m;
    scanf("%d%d", &n, &m);
    Dinic g;
    g.n = n;
    for (int i = 0; i < m; ++i) {
        int a, b, c;
        scanf("%d%d%d", &a, &b, &c);
        --a, --b;
        g.addEdge(a, b, c);
        g.addEdge(b, a, c);
    }
    D ans = g.dinic(0, n - 1);
    cout << ans << "\n";

    return 0;
}
