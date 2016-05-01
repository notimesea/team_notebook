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
    D c, cost, f;
    Edge() {}
    Edge(int fr, int to, D c, D cost) : fr(fr), to(to), c(c), cost(cost), f(0) {}
};

const int N = 5000;
struct G {
    vector<Edge> e;
    vector<int> g[N];
    char u[N];
    D d[N];
    int q[N];
    int p[N];
    int n, s, t;

    void addEdge(int fr, int to, D c, D cost) {
        g[fr].pb(int(e.size()));
        e.pb(Edge(fr, to, c, cost));
        g[to].pb(int(e.size()));
        e.pb(Edge(to, fr, 0, -cost));
    }

    bool fb() {
        memset(p, -1, n * sizeof(int));
        memset(d, 0x3f, n * sizeof(D));
        d[s] = 0;
        u[s] = 1;
        int q1 = 0, q2 = 0;
        q[q2++] = s;
        while (q1 != q2) {
            int v = q[q1++];
            if (q1 == N) q1 = 0;
            u[v] = 0;
            for (int id : g[v]) {
                if (e[id].f == e[id].c) continue;
                int to = e[id].to;
                D len = e[id].cost;
                if (d[to] > d[v] + len) {
                    d[to] = d[v] + len;
                    p[to] = id;
                    if (!u[to]) {
                        q[q2++] = to;
                        if (q2 == N) q2 = 0;
                        u[to] = 1;
                    }
                }
            }
        }
        return p[t] != -1;
    }

    D push(D pushed = INF) {
        int v = t;
        while (v != s) {
            pushed = min(pushed, e[p[v]].c - e[p[v]].f);
            v = e[p[v]].fr;
        }
        v = t;
        while (v != s) {
            e[p[v]].f += pushed;
            e[p[v] ^ 1].f -= pushed;
            v = e[p[v]].fr;
        }
        return pushed;
    }

    pair<D, D> minCostFlow(int S, int T, D need) {
        s = S, t = T;
        D flow = 0;
        D cost = 0;
        while (fb() && need > 0) {
            D add = push(need);
            cost += add * d[t];
            flow += add;
        }
        return make_pair(flow, cost);
    }
};

// TEST

int main() {

    //http://www.spoj.com/problems/FASTFLOW/
    /*
    int n, m;
    scanf("%d%d", &n, &m);
    G g;
    g.n = n;
    for (int i = 0; i < m; ++i) {
        int a, b, c;
        scanf("%d%d%d", &a, &b, &c);
        --a, --b;
        g.addEdge(a, b, c, 1);
        g.addEdge(b, a, c, 1);
    }
    D ans = g.minCostFlow(0, n - 1, INF).first;
    cout << ans << "\n";
    */

    return 0;
}
