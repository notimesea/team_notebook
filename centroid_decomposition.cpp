struct graph {
    vector <vector <int> > g;
    vector <int> size;
    vector <int> par;

    graph(int n) {
        g.resize(n);
        size.resize(n);
        par.assign(n, -1);
    }

    void addEdge(int u, int v) {
        g[u].pb(v);
        g[v].pb(u);
    }

    void dfs(int v, int p) {
        size[v] = 1;
        for (int to : g[v]) {
            if (to == p) continue;
            if (par[to] != -1) continue;
            dfs(to, v);
            size[v] += size[to];
        }
    }

    void centroid(int v, int p, int n, int prev) {
        for (int to : g[v]) {
            if (to == p) continue;
            if (par[to] != -1) continue;
            if (2 * size[to] > n) {
                centroid(to, v, n, prev);
                return;
            }
        }
        //parent of each center is another center, for biggest center its value equal to -2
        par[v] = prev;
        
        //I'm center =)

        for (int to : g[v]) {
            if (par[to] == -1) {
                dfs(to, to);
                centroid(to, v, size[to], v);
            }
        }
    }

    void init() {
        dfs(0, 0);
        centroid(0, 0, g.size(), -2);
    }

};
