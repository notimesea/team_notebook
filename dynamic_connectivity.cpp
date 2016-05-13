struct Dsu {
    vector <int *> ptrs;
    vector <int> vals;
    vector <int> p;
    vector <int> sz;
    int components;

    void init(int n) {
        p.resize(n);
        sz.resize(n);
        ptrs.reserve(n);
        vals.reserve(n);
        components = n;
        for (int i = 0; i < n; ++i) {
            p[i] = i, sz[i] = 1;
        }
    }

    void set(int *x, int y) {
        ptrs.push_back(x);
        vals.push_back(*x);
        *x = y;
    }

    void revert(int version) {
        while (ptrs.size() > version) {
            *ptrs.back() = vals.back();
            ptrs.pop_back();
            vals.pop_back();
        }
    }

    inline int get(int v) {
        if (p[v] != v) {
            return get(p[v]);
            //set(p + v, get(p[v])); //compressed path heuristic
        }
        return p[v];
    }

    inline void unite(int a, int b) {
        a = get(a);
        b = get(b);
        if (a == b) {
            return;
        }
        set(&components, components - 1);
        if (sz[a] > sz[b]) {
            swap(a, b);
        }
        set(&p[a], b);
        set(&sz[b], sz[a] + sz[b]);
    }
    inline int getVersion() {
        return vals.size();
    }
};


Dsu dsu;
enum Type {
    GET, UPDATE
};

struct Query {
    Type type;
    int l, r; //time of birth and death of each edge
    int u, v; //edge
};

void processQueries(vector <Query> & queries, int left, int right, vector <int> & answer) {
    if (left == right) {
        int version = dsu.getVersion();
        for (const auto & query : queries) {
            if (query.type == UPDATE && query.l < left && query.r >= right) {
                dsu.unite(query.u, query.v);
            }
        }
        for (const auto & query : queries) {
            if (query.type == GET && query.l == left) {
                answer[query.u] = (dsu.components == 1); //this example shows that id = query.u and checks if graph is connected at time moment = query.l
            }
        }
        dsu.revert(version);
        return;
    }

    int version = dsu.getVersion();

    vector <Query> queriesToSend;

    for (const auto & query : queries) {
        if (query.type == UPDATE && query.l < left && query.r >= right) {
            dsu.unite(query.u, query.v);
        } else if ((query.type == GET && query.l >= left && query.r <= right) || (query.type == UPDATE && query.l < right && query.r >= left)) {
            queriesToSend.push_back(query);
        }
    }

    int mid = (left + right) / 2;

    processQueries(queriesToSend, left, mid, answer);
    processQueries(queriesToSend, mid + 1, right, answer);

    dsu.revert(version);
}
