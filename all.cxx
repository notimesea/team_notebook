// 2 chinese
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

/////// AHO KORASIK
const int K = 26;
struct vertex {
	int next[K];
	bool leaf;
	int p;
	char pch;
	int link;
	int go[K];
};
 
vertex t[NMAX+1];
int sz;
 
void init() {
	t[0].p = t[0].link = -1;
	memset (t[0].next, 255, sizeof t[0].next);
	memset (t[0].go, 255, sizeof t[0].go);
	sz = 1;
}
 
void add_string (const string & s) {
	int v = 0;
	for (size_t i=0; i<s.length(); ++i) {
		char c = s[i]-'a'; //be careful here, it might be '0' instead of 'a'
		if (t[v].next[c] == -1) {
			memset (t[sz].next, 255, sizeof t[sz].next);
			memset (t[sz].go, 255, sizeof t[sz].go);
			t[sz].link = -1;
			t[sz].p = v;
			t[sz].pch = c;
			t[v].next[c] = sz++;
		}
		v = t[v].next[c];
	}
	t[v].leaf = true;
}
 
int go (int v, char c);
 
int get_link (int v) {
	if (t[v].link == -1)
		if (v == 0 || t[v].p == 0)
			t[v].link = 0;
		else
			t[v].link = go (get_link (t[v].p), t[v].pch);
	return t[v].link;
}
 
int go (int v, char c) {
	if (t[v].go[c] == -1)
		if (t[v].next[c] != -1)
			t[v].go[c] = t[v].next[c];
		else
			t[v].go[c] = v==0 ? 0 : go (get_link (v), c);
	return t[v].go[c];
}

//BRIDGES
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

// CUT POINTS
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

// BRIDGES ONLINE
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

// CARTESIAN
struct item {
        i64 key;
        int prior;
        int q;
        i64 sum;
        item * l, * r, * par;
        item() { }
        item (i64 key, int prior) : key(key), prior(prior), q(1), sum(key), l(NULL), r(NULL), par(NULL) { }
};
typedef item * pitem;

i64 sum(pitem t) {
	if (t) return t->sum;
	else return 0;
}

i64 q(pitem t) {
	if (t) return t->q;
	else return 0;
}

void upd(pitem t) {
	if (t)
		t->sum = t->key + sum(t->l) + sum(t->r), t->q = 1 + q(t->l) + q(t->r);
}

void setPar(pitem t, pitem p) {
	if (t) t->par = p;
}


void split (pitem t, i64 key, pitem & l, pitem & r) {
        if (!t) l = r = NULL;
        else if (key < t->key)
            split (t->l, key, l, t->l),  r = t;
        else
			split (t->r, key, t->r, r),  l = t;
		upd(t);
}

void merge (pitem & t, pitem l, pitem r) {
        if (!l || !r)
                t = l ? l : r;
        else if (l->prior > r->prior)
                merge (l->r, l->r, r),  t = l;
        else
			merge (r->l, l, r->l),  t = r;
		upd(t);
}

pitem merge (pitem l, pitem r) {
	if (!l || !r)
		return l ? l : r;
	else if (l->prior > r->prior) {
		l->r = merge (l->r, r);
		setPar(l->r, l);
		upd(l);
		return l;
	}
	else {
		r->l = merge (l, r->l);
		setPar(r->l, r);
		upd(r);
		return r;
	}
}

pair<pitem, pitem> split(pitem root, i64 key)
{
	if (!root)
		return make_pair((pitem)NULL, (pitem)NULL);
	
	if (key < root->key)
	{
		pair<pitem, pitem> splitted = split(root->l, key);
		root->l = splitted.second;
		setPar(splitted.second, root);
		setPar(splitted.first, NULL);
		upd(root);
		return make_pair(splitted.first, root);
	}
	else
	{
		pair<pitem, pitem> splitted = split(root->r, key);
		root->r = splitted.first;
		setPar(splitted.first, root);
		setPar(splitted.second, NULL);
		upd(root);
		return make_pair(root, splitted.second);
	}
}

void split2(pitem root, i64 key, pitem& l, pitem& r) {
	auto w = split(root, key);
	l = w.fi, r = w.se; setPar(l, NULL); setPar(r, NULL);
}

pitem par(pitem t) {
	if (t) return t->par;
	else return NULL;
}

void split_implicit(pitem t, pitem & l, pitem & r, int key) {
	if (!t)
		return void( l = r = 0 );
	if (q(t->l) <= key)
		split_implicit(t->l, l, t->l, key),  r = t;
	else
		split_implicit(t->r, t->r, r, key - 1 - q(t->l)),  l = t;
	upd(t);
}

//CENTROID DECOMPOSITION
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

// GEOMETRY
void intersect_with_line(double r, double a, double b, double c) {


    double x0 = -a*c/(a*a+b*b),  y0 = -b*c/(a*a+b*b);
    if (c*c > r*r*(a*a+b*b)+EPS)
        puts ("no points");
    else if (abs (c*c - r*r*(a*a+b*b)) < EPS) {
        puts ("1 point");
        cout << x0 << ' ' << y0 << '\n';
    }
    else {
        double d = r*r - c*c/(a*a+b*b);
        double mult = sqrt (d / (a*a+b*b));
        double ax,ay,bx,by;
        ax = x0 + b * mult;
        bx = x0 - b * mult;
        ay = y0 - a * mult;
        by = y0 + a * mult;
        puts ("2 points");
        cout << ax << ' ' << ay << '\n' << bx << ' ' << by << '\n';
    }
}

void intersect_two_circles() {
    double r1, x, y, r2;
    intersect_with_line(r1, -2x, -2y, x*x + y*y + r1*r1 - r2*r2);
}struct pt {
    int x, y, id;
};

inline bool cmp_x (const pt & a, const pt & b) {
    return a.x < b.x || a.x == b.x && a.y < b.y;
}

inline bool cmp_y (const pt & a, const pt & b) {
    return a.y < b.y;
}

pt a[MAXN];

double mindist;
int ansa, ansb;

inline void upd_ans (const pt & a, const pt & b) {
    double dist = sqrt ((a.x-b.x)*(a.x-b.x) + (a.y-b.y)*(a.y-b.y) + .0);
    if (dist < mindist)
        mindist = dist,  ansa = a.id,  ansb = b.id;
}

void rec (int l, int r) {
    if (r - l <= 3) {
        for (int i=l; i<=r; ++i)
            for (int j=i+1; j<=r; ++j)
                upd_ans (a[i], a[j]);
        sort (a+l, a+r+1, &cmp_y);
        return;
    }

    int m = (l + r) >> 1;
    int midx = a[m].x;
    rec (l, m),  rec (m+1, r);
    static pt t[MAXN];
    merge (a+l, a+m+1, a+m+1, a+r+1, t, &cmp_y);
    copy (t, t+r-l+1, a+l);

    int tsz = 0;
    for (int i=l; i<=r; ++i)
        if (abs (a[i].x - midx) < mindist) {
            for (int j=tsz-1; j>=0 && a[i].y - t[j].y < mindist; --j)
                upd_ans (a[i], t[j]);
            t[tsz++] = a[i];
        }
}


int main() {
    sort (a, a+n, &cmp_x);
    mindist = 1E20;
    rec (0, n-1);
}

// KOPELIOVICH 3D CONVEX HULL
#define forn(i, n) for (int i = 0; i < (int)(n); i++)
#define forit(i, a) for (__typeof((a).begin()) i = (a).begin(); i != (a).end(); i++)
#define sz(a) (int)(a).size()
#define all(a) (a).begin(), (a).end()
#define zero(a) memset(a, 0, sizeof(a))
#define pb push_back
#define mp make_pair

typedef long long ll;
typedef vector <int> vi;
typedef pair <int, int> pii;

using namespace std;

typedef long double dbl;

const dbl eps = 1e-9;

template <class T> T inline sqr(T x) { return x * x; }
template <class T> T inline sgn(T x) { return x < 0 ? -1 : 1; }

struct pnt
{
    dbl x, y, z;

    pnt( dbl _x, dbl _y, dbl _z ) : x(_x), y(_y), z(_z) { }
    pnt() { }

    pnt operator + ( pnt p ) const { return pnt(x + p.x, y + p.y, z + p.z); }
    pnt operator - ( pnt p ) const { return pnt(x - p.x, y - p.y, z - p.z); }
    pnt operator * ( dbl k ) const { return pnt(x * k, y * k, z * k); }
    pnt operator / ( dbl k ) const { return pnt(x / k, y / k, z / k); }
    pnt operator - () const { return pnt(-x, -y, -z); }

    pnt operator * ( pnt p ) const
    {
        return pnt(y * p.z - z * p.y,
                   z * p.x - x * p.z,
                   x * p.y - y * p.x);
    }
    dbl operator ^ ( pnt p ) const { return x * p.x + y * p.y + z * p.z; }

    bool operator == ( pnt p ) const { return fabs(x - p.x) + fabs(y - p.y) + fabs(z - p.z) < eps; }
    bool operator != ( pnt p ) const { return fabs(x - p.x) + fabs(y - p.y) + fabs(z - p.z) > eps; }
    bool operator < ( pnt p ) const
    {
        if (fabs(x - p.x) > eps) return x < p.x;
        if (fabs(y - p.y) > eps) return y < p.y;
        return z < p.z;
    }

    void read() { cin >> x >> y >> z; }

    dbl d2() const { return x * x + y * y + z * z; }
    dbl d() const { return sqrt(d2()); }

    pnt norm() const { return *this / d(); }
};

const int maxn = (int)1e3 + 10;

int n;
pnt p[maxn];

pnt getV( int i, int j, int k )
{
    pnt v = p[j] - p[i];
    dbl t = ((p[k] - p[i]) ^ v) / (v ^ v);
    return p[k] - (p[i] + v * t);
}

set < pair<int, pii> > m;

void go( int i, int j, int k )
{
    int t[] = {i, j, k};
    int x = 0;
    forn(y, 3)
        if (t[y] < t[x])
            x = y;
    pair<int, pii> state = mp(t[x], mp(t[(x + 1) % 3], t[(x + 2) % 3]));

    if (m.count(state))
        return;
    m.insert(state);

    forn(t, 2)
    {
        pnt no = getV(i, j, k).norm();
        dbl opt = 2;
        int ml = -1;
        forn(l, n)
            if (l != i && l != j && l != k)
            {
                pnt v = getV(i, j, l);
                dbl val = v ^ no;
                val = sgn(val) * sqr(val) / v.d2();
                if (val < opt)
                    ml = l, opt = val;
            }
        assert(ml != -1);
        go(i, ml, j);
        int tmp = i; i = j; j = k; k = tmp;
    }
}

int main()
{
    assert(freopen("e1.in", "r", stdin));
    assert(freopen("e1.out", "w", stdout));

    int tn;
    scanf("%d", &tn);
    while (scanf("%d", &n) == 1 && tn--)
    {
        forn(i, n)
            p[i].read();

        int mi = 0;
        forn(i, n)
            if (p[i] < p[mi])
                mi = i;

        int mj = !mi;
#define F(j) ((p[j] - p[mi]).norm().x)
        forn(j, n)
            if (j != mi)
            if (F(j) < F(mj))
                mj = j;

        int mk = -1;
        forn(i, n)
            if (i != mi && i != mj)
            {
                pnt no = ((p[mi] - p[i]) * (p[mj] - p[i])).norm();
                dbl d = no ^ p[i];
                int bad = 0;
                forn(j, n)
                    if (j != mi && j != mj && j != i)
                    {
                        assert(fabs((no ^ p[j]) - d) > eps);
                        if ((no ^ p[j]) < d)
                            bad = 1, j = n;
                    }
                if (!bad)
                    mk = i, i = n;
            }
        fprintf(stderr, "%d %d %d\n", mi, mj, mk);
        assert(mk != -1);
        swap(mk, mj);

        m.clear();
        go(mi, mj, mk);
        printf("%d\n", sz(m));
        forit(it, m)
            printf("3 %d %d %d\n", it->first, it->second.first, it->second.second);
    }
    return 0;
}

// CONVEX HULL
typedef ptl pt;

vector<pt> convex_hull_strict(vector<pt> a) {
    sort(a.begin(), a.end(), [&](const pt &l, const pt &r){
        if (l.y != r.y) return l.y < r.y;
        return l.x < r.x;
    });
    if (a.size() < 2) return a;
    int n = (int)a.size();
    pt p1 = a[0];
    pt p2 = a[n - 1];
    vector<pt> up(1, p1);
    vector<pt> down(1, p1);
    pt p21 = p2 - p1;
    for (int i = 1; i < n; ++i) {
        if (i == n - 1 || cross(p21, a[i] - p1) > 0) {
            while (up.size() > 1 && cross(up[up.size() - 1] - up[up.size() - 2], a[i] - up[up.size() - 1]) >= 0) {
                up.pop_back();
            }
            up.push_back(a[i]);
        }
        if (i == n - 1 || cross(p21, a[i] - p1) < 0) {
            while (down.size() > 1 &&
                   cross(down[down.size() - 1] - down[down.size() - 2], a[i] - down[down.size() - 1]) <= 0) {
                down.pop_back();
            }
            down.push_back(a[i]);
        }
    }
    a.clear();
    for (int i = 0; i < (int)up.size(); ++i) {
        a.push_back(up[i]);
    }
    for (int i = (int)down.size() - 2; i > 0; --i) {
        a.push_back(down[i]);
    }
    return a;
}


vector<pt> convex_hull_not_strict(vector<pt> a) {
    sort(a.begin(), a.end(), [&](const pt &l, const pt &r){
        if (l.x != r.x) return l.x < r.x;
        return l.y < r.y;
    });
    if (a.size() < 2) return a;
    int n = (int)a.size();
    pt p1 = a[0];
    pt p2 = a[n - 1];
    pt p21 = p2 - p1;
    bool is_on_line = true;
    for (int i = 0; i < (int)a.size(); ++i) {
        if (cross(p21, a[i] - p1) != 0) {
            is_on_line = false;
        }
    }
    if (is_on_line) return a;
    vector<pt> up(1, p1);
    vector<pt> down(1, p1);

    for (int i = 1; i < n; ++i) {
        if (i == n - 1 || cross(p21, a[i] - p1) >= 0) {
            while (up.size() > 1 && cross(up[up.size() - 1] - up[up.size() - 2], a[i] - up[up.size() - 1]) > 0) {
                up.pop_back();
            }
            up.push_back(a[i]);
        }
        if (i == n - 1 || cross(p21, a[i] - p1) <= 0) {
            while (down.size() > 1 &&
                   cross(down[down.size() - 1] - down[down.size() - 2], a[i] - down[down.size() - 1]) < 0) {
                down.pop_back();
            }
            down.push_back(a[i]);
        }
    }
    a.clear();
    for (int i = 0; i < (int)up.size(); ++i) {
        a.push_back(up[i]);
    }
    for (int i = (int)down.size() - 2; i > 0; --i) {
        a.push_back(down[i]);
    }
    return a;
}

// CONVEX HULL TRICK

struct Line {
    ll m, b;

    ll eval(ll x) {
        return m * x + b;
    }
};

struct Hull {
    vector <Line> lines;
    size_t ptr = 0;
    static int comp(Line a, Line b, Line c) {
        ld val = (b.m - a.m) * 1.0L * (c.b - a.b) + (b.b - a.b) * 1.0L * (a.m - c.m);
        if (val < -1e18) return -1;
        if (val > 1e18) return 1;

        ll v = (b.m - a.m) * (c.b - a.b) + (b.b - a.b) * (a.m - c.m);
        if (v < 0) return -1;
        if (v > 0) return 1;

        return 0;
    }
    void add(const Line &line) {
        if (!lines.empty()) assert(line.m > lines.back().m);
        while (lines.size() >= 2 && comp(lines[lines.size() - 2], lines[lines.size() - 1], line) >= 0 /* <= 0 for minima */) {
            lines.pop_back();
        }
        lines.push_back(line);
    }
    ll eval(long long x) {
        if (lines.empty()) {
            return LLONG_MIN;
        }
        ptr = min(ptr, lines.size() - 1);
        while (ptr + 1 < lines.size() && lines[ptr].eval(x) < lines[ptr + 1].eval(x)) {
            ++ptr;
        }
        return lines[ptr].eval(x);
    }
};

#define Line LineD

const ll is_query = -(1LL<<62);
struct Line {
    ll m, b;
    mutable function<const Line*()> succ;
    bool operator<(const Line& rhs) const {
        if (rhs.b != is_query) return m < rhs.m;
        const Line* s = succ();
        if (!s) return 0;
        ll x = rhs.m;
        return b - s->b < (s->m - m) * x;
    }
};
struct HullDynamic : public multiset<Line> { // will maintain upper hull for maximum
    bool bad(iterator y) {
        auto z = next(y);
        if (y == begin()) {
            if (z == end()) return 0;
            return y->m == z->m && y->b <= z->b;
        }
        auto x = prev(y);
        if (z == end()) return y->m == x->m && y->b <= x->b;
        ll delta = (x->b - y->b)*(z->m - y->m) - (y->b - z->b)*(y->m - x->m);

        ld d = (x->b - y->b) * 1.0L * (z->m - y->m) - (y->b - z->b) * 1.0L * (y->m - x->m);
        if (d > 1e18) {
            return true;
        }
        if (d < -1e18) {
            return false;
        }
        return delta >= 0;
    }

    void insert_line(ll m, ll b) {
        auto y = insert({ m, b });
        y->succ = [=] { return next(y) == end() ? 0 : &*next(y); };
        if (bad(y)) { erase(y); return; }
        while (next(y) != end() && bad(next(y))) erase(next(y));
        while (y != begin() && bad(prev(y))) erase(prev(y));
    }
    ll eval(ll x) {
        if (empty()) return LLONG_MIN;
        auto l = *lower_bound((Line) { x, is_query });
        return l.m * x + l.b;
    }
};


// FFT
const int N = 1 << 21;
const D PI = acosl(-1.0);
struct C {
    D x, y;
    C() {}
    C(D x, D y) : x(x), y(y) {}
    C operator + (const C &r) const {
        return C(x + r.x, y + r.y);
    }
    C operator - (const C &r) const {
        return C(x - r.x, y - r.y);
    }
    C operator * (const C &r) const {
        return C(x * r.x - y * r.y, y * r.x + x * r.y);
    }
    C &operator /= (const D &r) {
        x /= r;
        y /= r;
        return *this;
    }
};

C pw[N];
C ipw[N];
C cpw[N];

void initFFT() {
    for (int i = 0; i < N; ++i) {
        D ang = i * 2 * PI / N;
        D co = cosl(ang);
        D si = sinl(ang);
        pw[i] = C(co, si);
        ipw[i] = C(co, -si);
    }
}

void fft(C *a, int n, bool inv) {
    for (int i = 1, j = 0; i < n; ++i) {
        int bit = n >> 1;
        while (j >= bit) {
            j -= bit; bit >>= 1;
        }
        j += bit;
        if (i < j) swap(a[i], a[j]);
    }
    for (int len = 2, shift = N >> 1; len <= n; len <<= 1, shift >>= 1) {
        int len2 = len >> 1;
        for (int j = 0; j < len2; ++j) {
            cpw[j] = inv ? ipw[j * shift] : pw[j * shift];
        }
        for (int i = 0; i < n; i += len) {
            for (int j = 0; j < len2; ++j) {
                C u = a[i + j];
                C v = a[i + j + len2] * cpw[j];
                a[i + j] = u + v;
                a[i + j + len2] = u - v;
            }
        }
    }
    if (inv) {
        for (int i = 0; i < n; ++i) {
            a[i] /= n;
        }
    }
}

C A[N];
C B[N];

vector<ll> mulFFT(const vector<ll> &a, const vector<ll> &b) {
    int n = 1;
    while (n < (int)(a.size() + b.size()) - 1) {
        n <<= 1;
    }
    assert(n <= N);
    memset(A, 0, n * sizeof(C));
    memset(B, 0, n * sizeof(C));
    for (size_t i = 0; i < a.size(); ++i) {
        A[i].x = a[i];
    }
    for (size_t i = 0; i < b.size(); ++i) {
        B[i].x = b[i];
    }
    fft(A, n, false);
    fft(B, n, false);
    for (int i = 0; i < n; ++i) {
        A[i] = A[i] * B[i];
    }
    fft(A, n, true);
    vector<ll> c(n);
    for (int i = 0; i < n; ++i) {
        c[i] = llround(A[i].x);
    }
    while (!c.empty() && c.back() == 0) {
        c.pop_back();
    }
    return c;
}


// HEAVY-LIGHT DECOMPOSITION
struct Hld {
    vector <int> path;
    vector <vector <int> > g;
    vector <int> position;
    vector <int> parent;
    vector <int> size;
    vector <vector <int> > onPath;
    vector <int> rootId;
    vector <int> tin;
    vector <int> tout;
    int pathCount;
    int curPosition;
    int timer;

    Hld(int n) {
        path.resize(n);
        g.resize(n);
        position.resize(n);
        parent.resize(n);
        size.resize(n);
        tin.resize(n);
        tout.resize(n);
        pathCount = 0;
        curPosition = 0;
        timer = 0;
    }

    void addEdge(int u, int v) {
        g[u].push_back(v);
        g[v].push_back(u);
    }

    void dfs(int v, int p) {
        tin[v] = timer++;
        size[v] = 1;

        for (int to : g[v]) {
            if (to == p) {
                continue;
            }
            dfs(to, v);
            size[v] += size[to];
        }
        tout[v] = timer++;
    }

    void build(int v, int p, int pathId) {
        onPath[pathId].push_back(v);
        path[v] = pathId;
        position[v] = curPosition++;

        int best = -1;
        for (int to : g[v]) {
            if (to == p) {
                continue;
            }
            if (best == -1 || size[best] < size[to]) {
                best = to;
            }
        }

        if (best != -1) {
            parent[best] = v;
            build(best, v, pathId);
        }

        for (int to : g[v]) {
            if (to == p || to == best) {
                continue;
            }
            onPath.push_back(vector <int>());
            rootId.push_back(to);
            parent[to] = v;
            build(to, v, pathCount++);
        }
    }

    void build() {
        dfs(0, 0);
        onPath.push_back(vector <int>());
        rootId.push_back(0);
        parent[0] = -1;
        build(0, 0, pathCount++);
    }

    bool upper(int a, int b) {
        return tin[a] <= tin[b] && tout[a] >= tout[b];
    }

    int lca(int a, int b) {
        while (!upper(rootId[path[a]], b)) {
            a = parent[rootId[path[a]]];
        }
        while (!upper(rootId[path[b]], a)) {
            b = parent[rootId[path[b]]];
        }
        if (upper(a, b)) {
            return a;
        } else {
            return b;
        }
    }

    vector <pair <int, int> > query(int a, int b = 0) {
        vector <pair <int, int> > segs;

        if (!upper(b, a)) {
            return segs;
        }

        while (path[a] != path[b]) {
            segs.push_back(make_pair(position[rootId[path[a]]], position[a]));
            a = parent[rootId[path[a]]];
        }
        segs.push_back(make_pair(position[b], position[a]));
        return segs;
    }
};

// HUNGARIAN
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

// MAXIMAL INCIRCLE
}struct pt {
    double x, y;
    pt()  { }
    pt (double x, double y) : x(x), y(y)  { }
    pt operator- (const pt & p) const {
        return pt (x-p.x, y-p.y);
    }
};

double dist (const pt & a, const pt & b) {
    return sqrt ((a.x-b.x)*(a.x-b.x) + (a.y-b.y)*(a.y-b.y));
}

double get_ang (const pt & a, const pt & b) {
    double ang = abs (atan2 (a.y, a.x) - atan2 (b.y, b.x));
    return min (ang, 2*PI-ang);
}

struct line {
    double a, b, c;
    line (const pt & p, const pt & q) {
        a = p.y - q.y;
        b = q.x - p.x;
        c = - a * p.x - b * p.y;
        double z = sqrt (a*a + b*b);
        a/=z, b/=z, c/=z;
    }
};

double det (double a, double b, double c, double d) {
    return a * d - b * c;
}

pt intersect (const line & n, const line & m) {
    double zn = det (n.a, n.b, m.a, m.b);
    return pt (
            - det (n.c, n.b, m.c, m.b) / zn,
            - det (n.a, n.c, m.a, m.c) / zn
    );
}

bool parallel (const line & n, const line & m) {
    return abs (det (n.a, n.b, m.a, m.b)) < EPS;
}

double get_h (const pt & p1, const pt & p2,
              const pt & l1, const pt & l2, const pt & r1, const pt & r2)
{
    pt q1 = intersect (line (p1, p2), line (l1, l2));
    pt q2 = intersect (line (p1, p2), line (r1, r2));
    double l = dist (q1, q2);
    double alpha = get_ang (l2 - l1, p2 - p1) / 2;
    double beta = get_ang (r2 - r1, p1 - p2) / 2;
    return l * sin(alpha) * sin(beta) / sin(alpha+beta);
}

struct cmp {
    bool operator() (const pair<double,int> & a, const pair<double,int> & b) const {
        if (abs (a.first - b.first) > EPS)
            return a.first < b.first;
        return a.second < b.second;
    }
};

int main() {
    int n;
    vector<pt> p;

    vector<int> next (n), prev (n);
    for (int i=0; i<n; ++i) {
        next[i] = (i + 1) % n;
        prev[i] = (i - 1 + n) % n;
    }

    set < pair<double,int>, cmp > q;
    vector<double> h (n);
    for (int i=0; i<n; ++i) {
        h[i] = get_h (
                p[i], p[next[i]],
                p[i], p[prev[i]],
                p[next[i]], p[next[next[i]]]
        );
        q.insert (make_pair (h[i], i));
    }

    double last_time;
    while (q.size() > 2) {
        last_time = q.begin()->first;
        int i = q.begin()->second;
        q.erase (q.begin());

        next[prev[i]] = next[i];
        prev[next[i]] = prev[i];
        int nxt = next[i],   nxt1 = (nxt+1)%n,
                prv = prev[i],   prv1 = (prv+1)%n;
        if (parallel (line (p[nxt], p[nxt1]), line (p[prv], p[prv1])))
            break;

        q.erase (make_pair (h[nxt], nxt));
        q.erase (make_pair (h[prv], prv));

        h[nxt] = get_h (
                p[nxt], p[nxt1],
                p[prv1], p[prv],
                p[next[nxt]], p[(next[nxt]+1)%n]
        );
        h[prv] = get_h (
                p[prv], p[prv1],
                p[(prev[prv]+1)%n], p[prev[prv]],
                p[nxt], p[nxt1]
        );

        q.insert (make_pair (h[nxt], nxt));
        q.insert (make_pair (h[prv], prv));
    }

    cout << last_time << endl;
}

//KARATSUBA
template <class T>
struct karatsuba {
    static const int N = 1 << 22;
    T buffer[N];
    T* pos = buffer;
    T* allocate(const size_t & size) {
        memset(pos, 0, sizeof(T) * size);
        T * res = pos;
        pos += size;
        return res;
    }
    T* copyOf(const size_t & size, const T * src) {
        memcpy(pos, src, sizeof(T) * size);
        T * res = pos;
        pos += size;
        return res;
    }
    void clear(const size_t & size) {
        pos -= size;
    }
    T* multiply(const T* a, const T* b, size_t n) {
        assert(!(n & (n - 1)));
        if (n <= 32) {
            T * cur = allocate(2 * n);
            for (int i = 0; i < n; ++i) {
                for (int j = 0; j < n; ++j) {
                    cur[i + j] += a[i] * b[j];
                }
            }
            return cur;
        }

        const size_t n2 = n / 2;
        T* res = allocate(2 * n);

        T* l = multiply(a, b, n2);
        T* r = multiply(a + n2, b + n2, n2);

        T* la = copyOf(n2, a);
        T* lb = copyOf(n2, b);
        for (int i = 0; i < n2; ++i) la[i] += a[i + n2];
        for (int i = 0; i < n2; ++i) lb[i] += b[i + n2];

        T* m = multiply(la, lb, n2);
        for (int i = 0; i < n; ++i) res[i] += l[i];
        for (int i = 0; i < n; ++i) res[i + n] += r[i];
        for (int i = 0; i < n; ++i) res[i + n2] += m[i] - l[i] - r[i];
        clear(4 * n);
        return res;
    }
};


//MANACHER
string prepare(const string &s) {
    string t(2 * s.length() + 3, 'a');
    t[0] = '^';
    t[1] = '#';
    for (int i = 0; i < (int)s.length(); ++i) {
        t[2 * i + 2] = s[i];
        t[2 * i + 3] = '#';
    }
    t[2 * s.length() + 2] = '&';
    return t;
}

vector<int> manacher(const string &s) {
    assert(!s.empty());
    string t = prepare(s);
    vector<int> p(t.length(), 0);
    int c = 0, r = 0;
    for (int i = 1; i + 1 < (int)t.length(); ++i) {
        int j = 2 * c - i;
        p[i] = (r > i) ? min(r - i, p[j]) : 0;
        while (t[i - p[i] - 1] == t[i + p[i] + 1]) {
            ++p[i];
        }
        if (i + p[i] > r) {
            c = i;
            r = i + p[i];
        }
    }
    vector<int> ret(2 * s.length() - 1);
    for (int i = 0; i < int(2 * s.length() - 1); ++i) {
        ret[i] = p[i + 2];
    }
    return ret;
}


int main() {
    string s;
    cin >> s;
    auto res = manacher(s);
    int ans = 0;
    int id = 0;
    for (int i = 0; i < (int)res.size(); ++i) {
        if (res[i] > ans) {
            ans = res[i];
            id = i;
        }
    }
    string t = s.substr((id + 1) / 2 - (ans / 2), ans);
    cout << t << "\n";
    return 0;
}

//DINIC
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


// MINCOST
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

// NUMTHEORY
i64 ext_gcd (i64 a, i64 b, i64 & x, i64 & y) {
    if (a == 0) { x = 0; y = 1; return b; }
    i64 x1, y1, d = ext_gcd(b%a, a, x1, y1);
    x = y1 - (b / a) * x1; y = x1;
    return d;
}

int inv(i64 x, int p) {
	int s = p - 2;
	i64 r = 1;
	while (s) {
		if (s & 1) 
			r = (r * x) % p;
		x = (x * x) % p;
		s >>= 1;
	}
	return r;
}

bool isPrime(int x) {
	if (x < 2) return 0;
	for (int i = 2; i * i <= x; i++) if (x % i == 0) return 0;
	return 1;
}

vector<int> kto_p;
vector<vi> kto_r;
void init_kto(int n = 100) {
    kto_p.clear();
    for(int x = 1e9; (int)kto_p.size() < n; ++x)
        if( isPrime(x) )
            kto_p.pb(x);
    kto_r.assign(n, vi(n));
    forn(i,n)forn(j,n)if(j>i)
        kto_r[i][j]=inv(kto_p[i], kto_p[j]);
}
int get_by_kto(vi x) {
    int r(0), m(1);
    forn(i, x.size()) {
        forn(j, i) {
            i64 cur = (i64)(x[i] - x[j]) * kto_r[j][i];
            cur %= kto_p[i]; if( cur < 0 ) cur += kto_p[i];
            x[i] = cur;                  
        }
        r = r + (m * x[i]);
        m = m * kto_p[i];
    } return r;
}

int powmod (int a, int b, int p) {
	int res = 1;
	while (b)
		if (b & 1)
			res = int (res * 1ll * a % p),  --b;
		else
			a = int (a * 1ll * a % p),  b >>= 1;
	return res;
}
 
int generator (int p) {
	vector<int> fact;
	int phi = p-1,  n = phi;
	for (int i=2; i*i<=n; ++i)
		if (n % i == 0) {
			fact.push_back (i);
			while (n % i == 0)
				n /= i;
		}
	if (n > 1)
		fact.push_back (n);
 
	for (int res=2; res<=p; ++res) {
		bool ok = true;
		for (size_t i=0; i<fact.size() && ok; ++i)
			ok &= powmod (res, phi / fact[i], p) != 1;
		if (ok)  return res;
	}
	return -1;
}

int gauss (vector < vector<double> > a, vector<double> & ans) {
	const int inf = 1e9; const double eps = 1e-8;
	int n = (int) a.size();
	int m = (int) a[0].size() - 1;
 
	vector<int> where (m, -1);
	for (int col=0, row=0; col<m && row<n; ++col) {
		int sel = row;
		for (int i=row; i<n; ++i)
			if (abs (a[i][col]) > abs (a[sel][col]))
				sel = i;
		if (abs (a[sel][col]) < eps)
			continue;
		for (int i=col; i<=m; ++i)
			swap (a[sel][i], a[row][i]);
		where[col] = row;
 
		for (int i=0; i<n; ++i)
			if (i != row) {
				double c = a[i][col] / a[row][col];
				for (int j=col; j<=m; ++j)
					a[i][j] -= a[row][j] * c;
			}
		++row;
	}
 
	ans.assign (m, 0);
	for (int i=0; i<m; ++i)
		if (where[i] != -1)
			ans[i] = a[where[i]][m] / a[where[i]][i];
	for (int i=0; i<n; ++i) {
		double sum = 0;
		for (int j=0; j<m; ++j)
			sum += ans[j] * a[i][j];
		if (abs (sum - a[i][m]) > eps)
			return 0;
	}
 
	for (int i=0; i<m; ++i)
		if (where[i] == -1)
			return inf;
	return 1;
}

//N should be multiplied by two
template <class T>
double simpson(T f, double a, double b, int N = 1000 * 1000) {
	double s = 0;
	double h = (b - a) / N;
	for (int i=0; i<=N; ++i) {
		double x = a + h * i;
		s += f(x) * ((i==0 || i==N) ? 1 : ((i&1)==0) ? 2 : 4);
	}
	s *= h / 3;
	return s;
}


template <class T, class T2>
bool miller_rabin (T n, T2 b)
{
	if (n == 2)
		return true;
	if (n < 2 || even (n))
		return false;

	if (b < 2)
		b = 2;
	for (T g; (g = gcd (n, b)) != 1; ++b)
		if (n > g)
			return false;

	T n_1 = n;
	--n_1;
	T p, q;
	transform_num (n_1, p, q);
	T rem = powmod (T(b), q, n);
	if (rem == 1 || rem == n_1)
		return true;
	for (T i=1; i<p; i++)
	{
		mulmod (rem, rem, n);
		if (rem == n_1)
			return true;
	}
	return false;
}

template <class T>
T pollard_rho (T n, unsigned iterations_count = 100000)
{
	T
		b0 = rand() % n,
		b1 = b0,
		g;
	mulmod (b1, b1, n);
	if (++b1 == n)
		b1 = 0;
	g = gcd (abs (b1 - b0), n);
	for (unsigned count=0; count<iterations_count && (g == 1 || g == n); count++)
	{
		mulmod (b0, b0, n);
		if (++b0 == n)
			b0 = 0;
		mulmod (b1, b1, n);
		++b1;
		mulmod (b1, b1, n);
		if (++b1 == n)
			b1 = 0;
		g = gcd (abs (b1 - b0), n);
	}
	return g;
}

//PLANARITY
const int N = 300000;
int n, m;

namespace Planarity {

    int edge[N][2];
    bool oriented[N];
    vector<int> G[N];
    int height[N];
    int lowpt[N];
    int lowpt2[N];
    int nesting_depth[N];
    int parent_edge[N];
    const int INF = INT_MAX >> 2;


    void DFS1(int v) {
        int e = parent_edge[v];
        for (int edgeId : G[v]) {
            if (oriented[edgeId]) continue;
            if (edge[edgeId][1] == v) {
                swap(edge[edgeId][0], edge[edgeId][1]);
            }
            oriented[edgeId] = true;
            lowpt[edgeId] = lowpt2[edgeId] = height[v];

            int to = edge[edgeId][1];
            if (height[to] == INF) {
                parent_edge[to] = edgeId;
                height[to] = height[v] + 1;
                DFS1(to);
            } else {
                lowpt[edgeId] = height[to];
            }

            /**
             * determine nesting depth
             */
            nesting_depth[edgeId] = 2 * lowpt[edgeId];
            if (lowpt2[edgeId] < height[v]) {
                /* chordal */
                nesting_depth[edgeId]++;
            }

            /**
             * update lowpoints of parent edge e
             */
            if (e != -1) {
                if (lowpt[edgeId] < lowpt[e]) {
                    lowpt2[e] = min(lowpt[e], lowpt2[edgeId]);
                    lowpt[e]  = lowpt[edgeId];
                }
                else if (lowpt[edgeId] > lowpt[e]) {
                    lowpt2[e] = min(lowpt2[e], lowpt[edgeId]);
                }
                else {
                    lowpt2[e] = min(lowpt2[e], lowpt2[edgeId]);
                }
            }
        }
    }


    int ref[N];
    int side[N];
    class I {
    public:
        int low, high;
        I(int low = -1, int high = -1) : low(low), high(high) {}
        I(const I &i) : low(i.low), high(i.high) {}
        bool empty() const {
            return low == -1;
        }
        bool operator == (const I &i) const {
            return low == i.low && high == i.high;
        }
    };
    class P {
    public:
        I L;
        I R;
        P(){}
        P(const I &l, const I &r) : L(l), R(r) {}
        bool operator == (const P &p) const {
            return L == p.L && R == p.R;
        }
    };
    int lowpt_edge[N];
    P stack_bottom[N];
    vector<P> S;


    inline bool conflicting(I i, int edgeId) {
        return (!i.empty() && lowpt[i.high] > lowpt[edgeId]);
    }

    int addConstraints(int edgeId, int e) {
        P p = P(I(), I());
        /**
         * merge return edges of edgeId into p.R
         */
        do {
            P q = (S.empty() ? P(I(), I()) : S.back());
            if (!S.empty()) S.pop_back();
            if (!q.L.empty()) swap(q.L, q.R);
            if (!q.L.empty()) {
                /* not planar */
                return 0;
            }
            if (lowpt[q.R.low] > lowpt[e]) {
                /* merge intervals */
                if (p.R.empty()) p.R.high     = q.R.high;
                else             ref[p.R.low] = q.R.high;

                p.R.low = q.R.low;
            } else {
                /* make consistent */
                ref[q.R.low] = lowpt_edge[e];
            }
        } while (!((S.empty() ? P(I(), I()) : S.back()) == stack_bottom[edgeId]));
        /**
         * merge conflicting return edges of e1, . . . , eiâ��1 into p.L
         */
        P qw = (S.empty() ? P(I(), I()) : S.back());
        while (conflicting((S.empty() ? P(I(), I()) : S.back()).L, edgeId) ||
               conflicting((S.empty() ? P(I(), I()) : S.back()).R, edgeId)) {
            P q = (S.empty() ? P(I(), I()) : S.back());
            if (!S.empty()) S.pop_back();
            if (conflicting(q.R, edgeId)) swap(q.L, q.R);
            if (conflicting(q.R, edgeId)) {
                /* not planar */
                return 0;
            }
            /* merge interval below lowpt(edgeId) into p.R */
            ref[p.R.low] = q.R.high;
            if (q.R.low >= 0) p.R.low = q.R.low;

            if (p.L.empty()) {
                p.L.high = q.L.high;
            } else {
                ref[p.L.low] = q.L.high;
            }
            p.L.low = q.L.low;
        }
        if (!p.L.empty() || !p.R.empty()) {
            S.push_back(p);
        }
        return 1;
    }

    inline int lowest(P p) {
        if (p.L.empty()) return lowpt[p.R.low];
        if (p.R.empty()) return lowpt[p.L.low];
        return std::min(lowpt[p.L.low], lowpt[p.R.low]);
    }

    void trimBackEdges(int u) {
        while (!S.empty() && lowest(S.back()) == height[u]) {
            P p = S.back();
            S.pop_back();
            if (p.L.low >= 0) side[p.L.low] = -1;
        }
        if (!S.empty()) {
            P p = S.back();
            S.pop_back();
            /*
             * trim left interval
             */
            while (p.L.high != -1 && edge[p.L.high][1] == u) {
                p.L.high = ref[p.L.high];
            }
            if (p.L.high == -1 && p.L.low != -1) {
                ref[p.L.low] = p.R.low;
                side[p.L.low] = -1;
                p.L.low = -1;
            }

            /*
             * trim right interval
             */
            while (p.R.high != -1 && edge[p.R.high][1] == u) {
                p.R.high = ref[p.R.high];
            }
            if (p.R.high == -1 && p.R.low != -1) {
                ref[p.R.low] = p.L.low;
                side[p.R.low] = -1;
                p.R.low = -1;
            }
            S.push_back(p);
        }
    }

    int DFS2(int v) {
        int e = parent_edge[v];
        int E1 = -1;
        for (int edgeId : G[v]) {
            if (edge[edgeId][0] != v) continue;
            if (E1 == -1) {E1 = edgeId;}
            stack_bottom[edgeId] = (S.empty() ? P(I(), I()) : S.back());
            int to = edge[edgeId][1];
            if (edgeId == parent_edge[to]) {
                int res = DFS2(to);
                if (res == 0) {
                    return 0;
                }
            } else {
                lowpt_edge[edgeId] = edgeId;
                S.push_back(P(I(), I(edgeId, edgeId)));
            }
            if (lowpt[edgeId] < height[v]) {
                if (edgeId == E1) {
                    lowpt_edge[e] = lowpt_edge[edgeId];
                } else {
                    int res = addConstraints(edgeId, e);
                    if (res == 0) {
                        return 0;
                    }
                }
            }
        }

        /* v is not root */
        if (e != -1) {
            int u = edge[e][0];
            trimBackEdges(u);
            /**
             * side of e is side of a highest return edge
             */
            /* e has return edge */
            if (lowpt[e] < height[u]) {
                P back = (S.empty() ? P(I(), I()) : S.back());
                int hL = back.L.high;
                int hR = back.R.high;
                if (hL != -1 && (hR == -1 || lowpt[hL] > lowpt[hR])) {
                    ref[e] = hL;
                } else {
                    ref[e] = hR;
                }
            }
        }
        return 1;
    }

    int checkPlanarity() {
        /**
        * orientation
        */
        vector<int> roots;
        for (int i = 0; i < n; ++i) {
            height[i] = INF;
        }
        for (int i = 0; i < m; ++i) {
            oriented[i] = false;
        }
        for (int i = 0; i < n; ++i) {
            if (height[i] == INF) {
                height[i] = 0;
                parent_edge[i] = -1;
                roots.push_back(i);
                DFS1(i);
            }
        }

        /**
        * testing
        */
        // sort adjacency lists according to non-decreasing nesting depth
        for (int i = 0; i < n; ++i) {
            stable_sort(G[i].begin(), G[i].end(), [&](int edgeLeft, int edgeRight) {
                return nesting_depth[edgeLeft] < nesting_depth[edgeRight];
            });
        }
        for (int i = 0; i < m; ++i) {
            ref[i] = -1;
            side[i] = 1;
        }
        S.clear();
        for (int root : roots) {
            int res = DFS2(root);
            if (res == 0) {
                return 0;
            }
        }
        return 1;
    }
};

int solve() {
    n = nxt();
    m = nxt();
    for (int i = 0; i < n; ++i) {
        Planarity::G[i].clear();
    }
    for (int i = 0; i < m; ++i) {
        int a = nxt() - 1;
        int b = nxt() - 1;
        Planarity::edge[i][0] = a;
        Planarity::edge[i][1] = b;
        Planarity::G[a].push_back(i);
        Planarity::G[b].push_back(i);
    }

    if (n == 1) {
        return 1;
    }
    if (n == 2) {
        return m == 1;
    }

    if (m != 3 * (n - 2)) {
        return 0;
    }
    return Planarity::checkPlanarity();
}

// HALFPLANES INTERSECTION
const double eps = 1e-9;
using namespace std;

typedef ptd pt;

int sgn(const double &a) {
    if (a > eps)
        return 1;
    else if (a < -eps)
        return -1;
    else
        return 0;
}

struct halfPlane {
    pt a, dir;
    halfPlane() {}
    halfPlane(const pt &_a, const pt &_dir) : a(_a), dir(_dir) {}

    int type() const {
        if (sgn(dir.y) != 0)
            return sgn(dir.y);
        return sgn(dir.x);
    }

    bool operator<(const halfPlane &rhs) const {
        if (type() != rhs.type())
            return type() < rhs.type();
        double cr = cross(dir, rhs.dir);
        if (sgn(cr) != 0) {
            return sgn(cr) < 0;
        }
        return !include(rhs.a);
    }

    bool operator==(const halfPlane &rhs) const {
        if (type() != rhs.type())
            return false;
        double cr = cross(dir, rhs.dir);
        return sgn(cr) == 0;
    }

    bool include(const pt &p) const {
        return sgn(cross(dir, p - a)) >= 0;
    }
};

pt getIntersection(const halfPlane &lhs, const halfPlane &rhs) {
    double u = cross(lhs.a - rhs.a, rhs.dir);
    double v = cross(lhs.dir, rhs.dir);
    double t = -u / v;
    return lhs.a + lhs.dir * t;
}

const int N = 111111;

halfPlane hp[N];

int n;

int bot, top;

inline void intersect() {
    top = bot = 0;

    for (int i = 0; i < n; i++) {
        while (bot - top >= 2 && !hp[i].include(getIntersection(hp[bot - 1], hp[bot - 2])))
            bot--;
        while (bot - top >= 2 && !hp[i].include(getIntersection(hp[top], hp[top + 1])))
            top++;
        hp[bot++] = hp[i];
    }
    while (bot - top >= 2 && !hp[top].include(getIntersection(hp[bot - 1], hp[bot - 2])))
        bot--;
    while (bot - top >= 2 && !hp[bot - 1].include(getIntersection(hp[top], hp[top + 1])))
        top++;
}

int main() {
    while (scanf("%d", &n) == 1) {
        int sz = 0;
        for (int i = 0; i < n; ++i) {
            //double a = nxt(), b = nxt(), c = nxt(), d = nxt();
            double a, b, c, d;
            scanf("%lf%lf%lf%lf", &a, &b, &c, &d);
            pt A(a, b);
            pt dir(c - a, d - b);
            dir = dir * (1.0 / dir.len());
            hp[sz].a = A;
            hp[sz++].dir = dir;
        }
        n = sz;

        hp[n++] = halfPlane(pt(0, 0), pt(10000, 0));
        hp[n++] = halfPlane(pt(10000, 0), pt(0, 10000));
        hp[n++] = halfPlane(pt(10000, 10000), pt(-10000, 0));
        hp[n++] = halfPlane(pt(0, 10000), pt(0, -10000));

        sort(hp, hp + n);
        n = (int) (unique(hp, hp + n) - hp);

        intersect();

        if (bot - top < 2) {
            cout << "0.0" << endl;
            return 0;
        }

        double sq = 0;

        pt pts[bot - top];
        sz = 0;

        for (int i = top; i < bot; ++i) {
            int k = i + 1;
            if (k == bot) {
                k = top;
            }
            pts[sz++] = getIntersection(hp[i], hp[k]);
        }

        for (int i = 0; i < sz; ++i) {
            sq += cross(pts[i], pts[(i + 1) % sz]);
        }
        sq = abs(sq);
        sq /= 2;
        cout << setprecision(1) << fixed;
        printf("%.1f\n", sq);
    }
    return 0;
}

// SUFFIX ARRAY

struct SuffixArray {
    int *str;
    int *sa;
    int *ra;
    int *lcp;
    int *lgt;
    int **rmq;
    int n;
    SuffixArray() : str(0), sa(0), ra(0), lcp(0), lgt(0), rmq(0), n(0) {}

    void init(const vector<int> &s) {
        n = (int)s.size();
        str = new int[n];
        for (int i = 0; i < n; ++i) str[i] = s[i];
        __init();
    }
    void init(const string &s) {
        n = (int)s.length();
        str = new int[n];
        for (int i = 0; i < n; ++i) str[i] = s[i];
        __init();
    }
    void __init() {
        initSa();
        initLcp();
    }
    void initSa() {
        sa = new int[n];
        ra = new int[n];
        int *s = new int[n];
        int *c = new int[n];
        int *cl = new int[n];
        int *cnt = new int[n];
        for (int i = 0; i < n; ++i) {
            cl[i] = str[i];
            sa[i] = n - i - 1;
        }
        stable_sort(sa, sa + n, [&](int l, int r){ return str[l] < str[r]; });
        for (int len = 1; len < n; len <<= 1) {
            memcpy(c, cl, n * sizeof(int));
            memcpy(s, sa, n * sizeof(int));
            for (int i = 0; i < n; ++i) cnt[i] = i;
            int len2 = len >> 1;
            for (int i = 0; i < n; ++i) {
                if (i > 0 && c[sa[i]] == c[sa[i - 1]] && sa[i - 1] + len < n && c[sa[i] + len2] == c[sa[i - 1] + len2]) {
                    cl[sa[i]] = cl[sa[i - 1]];
                } else {
                    cl[sa[i]] = i;
                }
            }
            for (int i = 0; i < n; ++i) {
                int s1 = s[i] - len;
                if (s1 >= 0) {
                    sa[cnt[cl[s1]]++] = s1;
                }
            }
        }
        for (int i = 0; i < n; ++i) {
            ra[sa[i]] = i;
        }
        delete [] s;delete [] c;delete [] cl;delete [] cnt;
    }
    void initLcp() {
        lcp = new int[n];
        for (int i = 0, h = 0; i < n; ++i) {
            if (ra[i] < n - 1) {
                for (int j = sa[ra[i] + 1]; max(i, j) + h < n && str[i + h] == str[j + h]; ++h);
                lcp[ra[i]] = h;
                h = max(0, h - 1);
            } else {
                lcp[ra[i]] = 0;
            }
        }
        lgt = new int[n + 1];
        lgt[0] = lgt[1] = 0;
        for (int i = 2; i <= n; ++i) lgt[i] = lgt[i >> 1] + 1;
        rmq = new int*[lgt[n] + 1];
        for (int i = 0; i <= lgt[n]; ++i) {
            rmq[i] = new int[n];
        }
        for (int i = 0; i + 1 < n; ++i) {
            rmq[0][i] = lcp[i];
        }
        for (int l = 1; (1 << l) < n; ++l) {
            for (int i = 0; i + (1 << l) <= n; ++i) {
                rmq[l][i] = min(rmq[l - 1][i], rmq[l - 1][i + (1 << (l - 1))]);
            }
        }
    }
    int getLcpForSuffixes(int l, int r) {
        if (l == r) return n - l;
        l = ra[l]; r = ra[r];
        if (l > r) swap(l, r);
        --r;
        int j = lgt[r - l];
        return min(rmq[j][l], rmq[j][r - (1 << j) + 1]);
    }
    int getLcpForSuffixesIndexes(int l, int r) {
        if (l == r) return n - sa[l];
        if (l > r) swap(l, r);
        --r;
        int j = lgt[r - l];
        return min(rmq[j][l], rmq[j][r - (1 << j) + 1]);
    }
    pair<int, int> getIntervalForSuffixIndex(int sufIndex, int length) {
        int l, r;
        l = sufIndex;
        r = n;
        while (l < r) {
            int m = (l + r) >> 1;
            if (getLcpForSuffixesIndexes(sufIndex, m) >= length) {
                l = m + 1;
            } else {
                r = m;
            }
        }
        l = 0, r = sufIndex;
        while (l < r) {
            int m = (l + r) >> 1;
            if (getLcpForSuffixesIndexes(sufIndex, m) >= length) {
                r = m;
            } else {
                l = m + 1;
            }
        }
        return make_pair(l, r);
    }
    ~SuffixArray() {
        if (str) delete [] str;
        if (sa) delete [] sa;
        if (ra) delete [] ra;
        if (lgt) {
            for (int i = 0; i <= lgt[n]; ++i) {
                delete [] rmq[i];
            }
            delete [] lgt;
            delete [] rmq;
        }
    }
};

// SUFFIX AUTOMATA
const int MAXLEN = 2 * 1000 * 1000;
const int ALF = 26;

struct State {
    int len, link;
    int firstPos;
    int next[ALF];
    long long cntSubs;
    long long lenSubs;
    int occurences;
};

State st[MAXLEN * 2];
int sz, last;

void init(int v) {
    st[0].occurences = 0;
    st[0].cntSubs = -1;
    st[0].lenSubs = -1;
    memset(st[v].next, -1, sizeof(st[v].next));
}

void sa_init() {
    st[0].len = 0;
    st[0].link = -1;
    init(0);
    last = 0; sz = 1;
}

void sa_extend(int c) {
    c -= 'a';
    int cur = sz++;
    st[cur].len = st[last].len + 1;
    init(cur);
    st[cur].occurences = 1;
    st[cur].firstPos = st[cur].len;
    int p = last;
    for (;p != -1 && st[p].next[c] == -1; p = st[p].link) {
        st[p].next[c] = cur;
    }
    if (p == -1) {
        st[cur].link = 0;
    } else {
        int q = st[p].next[c];
        if (st[p].len + 1 == st[q].len) {
            st[cur].link = q;
        } else {
            int cl = sz++;
            init(cl);
            memcpy(st[cl].next, st[q].next, sizeof(st[q].next));
            st[cl].firstPos = st[q].firstPos;
            st[cl].link = st[q].link;
            st[cl].len = st[p].len + 1;
            for (; p != -1 && st[p].next[c] == q; p = st[p].link) {
                st[p].next[c] = cl;
            }
            st[q].link = st[cur].link = cl;
        }
    }
    last = cur;
}

void calcDifferentSubstrings(int v) {
    if (st[v].cntSubs != -1) return;
    st[v].cntSubs = 1;
    st[v].lenSubs = 0;
    for (int i = 0; i < ALF; ++i) {
        int to = st[v].next[i];
        if (to == -1) continue;
        calcDifferentSubstrings(to);
        st[v].cntSubs += st[to].cntSubs;
        st[v].lenSubs = st[to].cntSubs + st[to].lenSubs;
    }
}

void calcOccurences() {
    vector<vector<int> > states(sz);
    for (int i = 0; i < sz; ++i) {
        states[st[i].len].push_back(i);
    }
    for (int i = sz - 1; i > 0; --i) {
        for (int id: states[i]) {
            st[st[id].link].occurences += st[id].occurences;
        }
    }
}

int getOccurences(const string &s) {
    int v = 0;
    for (char c : s) {
        if (st[v].next[c - 'a'] == -1) return 0;
        v = st[v].next[c];
    }
    return st[v].occurences;
}


// TANGENTS
void tangents (pt c, double r1, double r2, vector<line> & ans) {
    double r = r2 - r1;
    double z = sqr(c.x) + sqr(c.y);
    double d = z - sqr(r);
    if (d < -EPS)  return;
    d = sqrt (abs (d));
    line l;
    l.a = (c.x * r + c.y * d) / z;
    l.b = (c.y * r - c.x * d) / z;
    l.c = r1;
    ans.push_back (l);
}

vector<line> tangents (circle a, circle b) {
    vector<line> ans;
    for (int i=-1; i<=1; i+=2)
        for (int j=-1; j<=1; j+=2)
            tangents (b-a, a.r*i, b.r*j, ans);
    for (size_t i=0; i<ans.size(); ++i)
        ans[i].c -= ans[i].a * a.x + ans[i].b * a.y;
    return ans;
}

//WALSH
template <class T>
void walshTransform(T * data, int n) {
    for (int len = 2; len <= n; len <<= 1) {
        int len2 = len >> 1;
        for (int r = 0; r < n; r += len) {
            int p1 = r, p2 = p1 + len2;
            for (int j = 0; j < len2; ++j, ++p1, ++p2) {
                T u = data[p1];
                T v = data[p2];
                data[p1] = u + v;
                data[p2] = u - v;
            }
        }
    }
}
template<class T>
void xorConvolution(T *left, T *right, int n) {
    walshTransform(left, n);
    walshTransform(right, n);
    for (int i = 0; i < n; ++i) {
        left[i] = left[i] * right[i];
    }
    walshTransform(left, n);
    for (int i = 0; i < n; ++i) {
        left[i] /= n;
    }
}

template<class T>
void xorConvolutionSquare(T *left, int n) {
    walshTransform(left, n);
    for (int i = 0; i < n; ++i) {
        left[i] = left[i] * left[i];
    }
    walshTransform(left, n);
    for (int i = 0; i < n; ++i) {
        left[i] /= n;
    }
}
