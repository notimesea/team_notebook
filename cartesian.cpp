#include <bits/stdc++.h>

#define fi first
#define se second
using namespace std;

typedef vector<int> vi;
typedef long long i64;


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

void testCartesian() {
	pitem root = NULL;
	vi a;
	vector<pitem> els;
    for (int i = 0; i < 5000; i++) {
		int type = rand() % 3;
		if (type == 0) {
			int x = rand();
			pitem l=0, r=0;
			split2(root, x, l, r);
			els.push_back(new item(x, rand()));
			l = merge(l, els.back());
			root = merge(l, r);
			a.push_back(x);
		}
		if (type == 1) {
			int x = rand();
			pitem l=0, r=0;
			split2(root, x, l, r);
			i64 sl = 0, sr = 0;
			int ql = 0, qr = 0;
			for (int y: a)
				if (y <= x) sl += y, ql++;
				else sr += y, qr++;
			assert(sl == sum(l) && sr == sum(r));
			assert(ql == q(l) && qr == q(r));
			for (auto p: els) {
				auto y = p;
				while (y) {
					assert((p->key <= x) == (y->key <= x));
					y = y->par;
				}
			}
			root = merge(l, r);
		}
	}	
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

// tested here hackerearth.com/april-circuits/algorithm/circ-bear-and-leaderboard-1/
int main() {
    testCartesian();
	return 0;
}
