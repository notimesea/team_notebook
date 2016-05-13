//
// Created by Sergey Kiyan on 01.05.16.
//

#include <bits/stdc++.h>
using namespace std;

const int MAXLEN = 2 * 1000 * 1000;

struct State {
    int len, link;
    int firstPos;
    map<char, int> next;
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
}

void sa_init() {
    for (int i = 0; i < sz; ++i) {
        st[i].next.clear();
    }
    st[0].len = 0;
    st[0].link = -1;
    init(0);
    last = 0; sz = 1;
}

void sa_extend(char c) {
    int cur = sz++;
    st[cur].len = st[last].len + 1;
    init(cur);
    st[cur].occurences = 1;
    st[cur].firstPos = st[cur].len;
    int p = last;
    for (;p != -1 && !st[p].next.count(c); p = st[p].link) {
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
            st[cl].next = st[q].next;
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
    for (pair<char, int> pa : st[v].next) {
        int to = pa.second;
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
        if (!st[v].next.count(c)) return 0;
        v = st[v].next[c];
    }
    return st[v].occurences;
}




// TEST
set<string> substringsTest(const string &s) {
    set<string> ret;
    for (size_t i = 0; i < s.length(); ++i) {
        for (size_t len = 0; i + len <= s.length(); ++i) {
            ret.insert(s.substr(i, len));
        }
    }
    return ret;
}

int getOccurencesTest(const string &s, const string &t) {
    int ret = 0;
    for (size_t i = 0; i + s.length() <= t.length(); ++i) {
        if (t.substr(i, s.length()) == s) {
            ret++;
        }
    }
    return ret;
}

void test() {
    const int N = 100;
    const int K = 2;
    const int IT = 100;
    int it = IT;
    while (it--) {
        string str(N, 'a');
        for (size_t i = 0; i < str.length(); ++i) {
            str[i] = char('a' + rand() % K);
        }
        sa_init();
        for (char c : str) {
            sa_extend(c);
        }
        calcDifferentSubstrings(0);
        calcOccurences();

        set<string> subs = substringsTest(str);
        assert(st[0].cntSubs == (int)subs.size());
        long long totalLength = 0;
        for (const string &s : subs) {
            totalLength += s.length();
        }
        assert(st[0].lenSubs == totalLength);
        for (const string &s : subs) {
            if (s.empty()) continue;
            assert(getOccurencesTest(s, str) == getOccurences(s));
        }
    }
}


void testPerformance() {
    vector<int> sizesForTest({10000, 100000, 200000, 400000, 500000, 1000000, 2000000});
    const int K = 2;
    const int IT = 4;

    for (int size : sizesForTest) {
        int it = IT;
        double avgTime = 0;
        while (it--) {
            string str((size_t) size, 'a');
            for (int i = 0; i < size; ++i) {
                str[i] = char('a' + (rand() % K));
            }
            double cl0 = clock();
            sa_init();
            for (char c : str) {
                sa_extend(c);
            }
            calcDifferentSubstrings(0);
            calcOccurences();
            double cl1 = clock();
            avgTime += (cl1 - cl0) / CLOCKS_PER_SEC;
        }
        avgTime /= IT;
        cerr << "Time for size " << size << " is " << avgTime << "s." << endl;
    }
}

int main() {
    test();
    testPerformance();
}