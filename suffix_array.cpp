//
// Created by Sergey Kiyan on 01.05.16.
//

#include <bits/stdc++.h>
using namespace std;

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


int getLcpTest(const string &l, const string &r) {
    size_t pos = 0;
    while (pos < l.length() && pos < r.length() && l[pos] == r[pos]) {
        ++pos;
    }
    return int(pos);
}

void testSa() {
    const int IT = 100;
    const int N = 100;
    const int K = 2;
    int it = IT;
    while (it--) {
        string str(N, 'a');
        for (int i = 0; i < N; ++i) {
            str[i] = char('a' + (rand() % K));
        }
        SuffixArray suffixArray;
        suffixArray.init(str);
        vector<int> sa(str.length());
        iota(sa.begin(), sa.end(), 0);
        sort(sa.begin(), sa.end(), [&](int l, int r) {
            return str.substr((size_t)l) < str.substr((size_t)r);
        });
        for (size_t i = 0; i < str.length(); ++i) {
            assert(suffixArray.sa[i] == sa[i]);
        }
        vector<int> ra(str.length());
        for (size_t i = 0; i < str.length(); ++i) {
            ra[sa[i]] = (int) i;
        }
        for (size_t i = 0; i < str.length(); ++i) {
            suffixArray.ra[i] = ra[i];
        }
        for (size_t i = 0; i < str.length(); ++i) {
            for (size_t j = 0; j < str.length(); ++j) {
                int lcp = getLcpTest(str.substr(i), str.substr(j));
                assert(lcp == suffixArray.getLcpForSuffixes((int)i, (int)j));
            }
        }
        for (size_t i = 0; i < str.length(); ++i) {
            for (size_t j = 0; j < str.length(); ++j) {
                int lcp = getLcpTest(str.substr((size_t)sa[i]), str.substr((size_t)sa[j]));
                assert(lcp == suffixArray.getLcpForSuffixesIndexes((int)i, (int)j));
            }
        }
        // maybe
        // TODO test for getIntervalForSuffixIndexs
    }
}

void testPerformance() {
    vector<int> sizesForTest({10000, 100000, 200000, 400000, 500000, 1000000});
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
            SuffixArray suffixArray;
            suffixArray.init(str);
            double cl1 = clock();
            avgTime += (cl1 - cl0) / CLOCKS_PER_SEC;
        }
        avgTime /= IT;
        cerr << "Time for size " << size << " is " << avgTime << "s." << endl;
    }
}



int main() {
    testSa();
    testPerformance();
}
