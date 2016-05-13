//
// Created by Sergey Kiyan on 01.05.16.
//

#include <bits/stdc++.h>
using namespace std;
typedef long long ll;


typedef long double D;
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



// TEST

vector<ll> mulTest(const vector<ll> &a, const vector<ll> &b) {
    if (a.empty()) return vector<ll>();
    if (b.empty()) return vector<ll>();
    vector<ll> c(a.size() + b.size() - 1, 0);
    for (size_t i = 0; i < a.size(); ++i) {
        for (size_t j = 0; j < b.size(); ++j) {
            c[i + j] += a[i] * b[j];
        }
    }
    while (!c.empty() && c.back() == 0) {
        c.pop_back();
    }
    return c;
}

void testFFT() {
    const int L = 100;
    const int V = 100;
    const int IT = 100;
    initFFT();
    int it = IT;
    while (it--) {
        vector<ll> a(rand() % L);
        vector<ll> b(rand() % L);
        for (size_t i = 0; i < a.size(); ++i) {
            a[i] = rand() % V;
        }
        for (size_t i = 0; i < b.size(); ++i) {
            b[i] = rand() % V;
        }
        assert(mulFFT(a, b) == mulTest(a, b));
    }
}

void testPerformance() {
    {
        double cl0 = clock();
        initFFT();
        double cl1 = clock();
        double time = (cl1 - cl0) / CLOCKS_PER_SEC;
        cerr << "Time for init is " << time << "s." << endl;
    }

    vector<int> sizesForTest({10000, 100000, 200000, 400000, 500000, 1000000});
    const int IT = 4;
    for (int size : sizesForTest) {
        int it = IT;
        double avgTime = 0;
        while (it--) {
            vector<ll> a(size);
            vector<ll> b(size);
            for (int i = 0; i < size; ++i) {
                a[i] = rand() % 10000;
            }
            for (int i = 0; i < size; ++i) {
                b[i] = rand() % 10000;
            }
            double cl0 = clock();
            a = mulFFT(a, b);
            double cl1 = clock();
            avgTime += (cl1 - cl0) / CLOCKS_PER_SEC;
        }
        avgTime /= IT;
        cerr << "Time for size " << size << " is " << avgTime << "s." << endl;
    }
}


int main() {
    testFFT();
    testPerformance();
}
