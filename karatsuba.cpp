#include <bits/stdc++.h>

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


typedef long long T;
karatsuba<T> k;

void testSpeed() {
    double start = clock();
    const int N = 1 << 18;
    T* a = k.allocate(N);
    T* b = k.allocate(N);
    for (int i = 0; i < N; ++i) {
        a[i] = rand();
        b[i] = rand();
    }
    T* c = k.multiply(a, b, N);
    k.clear(N + N + 2*N);
    assert(k.pos == k.buffer);
    double end = clock();
    printf("Elapsed time: %f\n", 1.0 * (end - start) / CLOCKS_PER_SEC); //1.4 on mac
}


void testQuality() {
    double start = clock();
    const int N = 1 << 11;
    T* a = k.allocate(N);
    T* b = k.allocate(N);
    for (int i = 0; i < N; ++i) {
        a[i] = rand();
        b[i] = rand();
    }
    T* c = k.multiply(a, b, N);
    T* d = k.allocate(2*N);
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            d[i + j] += a[i] * b[j];
        }
    }
    for (int i = 0; i < 2*N; ++i) {
        assert(c[i] == d[i]);
    }
    k.clear(N + N + 2*N + 2*N);
    assert(k.pos == k.buffer);
    double end = clock();
    printf("Elapsed time: %f\n", 1.0 * (end - start) / CLOCKS_PER_SEC);
}

int main() {
    testQuality();
    testSpeed();
    return 0;
}