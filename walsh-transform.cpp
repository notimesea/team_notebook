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
