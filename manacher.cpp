#include <bits/stdc++.h>
using namespace std;

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