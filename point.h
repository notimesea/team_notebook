#ifndef TEAM_NOTEBOOK_INNOPOLIS_U_POINT_H
#define TEAM_NOTEBOOK_INNOPOLIS_U_POINT_H


template <class T>
struct pt {
    T x, y;
    pt() : x(), y() {}
    pt(T x, T y) : x(x), y(y) {}
    pt operator + (const pt &r) const {
        return pt(x + r.x, y + r.y);
    }
    pt operator - (const pt &r) const {
        return pt(x - r.x, y - r.y);
    }
    pt operator * (const T &r) const {
        return pt(x * r, y * r);
    }
    pt rot(T co, T si) const {
        return pt(x * co - y * si, x * si + y * co);
    }
    pt rot(T ang) const {
        return rot(cos(ang), sin(ang));
    }
};

template <class T>
T dot(const pt<T> &l, const pt<T> &r) {
    return l.x * r.x + l.y * r.y;
}

template <class T>
T cross(const pt<T> &l, const pt<T> &r) {
    return l.x * r.y - l.y * r.x;
}

typedef pt<int> pti;
typedef pt<long long> ptl;
typedef pt<double> ptd;
typedef pt<long double> ptld;


#endif //TEAM_NOTEBOOK_INNOPOLIS_U_POINT_H
