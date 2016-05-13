#include <set>  
#include <map>  
#include <queue>  
#include <math.h>  
#include <vector>  
#include <string>  
#include <stdio.h>  
#include <string.h>  
#include <stdlib.h>  
#include <iostream>  
#include <algorithm>  

#define eps 1e-10
#define pi acos(-1.0)
#define inf 107374182
#define inf64 1152921504606846976
#define lc l,m,tr<<1
#define rc m + 1,r,tr<<1|1
#define zero(a) fabs(a)<eps
#define iabs(x)  ((x) > 0 ? (x) : -(x))
#define clear1(A, X, SIZE) memset(A, X, sizeof(A[0]) * (SIZE))
#define clearall(A, X) memset(A, X, sizeof(A))
#define memcopy1(A , X, SIZE) memcpy(A , X ,sizeof(X[0])*(SIZE))
#define memcopyall(A, X) memcpy(A , X ,sizeof(X))
#define max( x, y )  ( ((x) > (y)) ? (x) : (y) )
#define min( x, y )  ( ((x) < (y)) ? (x) : (y) )

using namespace std;

const int maxn = 30005;
int dq[maxn], top, bot, pn, order[maxn], ln,num;
struct Point
{
    double x, y;
} p[maxn];
struct Line
{
    Point a, b;
    double angle;
} l[maxn];
int dblcmp(double k)
{
    if(fabs(k) < eps) return 0;
    return k > 0 ? 1 : -1;
}
double multi(Point p0, Point p1, Point p2)
{
    return (p1.x-p0.x)*(p2.y-p0.y) - (p1.y-p0.y)*(p2.x-p0.x);
}
bool cmp(int u, int v)
{
    int d = dblcmp(l[u].angle-l[v].angle);
    if(!d) return dblcmp(multi(l[u].a, l[v].a, l[v].b)) > 0;
    //大于0取向量左半部分为半平面，小于0，取右半部分  
    return d < 0;
}
void getIntersect(Line l1, Line l2, Point& p)
{
    double dot1, dot2;
    dot1 = multi(l2.a, l1.b, l1.a);
    dot2 = multi(l1.b, l2.b, l1.a);
    p.x = (l2.a.x * dot2 + l2.b.x * dot1) / (dot2 + dot1);
    p.y = (l2.a.y * dot2 + l2.b.y * dot1) / (dot2 + dot1);
}
bool judge(Line l0, Line l1, Line l2)
{
    Point p;
    getIntersect(l1, l2, p);
    return dblcmp(multi(p, l0.a, l0.b)) < 0;
    //大于小于符号与上面cmp（）中注释处相反  
}
void addLine(double x1, double y1, double x2, double y2)
{
    l[ln].a.x = x1;
    l[ln].a.y = y1;
    l[ln].b.x = x2;
    l[ln].b.y = y2;
    l[ln].angle = atan2(y2-y1, x2-x1);
    order[ln] = ln;
    ln++;
}
void halfPlaneIntersection()
{
    int i, j;
    sort(order, order+ln, cmp);
    for(i = 1, j = 0; i < ln; i++)
        if(dblcmp(l[order[i]].angle-l[order[j]].angle) > 0)
            order[++j] = order[i];
    ln = j + 1;
    dq[0] = order[0];
    dq[1] = order[1];
    bot = 0;
    top = 1;
    for(i = 2; i < ln; i++)
    {
        while(bot < top && judge(l[order[i]], l[dq[top-1]], l[dq[top]]))
            top--;
        while(bot < top && judge(l[order[i]], l[dq[bot+1]], l[dq[bot]]))
            bot++;
        dq[++top] = order[i];
    }
    while(bot < top && judge(l[dq[bot]], l[dq[top-1]], l[dq[top]])) top--;
    while(bot < top && judge(l[dq[top]], l[dq[bot+1]], l[dq[bot]])) bot++;
    num=0;
    for(i=bot; i<top; i++)
    {
        getIntersect(l[dq[i]],l[dq[i+1]],p[num++]);
    }
    if(top > bot+1)getIntersect(l[dq[top]], l[dq[bot]],p[num++]);

}
double xmul(Point p0,Point p1,Point p2)
{
    return (p1.x-p0.x)*(p2.y-p0.y)-(p2.x-p0.x)*(p1.y-p0.y);
}

double Get_area()
{
    double area = 0;
    for(int i = 1; i < num-1; i++)
        area += xmul(p[0], p[i], p[i+1]);
    return fabs(area)/2.0;
}
int main()
{
    int i;
    double x1,x2,y1,y2;
    scanf ("%d", &pn);
    for(ln = i = 0; i < pn; i++)
    {
        scanf("%lf%lf%lf%lf",&x1,&y1,&x2,&y2);
        addLine(x1, y1, x2,y2);
    }
    addLine(0,0,10000,0);
    addLine(10000,0,10000,10000);
    addLine(10000,10000,0,10000);
    addLine(0,10000,0,0);
    halfPlaneIntersection();
    printf("%.1f\n",Get_area());
    return 0;
}
