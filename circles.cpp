//
// Created by Сергей Киян on 13.05.16.
//

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
}