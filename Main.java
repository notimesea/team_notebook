import java.awt.geom.*;
import java.io.*;
import java.util.*;
public class Main {
    public static void main(String[] args) throws Exception {
        new Solution().solve(1, new InputReader(System.in), new PrintWriter(System.out));
    }
}

class Solution {
    int N = 1000;

    Area getEllipse(double xc, double yc, double rx, double ry) {
        Path2D res = new Path2D.Double();
        Rectangle2D rect = new Rectangle2D.Double(xc - rx, yc - ry, 2 * rx, 2 * ry);
        double init = 0;
        for (int i = 0; i < N; ++i) {
            Arc2D arc = new Arc2D.Double(rect, 360.0 / N * i + init, 360.0 / N, Arc2D.OPEN);
            res.append(arc, true);
        }
        res.closePath();
        return new Area(res);
    }

    double square(Area area) {
        double arr[] = new double[6];
        double sx = 0, px = 0, sy = 0, py = 0;
        double ans = 0;
        double res = 0;
        for (PathIterator it = area.getPathIterator(null); !it.isDone(); it.next()) {
            if (it.currentSegment(arr) == PathIterator.SEG_MOVETO) {
                sx = px = arr[0];
                sy = py = arr[1];
                continue;
            }
            if (it.currentSegment(arr) == PathIterator.SEG_CLOSE) {
                res += sx * py - sy * px;
                ans += res;
                res = 0;
                continue;
            }
            if (it.currentSegment(arr) == PathIterator.SEG_CUBICTO) {
                double x0 = px, y0 = py;
                double x1 = arr[0], y1 = arr[1];
                double x2 = arr[2], y2 = arr[3];
                double x3 = arr[4], y3 = arr[5];
                res -= 0.6 * x0 * y1 + 0.3 * x0 * y2 + 0.1 * x0 * y3 +
                        x1 * (-0.6 * y0 + 0.3 * y2 + 0.3 * y3) +
                        x2 * (-0.3 * y0 - 0.3 * y1 + 0.6 * y3) -
                        0.1 * x3 * y0 - 0.3 * x3 * y1 - 0.6 * x3 * y2;
                px = x3;
                py = y3;
                continue;
            }
            if (it.currentSegment(arr) == PathIterator.SEG_LINETO) {
                double x = arr[0];
                double y = arr[1];
                res += x * py - y * px;
                px = x;
                py = y;
                continue;
            }
            throw new Error();
        }
        return ans / 2;
    }


    public void solve(int testNumber, InputReader in, PrintWriter out) {
        int m = in.nextInt();
        Area area = new Area();

        for (int i = 0; i < m; ++i) {
            int x = in.nextInt();
            int y = in.nextInt();
            int r = in.nextInt();
            area.add(getEllipse(x, y, r, r));
        }

        int n = in.nextInt();
        Path2D poly = new Path2D.Double();
        for (int i = 0; i < n; ++i) {
            int x = in.nextInt();
            int y = in.nextInt();
            if (i == 0) {
                poly.moveTo(x, y);
            } else {
                poly.lineTo(x, y);
            }
        }
        poly.closePath();
        area.intersect(new Area(poly));
        out.printf(Locale.US, "%.10f\n", square(area));
        out.close();
    }
}

class InputReader {
    public BufferedReader reader;
    public StringTokenizer tokenizer;

    public InputReader(InputStream stream) {
        reader = new BufferedReader(new InputStreamReader(stream), 32768);
        tokenizer = null;
    }

    public String next() {
        while (tokenizer == null || !tokenizer.hasMoreTokens()) {
            try {
                tokenizer = new StringTokenizer(reader.readLine());
            } catch (IOException e) {
                throw new RuntimeException(e);
            }
        }
        return tokenizer.nextToken();
    }

    public int nextInt() {
        return Integer.parseInt(next());
    }

    public double nextDouble() {
        return Double.parseDouble(next());
    }
}