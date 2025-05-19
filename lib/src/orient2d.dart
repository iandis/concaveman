import 'dart:math';
import 'dart:typed_data';

// Constants
const double _epsilon = 1.1102230246251565e-16;
const double _resulterrbound = (3 + 8 * _epsilon) * _epsilon;

const double _ccwerrboundA = (3 + 16 * _epsilon) * _epsilon;
const double _ccwerrboundB = (2 + 12 * _epsilon) * _epsilon;
const double _ccwerrboundC = (9 + 64 * _epsilon) * _epsilon * _epsilon;

Float64List _vec(int n) => Float64List(n);

// Estimate the sum of an expansion
double _estimate(int elen, Float64List e) {
  double Q = e[0];
  for (int i = 1; i < elen; i++) Q += e[i];
  return Q;
}

// Fast expansion sum with zero elimination
int _sum(int elen, Float64List e, int flen, Float64List f, Float64List h) {
  double Q, Qnew, hh;
  int eindex = 0, findex = 0, hindex = 0;
  double enow = e[eindex], fnow = f[findex];

  if ((fnow > enow) == (fnow > -enow)) {
    Q = enow;
    enow = (++eindex < elen) ? e[eindex] : 0;
  } else {
    Q = fnow;
    fnow = (++findex < flen) ? f[findex] : 0;
  }

  while (eindex < elen && findex < flen) {
    if ((fnow > enow) == (fnow > -enow)) {
      Qnew = Q + enow;
      hh = enow - (Qnew - Q);
      enow = (++eindex < elen) ? e[eindex] : 0;
    } else {
      Qnew = Q + fnow;
      hh = fnow - (Qnew - Q);
      fnow = (++findex < flen) ? f[findex] : 0;
    }
    Q = Qnew;
    if (hh != 0) h[hindex++] = hh;
  }

  while (eindex < elen) {
    Qnew = Q + enow;
    hh = enow - (Qnew - Q);
    Q = Qnew;
    if (hh != 0) h[hindex++] = hh;
    enow = (++eindex < elen) ? e[eindex] : 0;
  }

  while (findex < flen) {
    Qnew = Q + fnow;
    hh = fnow - (Qnew - Q);
    Q = Qnew;
    if (hh != 0) h[hindex++] = hh;
    fnow = (++findex < flen) ? f[findex] : 0;
  }

  if (Q != 0 || hindex == 0) h[hindex++] = Q;
  return hindex;
}

// Core orient2d function
double orient2d(Point<double> a, Point<double> b, Point<double> c, ) {
  final detleft = (a.y - c.y) * (b.x - c.x);
  final detright = (a.x - c.x) * (b.y - c.y);
  final det = detleft - detright;

  final detsum = (detleft.abs() + detright.abs());
  if (det.abs() >= _ccwerrboundA * detsum) return det;

  return -_orient2dAdapt(a.x, a.y, b.x, b.y, c.x, c.y, detsum);
}

// Full adaptive precision version
double _orient2dAdapt(
  double ax,
  double ay,
  double bx,
  double by,
  double cx,
  double cy,
  double detsum,
) {
  final acx = ax - cx;
  final bcx = bx - cx;
  final acy = ay - cy;
  final bcy = by - cy;

  final B = _vec(4);
  _crossProduct(acx, bcx, acy, bcy, B);
  final det = _estimate(4, B);

  if (det.abs() >= _ccwerrboundB * detsum) return det;

  final acxtail = _twoDiffTail(ax, cx, acx);
  final bcxtail = _twoDiffTail(bx, cx, bcx);
  final acytail = _twoDiffTail(ay, cy, acy);
  final bcytail = _twoDiffTail(by, cy, bcy);

  if (acxtail == 0 && acytail == 0 && bcxtail == 0 && bcytail == 0) return det;

  final errbound = _ccwerrboundC * detsum + _resulterrbound * det.abs();
  double detAdj = det + (acx * bcytail + bcy * acxtail) - (acy * bcxtail + bcx * acytail);
  if (detAdj.abs() >= errbound) return detAdj;

  final u = _vec(4);
  final C1 = _vec(8);
  final C2 = _vec(12);
  final D = _vec(16);

  _crossProduct(acxtail, bcx, acytail, bcy, u);
  final C1len = _sum(4, B, 4, u, C1);

  _crossProduct(acx, bcxtail, acy, bcytail, u);
  final C2len = _sum(C1len, C1, 4, u, C2);

  _crossProduct(acxtail, bcxtail, acytail, bcytail, u);
  final Dlen = _sum(C2len, C2, 4, u, D);

  return D[Dlen - 1];
}

// Computes cross product into expansion
void _crossProduct(double ax, double bx, double ay, double by, Float64List out) {
  out[0] = ay * bx;
  out[1] = -(ax * by);
  out[2] = 0;
  out[3] = 0;
}

// Computes tail of a difference
double _twoDiffTail(double a, double b, double x) {
  double bvirt = a - x;
  double avirt = x + bvirt;
  double bround = bvirt - b;
  double around = a - avirt;
  return around + bround;
}
