import 'dart:math';
import 'dart:typed_data';

const epsilon = 1.1102230246251565e-16;
const splitter = 134217729;
const resulterrbound = (3 + 8 * epsilon) * epsilon;

const ccwerrboundA = (3 + 16 * epsilon) * epsilon;
const ccwerrboundB = (2 + 12 * epsilon) * epsilon;
const ccwerrboundC = (9 + 64 * epsilon) * epsilon * epsilon;

double _orient2dadapt(Point<double> a, Point<double> b, Point<double> c, double detsum) {
  final acx = a.x - c.x;
  final bcx = b.x - c.x;
  final acy = a.y - c.y;
  final bcy = b.y - c.y;

  final C1 = Float64List(8);
  final C2 = Float64List(12);
  final D = Float64List(16);

  final B = _crossProduct(acx, acy, bcx, bcy);

  double det = _estimate(4, B);
  double errbound = ccwerrboundB * detsum;
  if (det >= errbound || -det >= errbound) {
    return det;
  }

  final acxtail = _twoDiffTail(a.x, c.x, acx);
  final bcxtail = _twoDiffTail(b.x, c.x, bcx);
  final acytail = _twoDiffTail(a.y, c.y, acy);
  final bcytail = _twoDiffTail(b.y, c.y, bcy);

  if (acxtail == 0 && acytail == 0 && bcxtail == 0 && bcytail == 0) {
    return det;
  }

  errbound = ccwerrboundC * detsum + resulterrbound * det.abs();
  det += (acx * bcytail + bcy * acxtail) - (acy * bcxtail + bcx * acytail);
  if (det >= errbound || -det >= errbound) return det;

  final u1 = _crossProduct(acxtail, acytail, bcx, bcy);
  final C1len = _sum(4, B, 4, u1, C1);

  final u2 = _crossProduct(acx, acy, bcxtail, bcytail);
  final C2len = _sum(C1len, C1, 4, u2, C2);

  final u3 = _crossProduct(acxtail, acytail, bcxtail, bcytail);
  final Dlen = _sum(C2len, C2, 4, u3, D);

  return D[Dlen - 1];
}

double orient2d(Point<double> a, Point<double> b, Point<double> c) {
  final detleft = (a.y - c.y) * (b.x - c.x);
  final detright = (a.x - c.x) * (b.y - c.y);
  final det = detleft - detright;

  final detsum = (detleft + detright).abs();
  if (det.abs() >= ccwerrboundA * detsum) return det;

  return -_orient2dadapt(a, b, c, detsum);
}

double _orient2dfast(Point<double> a, Point<double> b, Point<double> c) {
  return (a.y - c.y) * (b.x - c.x) - (a.x - c.x) * (b.y - c.y);
}

double _estimate(int elen, Float64List e) {
  var Q = e[0];
  for (int i = 1; i < elen; i++) {
    Q += e[i];
  }
  return Q;
}

Float64List _crossProduct(double acx, double acy, double bcx, double bcy) {
  final vector = Float64List(4);
  var bvirt, c, ahi, alo, bhi, blo, _i, _j, _0, s1, s0, t1, t0, u3;
  s1 = acx * bcy;

  c = splitter * acx;

  ahi = c - (c - acx);

  alo = acx - ahi;

  c = splitter * bcy;

  bhi = c - (c - bcy);

  blo = bcy - bhi;

  s0 = alo * blo - (s1 - ahi * bhi - alo * bhi - ahi * blo);

  t1 = acy * bcx;

  c = splitter * acy;

  ahi = c - (c - acy);

  alo = acy - ahi;

  c = splitter * bcx;

  bhi = c - (c - bcx);

  blo = bcx - bhi;

  t0 = alo * blo - (t1 - ahi * bhi - alo * bhi - ahi * blo);

  _i = s0 - t0;

  bvirt = s0 - _i;

  vector[0] = s0 - (_i + bvirt) + (bvirt - t0);

  _j = s1 + _i;

  bvirt = _j - s1;

  _0 = s1 - (_j - bvirt) + (_i - bvirt);

  _i = _0 - t1;

  bvirt = _0 - _i;

  vector[1] = _0 - (_i + bvirt) + (bvirt - t1);

  u3 = _j + _i;

  bvirt = u3 - _j;

  vector[2] = _j - (u3 - bvirt) + (_i - bvirt);

  vector[3] = u3;
  return vector;
}

// fast_expansion_sum_zeroelim routine from oritinal code
int _sum(int elen, Float64List e, int flen, Float64List f, Float64List h) {
  double Q;
  double Qnew;
  double hh;
  double bvirt;

  double enow = e[0];

  double fnow = f[0];

  int eindex = 0;

  int findex = 0;

  if ((fnow > enow) == (fnow > -enow)) {
    Q = enow;

    enow = e[++eindex];
  } else {
    Q = fnow;

    fnow = f[++findex];
  }

  int hindex = 0;

  if (eindex < elen && findex < flen) {
    if ((fnow > enow) == (fnow > -enow)) {
      Qnew = enow + Q;

      hh = Q - (Qnew - enow);

      enow = e[++eindex];
    } else {
      Qnew = fnow + Q;

      hh = Q - (Qnew - fnow);

      fnow = f[++findex];
    }

    Q = Qnew;

    if (hh != 0) {
      h[hindex++] = hh;
    }

    while (eindex < elen && findex < flen) {
      if ((fnow > enow) == (fnow > -enow)) {
        Qnew = Q + enow;

        bvirt = Qnew - Q;

        hh = Q - (Qnew - bvirt) + (enow - bvirt);

        enow = e[++eindex];
      } else {
        Qnew = Q + fnow;

        bvirt = Qnew - Q;

        hh = Q - (Qnew - bvirt) + (fnow - bvirt);

        fnow = f[++findex];
      }

      Q = Qnew;

      if (hh != 0) {
        h[hindex++] = hh;
      }
    }
  }

  while (eindex < elen) {
    Qnew = Q + enow;

    bvirt = Qnew - Q;

    hh = Q - (Qnew - bvirt) + (enow - bvirt);

    enow = e[++eindex];

    Q = Qnew;

    if (hh != 0) {
      h[hindex++] = hh;
    }
  }

  while (findex < flen) {
    Qnew = Q + fnow;

    bvirt = Qnew - Q;

    hh = Q - (Qnew - bvirt) + (fnow - bvirt);

    fnow = f[++findex];

    Q = Qnew;

    if (hh != 0) {
      h[hindex++] = hh;
    }
  }

  if (Q != 0 || hindex == 0) {
    h[hindex++] = Q;
  }

  return hindex;
}

double _twoDiffTail(double a, double b, double x) {
  double bvirt = a - x;
  return a - (x + bvirt) + (bvirt - b);
}
