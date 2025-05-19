import 'dart:math';
import 'package:concaveman/src/orient2d.dart';
import 'package:concaveman/src/point_in_polygon.dart';
import 'package:rbush/rbush.dart';

List<Point<double>> concaveman(List<Point<double>> points, [double concavity = 2, double lengthThreshold = 0]) {
  // a relative measure of concavity; higher value means simpler hull
  concavity = max(0, concavity);

  final sqConcavity = concavity * concavity;
  final sqLengthThreshold = lengthThreshold * lengthThreshold;

  final hull = _fastConvexHull(points);
  // index the points with an R-tree
  final pointTree = RBushBase<_PointBox>(
    maxEntries: 16,
    toBBox: (p) => p,
    getMinX: (p) => p.minX,
    getMinY: (p) => p.minY,
  )..load(points.map(_PointBox.new).toList());

  // index the segments with an R-tree (for intersection checks)
  final edgeTree = RBushBase<_EdgeBox>(
    maxEntries: 16,
    toBBox: (e) => e,
    getMinX: (e) => e.minX,
    getMinY: (e) => e.minY,
  );

  _LinkedNode? last;
  final queue = <_LinkedNode>[];
  for (final p in hull) {
    pointTree.remove(_PointBox(p));
    last = _insertNode(p, last);
    queue.add(last);
  }

  for (final node in queue) {
    edgeTree.insert(_EdgeBox(node));
  }

  // process edges one by one
  while (queue.isNotEmpty) {
    final node = queue.removeAt(0);
    final a = node.p;
    final b = node.next!.p;

    // skip the edge if it's already short enough
    final sqLen = _getSqDist(a, b);
    if (sqLen < sqLengthThreshold) continue;

    final maxSqLen = sqLen / sqConcavity;

    // find the best connection point for the current edge to flex inward to
    final candidate = _findCandidate(
      pointTree,
      node.prev!.p,
      a,
      b,
      node.next!.next!.p,
      maxSqLen,
      edgeTree,
    );

    // if we found a connection and it satisfies our concavity measure
    if (candidate != null && min(_getSqDist(candidate, a), _getSqDist(candidate, b)) <= maxSqLen) {
      // connect the edge endpoints through this point and add 2 new edges to the queue
      queue.add(node);
      queue.add(_insertNode(candidate, node));

      // update point and segment indexes
      pointTree.remove(_PointBox(candidate));
      edgeTree.remove(_EdgeBox(node));
      edgeTree.insert(_EdgeBox(node));
      edgeTree.insert(_EdgeBox(node.next!));
    }
  }

  // convert the resulting hull linked list to an array of points
  _LinkedNode? node = last;
  final concave = <Point<double>>[];
  do {
    concave.add(node!.p);
    node = node.next;
  } while (node != last);

  concave.add(node!.p);

  return concave;
}

// --- Core Structures ---

class _LinkedNode {
  _LinkedNode(this.p);

  final Point<double> p;
  _LinkedNode? prev;
  _LinkedNode? next;
}

class _PointBox extends RBushBox {
  _PointBox(this.point)
      : super(
          minX: point.x,
          minY: point.y,
          maxX: point.x,
          maxY: point.y,
        );

  final Point<double> point;
}

class _EdgeBox extends RBushBox {
  _EdgeBox(this.node)
      : super(
          minX: min(node.p.x, node.next!.p.x),
          minY: min(node.p.y, node.next!.p.y),
          maxX: max(node.p.x, node.next!.p.x),
          maxY: max(node.p.y, node.next!.p.y),
        );

  final _LinkedNode node;
}

// --- Utilities ---

_LinkedNode _insertNode(Point<double> p, _LinkedNode? prev) {
  final node = _LinkedNode(p);
  if (prev == null) {
    node.prev = node;
    node.next = node;
  } else {
    final next = prev.next!;
    node.next = next;
    node.prev = prev;
    prev.next = node;
    next.prev = node;
  }
  return node;
}

double _getSqDist(Point<double> a, Point<double> b) {
  final dx = a.x - b.x;
  final dy = a.y - b.y;
  return dx * dx + dy * dy;
}

bool _intersects(Point<double> p1, Point<double> q1, Point<double> p2, Point<double> q2) {
  if (p1 == q2 || q1 == p2) return false;
  final o1 = orient2d(p1, q1, p2);
  final o2 = orient2d(p1, q1, q2);
  final o3 = orient2d(p2, q2, p1);
  final o4 = orient2d(p2, q2, q1);
  return (o1 > 0) != (o2 > 0) && (o3 > 0) != (o4 > 0);
}

double _cross(Point<double> p1, Point<double> p2, Point<double> p3) {
  return orient2d(p1, p2, p3);
}

bool _noIntersections(Point<double> a, Point<double> b, RBushBase<_EdgeBox> segTree) {
  final searchBox = RBushBox(
    minX: min(a.x, b.x),
    minY: min(a.y, b.y),
    maxX: max(a.x, b.x),
    maxY: max(a.y, b.y),
  );
  final candidates = segTree.search(searchBox);

  for (final edge in candidates) {
    final p1 = edge.node.p;
    final p2 = edge.node.next!.p;
    if (_intersects(p1, p2, a, b)) return false;
  }

  return true;
}

Point<double>? _findCandidate(
  RBushBase<_PointBox> tree,
  Point<double> a,
  Point<double> b,
  Point<double> c,
  Point<double> d,
  double maxSqDist,
  RBushBase<_EdgeBox> segTree,
) {
  final queue = TinyQueue<_QueueItem>(
    [],
    (a, b) => a.dist.compareTo(b.dist),
  );

  void traverse(dynamic node) {
    if (node.leaf == true) {
      for (final _PointBox box in node.leafChildren) {
        final dist = _sqSegDist(box.point, b, c);
        if (dist <= maxSqDist) {
          queue.push(_QueueItem(point: box.point, dist: dist));
        }
      }
    } else {
      for (final RBushBox child in node.children) {
        final dist = _sqSegBoxDist(b, c, child);
        if (dist <= maxSqDist) {
          traverse(child);
        }
      }
    }
  }

  traverse(tree.data);

  while (queue.isNotEmpty) {
    final item = queue.pop();
    final p = item.point!;
    final d0 = _sqSegDist(p, a, b);
    final d1 = _sqSegDist(p, c, d);

    if (item.dist < d0 && item.dist < d1 && _noIntersections(b, p, segTree) && _noIntersections(c, p, segTree)) {
      return p;
    }
  }

  return null;
}

class _QueueItem {
  _QueueItem({this.point, required this.dist});

  final Point<double>? point;
  final double dist;
}

// Squared distance from point to segment
double _sqSegDist(Point<double> p, Point<double> p1, Point<double> p2) {
  final x = p1.x;
  final y = p1.y;
  final dx = p2.x - x;
  final dy = p2.y - y;

  double px = p.x, py = p.y;
  double t = 0;
  if (dx != 0 || dy != 0) {
    t = ((px - x) * dx + (py - y) * dy) / (dx * dx + dy * dy);
    if (t > 1) {
      px = p2.x;
      py = p2.y;
    } else if (t > 0) {
      px = x + dx * t;
      py = y + dy * t;
    } else {
      px = x;
      py = y;
    }
  }

  final dx1 = p.x - px;
  final dy1 = p.y - py;
  return dx1 * dx1 + dy1 * dy1;
}

// Squared distance from a segment to an RBushBox
double _sqSegBoxDist(Point<double> a, Point<double> b, RBushBox box) {
  if (_pointInBox(a, box) || _pointInBox(b, box)) return 0;
  final d1 = _sqSegSegDist(a, b, Point<double>(box.minX, box.minY), Point<double>(box.maxX, box.minY));
  if (d1 == 0) return 0;
  final d2 = _sqSegSegDist(a, b, Point<double>(box.minX, box.minY), Point<double>(box.minX, box.maxY));
  if (d2 == 0) return 0;
  final d3 = _sqSegSegDist(a, b, Point<double>(box.maxX, box.minY), Point<double>(box.maxX, box.maxY));
  if (d3 == 0) return 0;
  final d4 = _sqSegSegDist(a, b, Point<double>(box.minX, box.maxY), Point<double>(box.maxX, box.maxY));
  if (d4 == 0) return 0;
  return [d1, d2, d3, d4].reduce(min);
}

bool _pointInBox(Point<double> p, RBushBox b) {
  return p.x >= b.minX && p.x <= b.maxX && p.y >= b.minY && p.y <= b.maxY;
}

// Squared segment to segment distance
// Ported from geomalgorithms.com
double _sqSegSegDist(Point<double> p1, Point<double> p2, Point<double> q1, Point<double> q2) {
  final u = Point<double>(p2.x - p1.x, p2.y - p1.y);
  final v = Point<double>(q2.x - q1.x, q2.y - q1.y);
  final w = Point<double>(p1.x - q1.x, p1.y - q1.y);

  final a = u.x * u.x + u.y * u.y;
  final b = u.x * v.x + u.y * v.y;
  final c = v.x * v.x + v.y * v.y;
  final d = u.x * w.x + u.y * w.y;
  final e = v.x * w.x + v.y * w.y;
  final D = a * c - b * b;

  double sc = 0, sN = 0, sD = D;
  double tc = 0, tN = 0, tD = D;

  if (D == 0) {
    sN = 0;
    sD = 1;
    tN = e;
    tD = c;
  } else {
    sN = b * e - c * d;
    tN = a * e - b * d;
    if (sN < 0) {
      sN = 0;
      tN = e;
      tD = c;
    } else if (sN > sD) {
      sN = sD;
      tN = e + b;
      tD = c;
    }
  }

  if (tN < 0) {
    tN = 0;
    if (-d < 0) {
      sN = 0;
    } else if (-d > a) {
      sN = sD;
    } else {
      sN = -d;
      sD = a;
    }
  } else if (tN > tD) {
    tN = tD;
    if ((-d + b) < 0) {
      sN = 0;
    } else if ((-d + b) > a) {
      sN = sD;
    } else {
      sN = -d + b;
      sD = a;
    }
  }

  sc = sN == 0 ? 0 : sN / sD;
  tc = tN == 0 ? 0 : tN / tD;

  final dp = Point<double>(
    w.x + (sc * u.x) - (tc * v.x),
    w.y + (sc * u.y) - (tc * v.y),
  );

  return dp.x * dp.x + dp.y * dp.y;
}

List<Point<double>> _convexHull(List<Point<double>> points) {
  if (points.length <= 1) return List.from(points);
  final sorted = List<Point<double>>.from(points)..sort((a, b) => a.x == b.x ? a.y.compareTo(b.y) : a.x.compareTo(b.x));

  final lower = <Point<double>>[];
  for (final p in sorted) {
    while (lower.length >= 2 && _cross(lower[lower.length - 2], lower[lower.length - 1], p) <= 0) {
      lower.removeLast();
    }
    lower.add(p);
  }

  final upper = <Point<double>>[];
  for (int i = sorted.length - 1; i >= 0; i--) {
    final p = sorted[i];
    while (upper.length >= 2 && _cross(upper[upper.length - 2], upper[upper.length - 1], p) <= 0) {
      upper.removeLast();
    }
    upper.add(p);
  }

  upper.removeLast();
  lower.removeLast();
  return lower + upper;
}

List<Point<double>> _fastConvexHull(List<Point<double>> points) {
  var left = points[0];
  var right = points[0];
  var top = points[0];
  var bottom = points[0];

  for (final p in points) {
    if (p.x < left.x) left = p;
    if (p.x > right.x) right = p;
    if (p.y < top.y) top = p;
    if (p.y > bottom.y) bottom = p;
  }

  final cull = [left, top, right, bottom];
  final filtered = points.where((p) => !isPointInPolygon(p, cull)).toList();
  return _convexHull([...cull, ...filtered]);
}
