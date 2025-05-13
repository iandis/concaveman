// speed up convex hull by filtering out points inside quadrilateral formed by 4 extreme points
import 'dart:math';

import 'package:concaveman/src/orient2d.dart';
import 'package:concaveman/src/point_in_polygon.dart';
import 'package:rbush/rbush.dart';

List<Point<double>> concaveman(List<Point<double>> points, [double? concavity, double? lengthThreshold]) {
  // a relative measure of concavity; higher value means simpler hull
  concavity = max(0, concavity ?? 2);

  // when a segment goes below this length threshold, it won't be drilled down further
  lengthThreshold = lengthThreshold ?? 0;

  // start with a convex hull of the points
  var hull = fastConvexHull(points);

  // index the points with an R-tree
  final tree = RBushPoint(16);
  tree.load(points);

  // turn the convex hull into a linked list and populate the initial edge queue with the nodes
  var queue = <SegPoint>[];
  SegPoint? last;
  for (var i = 0; i < hull.length; i++) {
    final p = hull[i];
    tree.remove(p);
    last = insertNode(p, last);
    queue.add(last);
  }

  // index the segments with an R-tree (for intersection checks)
  var segTree = RBushSegPoint(16);
  for (int i = 0; i < queue.length; i++) {
    segTree.insert(updateBBox(queue[i]));
  }

  var sqConcavity = concavity * concavity;
  var sqLenThreshold = lengthThreshold * lengthThreshold;

  // process edges one by one
  while (queue.isNotEmpty) {
    var node = queue.removeAt(0);
    var a = node.p;
    var b = node.next!.p;

    // skip the edge if it's already short enough
    var sqLen = getSqDist(a, b);
    if (sqLen < sqLenThreshold) continue;

    var maxSqLen = sqLen / sqConcavity;

    // find the best connection point for the current edge to flex inward to
    final p = findCandidate(tree, node.prev!.p, a, b, node.next!.next!.p, maxSqLen, segTree);

    // if we found a connection and it satisfies our concavity measure
    if (p != null && min(getSqDist(p, a), getSqDist(p, b)) <= maxSqLen) {
      // connect the edge endpoints through this point and add 2 new edges to the queue
      queue.add(node);
      queue.add(insertNode(p, node));

      // update point and segment indexes
      tree.remove(p);
      segTree.remove(node);
      segTree.insert(updateBBox(node));
      segTree.insert(updateBBox(node.next!));
    }
  }

  // convert the resulting hull linked list to an array of points
  SegPoint? node = last;
  var concave = <Point<double>>[];
  do {
    concave.add(node!.p);
    node = node.next;
  } while (node != last);

  concave.add(node!.p);

  return concave;
}

List<Point<double>> fastConvexHull(List<Point<double>> points) {
  var left = points[0];
  var top = points[0];
  var right = points[0];
  var bottom = points[0];

  // find the leftmost, rightmost, topmost and bottommost points
  for (var i = 0; i < points.length; i++) {
    var p = points[i];
    if (p.x < left.x) left = p;
    if (p.x > right.x) right = p;
    if (p.y < top.y) top = p;
    if (p.y > bottom.y) bottom = p;
  }

  // filter out points that are inside the resulting quadrilateral
  var cull = [left, top, right, bottom];
  var filtered = cull.toList();
  for (var i = 0; i < points.length; i++) {
    if (!isPointInPolygon(points[i], cull)) filtered.add(points[i]);
  }

  // get convex hull around the filtered points
  return convexHull(filtered);
}

// TODO: check again
int compareByX(Point<double> a, Point<double> b) {
  return a.x == b.x ? a.y.compareTo(b.y) : a.x.compareTo(b.x);
}

List<Point<double>> convexHull(List<Point<double>> points) {
  points.sort(compareByX);

  final lower = <Point<double>>[];
  for (var i = 0; i < points.length; i++) {
    while (lower.length >= 2 && cross(lower[lower.length - 2], lower[lower.length - 1], points[i]) <= 0) {
      lower.removeLast();
    }
    lower.add(points[i]);
  }

  final upper = <Point<double>>[];
  for (var ii = points.length - 1; ii >= 0; ii--) {
    while (upper.length >= 2 && cross(upper[upper.length - 2], upper[upper.length - 1], points[ii]) <= 0) {
      upper.removeLast();
    }
    upper.add(points[ii]);
  }

  upper.removeLast();
  lower.removeLast();
  return lower + upper;
}

double cross(Point<double> p1, Point<double> p2, Point<double> p3) {
  return orient2d(p1, p2, p3);
}

class _QueueItem {
  const _QueueItem(this.node, this.dist);
  final dynamic node;
  final double dist;
}

Point<double>? findCandidate(
  RBushPoint tree,
  Point<double> a,
  Point<double> b,
  Point<double> c,
  Point<double> d,
  double maxDist,
  RBushSegPoint segTree,
) {
  var queue = TinyQueue(<_QueueItem>[], _compareDist);
  var node = tree.data;
  // search through the point R-tree with a depth-first search using a priority queue
  // in the order of distance to the edge (b, c)
  while (true) {
    for (var i = 0; i < node.children.length; i++) {
      final child = node.leaf ? node.leafChildren[i] : node.children[i];
      var dist = node.leaf ? sqSegDist(child as Point<double>, b, c) : sqSegBoxDist(b, c, child as RBushBox);
      if (dist > maxDist) continue; // skip the node if it's farther than we ever need

      queue.push(_QueueItem(child, dist));
    }

    while (queue.isNotEmpty && queue.peek().node is Point<double>) {
      var item = queue.pop();
      var p = item.node as Point<double>;

      // skip all points that are as close to adjacent edges (a,b) and (c,d),
      // and points that would introduce self-intersections when connected
      var d0 = sqSegDist(p, a, b);
      var d1 = sqSegDist(p, c, d);
      if (item.dist < d0 && item.dist < d1 && noIntersections(b, p, segTree) && noIntersections(c, p, segTree)) {
        return p;
      }
    }

    if (queue.isEmpty) break;

    node = queue.pop().node;
  }

  return null;
}

// segment to segment distance, ported from http://geomalgorithms.com/a07-_distance.html by Dan Sunday
double sqSegSegDist(double x0, double y0, double x1, double y1, double x2, double y2, double x3, double y3) {
  var ux = x1 - x0;
  var uy = y1 - y0;
  var vx = x3 - x2;
  var vy = y3 - y2;
  var wx = x0 - x2;
  var wy = y0 - y2;
  var a = ux * ux + uy * uy;
  var b = ux * vx + uy * vy;
  var c = vx * vx + vy * vy;
  var d = ux * wx + uy * wy;
  var e = vx * wx + vy * wy;
  var D = a * c - b * b;

  double sc, sN, tc, tN;
  var sD = D;
  var tD = D;

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

  if (tN < 0.0) {
    tN = 0.0;
    if (-d < 0.0) {
      sN = 0.0;
    } else if (-d > a)
      sN = sD;
    else {
      sN = -d;
      sD = a;
    }
  } else if (tN > tD) {
    tN = tD;
    if ((-d + b) < 0.0) {
      sN = 0;
    } else if (-d + b > a)
      sN = sD;
    else {
      sN = -d + b;
      sD = a;
    }
  }

  sc = sN == 0 ? 0 : sN / sD;
  tc = tN == 0 ? 0 : tN / tD;

  var cx = (1 - sc) * x0 + sc * x1;
  var cy = (1 - sc) * y0 + sc * y1;
  var cx2 = (1 - tc) * x2 + tc * x3;
  var cy2 = (1 - tc) * y2 + tc * y3;
  var dx = cx2 - cx;
  var dy = cy2 - cy;

  return dx * dx + dy * dy;
}

int _compareDist(_QueueItem a, _QueueItem b) {
  return a.dist.compareTo(b.dist);
}

// square distance from a segment bounding box to the given one
double sqSegBoxDist(Point<double> a, Point<double> b, RBushBox bbox) {
  if (inside(a, bbox) || inside(b, bbox)) return 0;
  var d1 = sqSegSegDist(a.x, a.y, b.x, b.y, bbox.minX, bbox.minY, bbox.maxX, bbox.minY);
  if (d1 == 0) return 0;
  var d2 = sqSegSegDist(a.x, a.y, b.x, b.y, bbox.minX, bbox.minY, bbox.minX, bbox.maxY);
  if (d2 == 0) return 0;
  var d3 = sqSegSegDist(a.x, a.y, b.x, b.y, bbox.maxX, bbox.minY, bbox.maxX, bbox.maxY);
  if (d3 == 0) return 0;
  var d4 = sqSegSegDist(a.x, a.y, b.x, b.y, bbox.minX, bbox.maxY, bbox.maxX, bbox.maxY);
  if (d4 == 0) return 0;
  // TODO: check again
  return min(
    min(min(d1, d2), d3),
    d4,
  );
}

bool inside(Point<double> a, RBushBox bbox) {
  return a.x >= bbox.minX && a.x <= bbox.maxX && a.y >= bbox.minY && a.y <= bbox.maxY;
}

class SegPoint {
  SegPoint({required this.p, this.prev, this.next});
  Point<double> p;
  SegPoint? prev;
  SegPoint? next;
  RBushBox bbox = RBushBox(minX: 0, minY: 0, maxX: 0, maxY: 0);
}

// check if the edge (a,b) doesn't intersect any other edges
bool noIntersections(Point<double> a, Point<double> b, RBushSegPoint segTree) {
  var minX = min(a.x, b.x);
  var minY = min(a.y, b.y);
  var maxX = max(a.x, b.x);
  var maxY = max(a.y, b.y);

  var edges = segTree.search(RBushBox(minX: minX, minY: minY, maxX: maxX, maxY: maxY));
  for (var i = 0; i < edges.length; i++) {
    if (intersects(edges[i].p, edges[i].next!.p, a, b)) return false;
  }
  return true;
}

// check if the edges (p1,q1) and (p2,q2) intersect
bool intersects(Point<double> p1, Point<double> q1, Point<double> p2, Point<double> q2) {
  return p1 != q2 &&
      q1 != p2 &&
      cross(p1, q1, p2) > 0 != cross(p1, q1, q2) > 0 &&
      cross(p2, q2, p1) > 0 != cross(p2, q2, q1) > 0;
}

// update the bounding box of a node's edge
SegPoint updateBBox(SegPoint node) {
  var p1 = node.p;
  var p2 = node.next!.p;
  node.bbox
    ..minX = min(p1.x, p2.x)
    ..minY = min(p1.y, p2.y)
    ..maxX = max(p1.x, p2.x)
    ..maxY = max(p1.y, p2.y);
  return node;
}

// create a new node in a doubly linked list
SegPoint insertNode(Point<double> p, SegPoint? prev) {
  var node = SegPoint(p: p, prev: null, next: null);

  if (prev == null) {
    node.prev = node;
    node.next = node;
  } else {
    node.next = prev.next;
    node.prev = prev;
    prev.next?.prev = node;
    prev.next = node;
  }
  return node;
}

// square distance between 2 points
double getSqDist(Point<double> p1, Point<double> p2) {
  var dx = p1.x - p2.x, dy = p1.y - p2.y;

  return dx * dx + dy * dy;
}

// square distance from a point to a segment
double sqSegDist(Point<double> p, Point<double> p1, Point<double> p2) {
  var x = p1.x, y = p1.y, dx = p2.x - x, dy = p2.y - y;

  if (dx != 0 || dy != 0) {
    var t = ((p.x - x) * dx + (p.y - y) * dy) / (dx * dx + dy * dy);

    if (t > 1) {
      x = p2.x;
      y = p2.y;
    } else if (t > 0) {
      x += dx * t;
      y += dy * t;
    }
  }

  dx = p.x - x;
  dy = p.y - y;

  return dx * dx + dy * dy;
}

class RBushPoint extends RBushBase<Point<double>> {
  RBushPoint(int maxEntries)
      : super(
          maxEntries: maxEntries,
          toBBox: (item) {
            return RBushBox(
              minX: item.x,
              minY: item.y,
              maxX: item.x,
              maxY: item.y,
            );
          },
          getMinX: (item) => item.x,
          getMinY: (item) => item.y,
        );
}

class RBushSegPoint extends RBushBase<SegPoint> {
  RBushSegPoint(int maxEntries)
      : super(
          maxEntries: maxEntries,
          toBBox: (item) {
            final point = item.p;
            return RBushBox(
              minX: point.x,
              minY: point.y,
              maxX: point.x,
              maxY: point.y,
            );
          },
          getMinX: (item) => item.bbox.minX,
          getMinY: (item) => item.bbox.minY,
        );
}
