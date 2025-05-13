import 'dart:math';

bool isPointInPolygon(Point<double> point, List<Point<double>> vs, [int start = 0, int? end]) {
  final Point(:x, :y) = point;
  var inside = false;

  end ??= vs.length;
  var len = end - start;
  for (var i = 0, j = len - 1; i < len; j = i++) {
    final xi = vs[i + start].x;
    final yi = vs[i + start].y;

    final xj = vs[j + start].x;
    final yj = vs[j + start].y;

    final intersect = ((yi > y) != (yj > y)) && (x < (xj - xi) * (y - yi) / (yj - yi) + xi);
    if (intersect) inside = !inside;
  }

  return inside;
}
