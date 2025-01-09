import 'dart:math';

import 'package:concaveman/concaveman.dart';

void main() {
  final points = <Point<double>>[
    Point(0, 0),
    Point(0, 1),
    Point(1, 1),
    Point(1, 0),
    Point(0.5, 0.5),
  ];
  print(concaveman(points));
}
