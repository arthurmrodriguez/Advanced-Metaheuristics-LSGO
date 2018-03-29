#include "Circle.h"

Circle::Circle(int _radius, int _width) : radius(_radius), width(_width) {
  int dim = 2*radius + 1;

  map.resize(dim);

  for (int i=0; i<dim; i++) {
    map[i].resize(dim);
    for (int j=0; j<dim; j++)
      map[i][j] = 0;
  }

  circleMidpoint(radius, radius, radius);
  circleMidpoint(radius, radius, radius-width);
}

Circle::~Circle() {
}

void Circle::circlePoints(int cx, int cy, int x, int y) {
  if (x==0) {
    map[cx][cy+y] += 1;
    map[cx][cy-y] += 1;
    map[cx+y][cy] += 1;
    map[cx-y][cy] += 1;
  }
  else if (x==y) {
    map[cx+x][cy+y] += 1;
    map[cx-x][cy+y] += 1;
    map[cx+x][cy-y] += 1;
    map[cx-x][cy-y] += 1;
  }
  else  if (x<y) {
    map[cx+x][cy+y] += 1;
    map[cx-x][cy+y] += 1;
    map[cx+x][cy-y] += 1;
    map[cx-x][cy-y] += 1;
    map[cx+y][cy+x] += 1;
    map[cx-y][cy+x] += 1;
    map[cx+y][cy-x] += 1;
    map[cx-y][cy-x] += 1;
  }
}

void Circle::circleMidpoint(int xCenter, int yCenter, int radius) {
  int x = 0;
  int y = radius;
  int p = (5 - radius*4)/4;

  circlePoints(xCenter, yCenter, x, y);

  while (x < y) {
    x++;

    if (p < 0)
      p += 2*x+1;
    else {
      y--;
      p += 2*(x-y)+1;
    }

    circlePoints(xCenter, yCenter, x, y);
  }
}
