#include "lattice2grid.h"

Vector2D Vector2D::operator+(const Vector2D &op) const {
  Vector2D result;
  result.x = x + op.x;
  result.y = y + op.y;
  return result;
}

Vector2D Vector2D::operator-(const Vector2D &op) const {
  Vector2D result;
  result.x = x - op.x;
  result.y = y - op.y;
  return result;
}

Vector2D Vector2D::operator-() const {
  Vector2D result;
  result.x = -x;
  result.y = -y;
  return result;
}

double Vector2D::operator*(const Vector2D &op) const {
  return op.x*x + op.y*y;
}

Vector2D Vector2D::operator*(const double &s) const {
  Vector2D result;
  result.x = x*s;
  result.y = y*s;
  return result;
}


Vector2D Vector2D::operator/(const double &s) const {
  Vector2D result;
  result.x = x/s;
  result.y = y/s;
  return result;
}

double Vector2D::crossProduct(const Vector2D &op) const {
  return x*op.y - y*op.x;
}

bool BilinearCoord(const Vector2D &query,
  const Vector2D &p00,
  const Vector2D &p01,
  const Vector2D &p10,
  const Vector2D &p11,
  double *x,
  double *y) {
  Vector2D ex0 = p10 - p00;
  Vector2D ex1 = p11 - p01;
  Vector2D qx0 = query - p00;
  Vector2D qx1 = query - p01;
  double a = ex0.crossProduct(ex1);
  double b = ex0.crossProduct(qx1) + qx0.crossProduct(ex1);
  double c = qx0.crossProduct(qx1);

  double d = b*b - 4*a*c;
  if (d < 0.0)
    return false;

  double d1 = (b - sqrt(d))/2;
  double d2 = (b + sqrt(d))/2;
  *x = 0.0;
  if (c >= 1e-13) {
    if (max(d1, d2)<= 0.999999*c)
      return false;
    *x = c/max(d1, d2);
  } else if (c <= -1e-13) {
    if (min(d1, d2) >= 0.999999*c)
      return false;
    *x = c/min(d1, d2);
  } else {
    if (max(fabs(d1), fabs(d2)) <= 0.999999*fabs(c))
      return false;
    *x = fabs(c)/max(fabs(d1), fabs(d2));
  }

  Vector2D q0y = p00*(1-(*x)) + p10*(*x);
  Vector2D q1y = p01*(1-(*x)) + p11*(*x);

  Vector2D qy = query - q0y;
  Vector2D q10 = q1y - q0y;

  *y = (qy*q10)/(q10*q10);
  if (*y < 1.000001 && *y > -0.000001)
    return true;
  else
    return false;
}

void FFD2Image(const vector<ColorPT> &lattice,
  const ParaFI &para,
  vector<Color3f> *image) {
  unsigned numPixels = para.nCols*para.nRows;
  image->resize(numPixels);

  for (unsigned i = 0; i < para.dimLatticeX-1; ++i) {
    for (unsigned j = 0; j < para.dimLatticeY-1; ++j) {
      unsigned id00 = i*para.dimLatticeY + j;
      unsigned id01 = i*para.dimLatticeY + j+1;
      unsigned id10 = (i+1)*para.dimLatticeY + j;
      unsigned id11 = (i+1)*para.dimLatticeY + j+1;

      const ColorPT &p00 = lattice[id00];
      const ColorPT &p01 = lattice[id01];
      const ColorPT &p10 = lattice[id10];
      const ColorPT &p11 = lattice[id11];

      double lowerx = min(min(p00.pos.x, p01.pos.x), min(p10.pos.x, p11.pos.x)) - 1e-6;
      double lowery = min(min(p00.pos.y, p01.pos.y), min(p10.pos.y, p11.pos.y)) - 1e-6;
      double upperx = max(max(p00.pos.x, p01.pos.x), max(p10.pos.x, p11.pos.x)) + 1e-6;
      double uppery = max(max(p00.pos.y, p01.pos.y), max(p10.pos.y, p11.pos.y)) + 1e-6;

      int x_left = int((lowerx - para.lowerCorner.x)/para.pixelGap - 0.5);
      int y_left = int((lowery - para.lowerCorner.y)/para.pixelGap - 0.5);
      
      int x_right = int((upperx - para.lowerCorner.x)/para.pixelGap - 0.5);
      int y_right = int((uppery - para.lowerCorner.y)/para.pixelGap - 0.5);

      for (int i = max(0, x_left); i <= min(int(para.nCols)-1, x_right); ++i) {
        for (int j = max(0, y_left); j <= min(int(para.nRows)-1, y_right); ++j) {
          Vector2D pixelPos;
          pixelPos.x = para.lowerCorner.x + (i+0.5)*para.pixelGap;
          pixelPos.y = para.lowerCorner.y + (j+0.5)*para.pixelGap;
          if (j == 99 && i== 0)
            int hqx = 10;
          double x, y;
          if (BilinearCoord(pixelPos, p00.pos, p01.pos, p10.pos, p11.pos, &x, &y)) {
            int pixelId = i*para.nRows + j;
            Color3f *c = &(*image)[pixelId];
            for (int i = 0; i < 3; ++i) {
              c->color[i] = (p00.color.color[i]*(1-x) + p10.color.color[i]*x)*(1-y)
                + (p01.color.color[i]*(1-x) + p11.color.color[i]*x)*y;
            }
          }
        }
      }
    }
  }
}