#ifndef lattice2grid_h_
#define lattice2grid_h_

#include <vector>
using namespace std;

class Vector2D {
public:
  Vector2D() {
  }
  ~Vector2D() {
  }
  Vector2D operator+(const Vector2D &op) const;
  Vector2D operator-(const Vector2D &op) const;
  Vector2D operator-() const;
  double operator*(const Vector2D &op) const;
  Vector2D operator*(const double &s) const;
  Vector2D operator/(const double &s) const;
  double crossProduct(const Vector2D &op) const;

  double x;
  double y;
};

class Color3f {
 public:
  Color3f() {
    color[0] = color[1] = color[2] = 0.85;
  }
  ~Color3f() {
  }
  float color[3];
};

class ColorPT {
 public:
  ColorPT() {
  }
  ~ColorPT() {
  }
  Color3f color;
  Vector2D pos;
};

struct ParaFI {
 public:
  ParaFI() {
    dimLatticeX = 0;
    dimLatticeY = 0;

    nCols = 128;
    nRows = 128;
    pixelGap = 0.3;
  }
  ~ParaFI() {
  }
  unsigned dimLatticeX;
  unsigned dimLatticeY;

  Vector2D lowerCorner;
  unsigned nCols;
  unsigned nRows;
  double pixelGap;
};

bool BilinearCoord(const Vector2D &query,
  const Vector2D &p00,
  const Vector2D &p01,
  const Vector2D &p10,
  const Vector2D &p11,
  double *x,
  double *y);

void FFD2Image(const vector<ColorPT> &lattice,
  const ParaFI &para,
  vector<Color3f> *image);

#endif