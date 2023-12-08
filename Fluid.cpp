#include "Fluid.h"
#include <cmath>
#include <iomanip>
#include <iostream> // debug!!

using namespace std;

Fluid::Fluid() {
  _pDensity = _pOldDensity = _pDensityStorage = _pCurl = nullptr;
  _pColor = nullptr;
  _pVelocity = _pOldVelocity = _pVelocityStorage = nullptr;
  _invTrueWidth = _invTrueHeight = _invWidth = _invHeight = _invNumCells = _dt =
      0.0;
  _width = _height = _numCells = _trueWidth = _trueHeight = _drawNumCells = 0;
  _diffusion = 0.002;
  _densityFade = 1;
  _viscosity = 1;
  _solverIterations = 10;
}

Fluid &Fluid::setSize(int xSize, int ySize, int scale) {
  _width = xSize;
  _height = ySize;
  _drawNumCells = _width * _height;
  _drawScale = scale;

  _invWidth = 1.0 / _width;
  _invHeight = 1.0 / _height;

  _trueWidth = _width + 2;
  _trueHeight = _height + 2;

  _invTrueWidth = 1.0 / _trueWidth;
  _invTrueHeight = 1.0 / _trueHeight;

  _numCells = _trueWidth * _trueHeight;
  _invNumCells = 1.0f / _numCells;

  //   _densityFade /= xSize;

  reinitialize();

  return *this;
}

void Fluid::reinitialize() {
  clear();

  _pDensity = new double[_numCells];
  _pOldDensity = new double[_numCells];
  _pDensityStorage = new double[_numCells];
  _pColor = new sf::Color[_numCells];
  _pVelocity = new sf::Vector2f[_numCells];
  _pOldVelocity = new sf::Vector2f[_numCells];
  _pVelocityStorage = new sf::Vector2f[_numCells];
  _pCurl = new double[_numCells];

  for (int i = 0; i < _numCells; i++) {
    _pDensity[i] = 0.0f;
    _pOldDensity[i] = 0.0f;
    _pColor[i] = sf::Color::Black;
    _pVelocity[i] = sf::Vector2f(0.0f, 0.0f);
    _pOldVelocity[i] = sf::Vector2f(0.0f, 0.0f);
    _pCurl[i] = 0.0f;
  }
  _image.create(_width * _drawScale, _height * _drawScale);
  _tex.loadFromImage(_image);
  _sprite.setTexture(_tex);
}

void Fluid::clear() {
  if (_pDensity)
    delete[] _pDensity;
  if (_pOldDensity)
    delete[] _pOldDensity;
  if (_pColor)
    delete[] _pColor;
  if (_pVelocity)
    delete[] _pVelocity;
  if (_pOldVelocity)
    delete[] _pOldVelocity;
  if (_pCurl)
    delete[] _pCurl;
}

uint16_t Fluid::index(const uint16_t X, const uint16_t Y) const {
  return Y * _trueWidth + X;
}

void Fluid::color() {
  double density;
  for (int i = 0; i < _numCells; i++) {
    density = _pDensity[i];
    _pColor[i] = sf::Color(density * 25.25 * 6, density * 25.25 * 8,
                           density * 25.25 * 12);
  }
}

void Fluid::draw(sf::RenderWindow &window) {
  sf::Color cellColor;
  for (uint16_t y = 0; y < _height; y++) {
    for (uint16_t x = 0; x < _width; x++) {
      cellColor = _pColor[index(x + 1, y + 1)];
      for (uint16_t i = 0; i < _drawScale; i++) {
        for (uint16_t j = 0; j < _drawScale; j++) {
          _image.setPixel((x * _drawScale) + j, (y * _drawScale) + i,
                          cellColor);
        }
      }
    }
    _tex.loadFromImage(_image);
    window.draw(_sprite);
  }
}

void Fluid::addVelocity(const sf::Vector2f velocity, const uint16_t x,
                        const uint16_t y) {
  _pVelocity[index(x, y)] += velocity;
}

void Fluid::addDensity(const double density, const uint16_t windowX,
                       const uint16_t windowY) {
  uint16_t x = (windowX / _drawScale) + 1;
  uint16_t y = (windowY / _drawScale) + 1;
  _pDensity[index(x, y)] = clamp(_pDensity[index(x, y)] + density, 0., 1.);
}

void Fluid::curl(const uint16_t x, const uint16_t y) {
  double du_dy = _pVelocity[index(x, y + 1)].x - _pVelocity[index(x, y - 1)].x;
  double dv_dx = _pVelocity[index(y + 1, x)].y - _pVelocity[index(y - 1, x)].y;
  _pCurl[index(x, y)] = (du_dy - dv_dx) * 0.5f;
}

float clamp(const float VALUE, const float MIN, const float MAX) {
  if (VALUE < MIN)
    return MIN;
  if (VALUE > MAX)
    return MAX;
  return VALUE;
}

void Fluid::setDensityBoundary() {
  int destination1Index, destination2Index, source1Index, source2Index;
  int step = _trueWidth;

  destination1Index = index(0, 1);
  source1Index = index(1, 1);
  destination2Index = index(_width + 1, 1);
  source2Index = index(_width, 1);
  for (int i = 1; i <= _height; i++) {
    _pDensity[destination1Index] = _pDensity[source1Index];
    destination1Index += step;
    source1Index += step;
    _pDensity[destination2Index] = _pDensity[source2Index];
    destination2Index += step;
    source2Index += step;
  }

  destination1Index = index(1, 0);
  source1Index = index(1, 1);
  destination2Index = index(1, _height + 1);
  source2Index = index(1, _height);
  for (int i = 1; i <= _width; i++) {
    _pDensity[destination1Index++] = _pDensity[source1Index++];
    _pDensity[destination2Index++] = _pDensity[source2Index++];
  }

  _pDensity[index(0, 0)] =
      0.5f * (_pDensity[index(1, 0)] + _pDensity[index(0, 1)]);
  _pDensity[index(0, _height + 1)] =
      0.5f * (_pDensity[index(1, _height + 1)] + _pDensity[index(0, _height)]);
  _pDensity[index(_width + 1, 0)] =
      0.5f * (_pDensity[index(_width, 0)] + _pDensity[index(_width + 1, 1)]);
  _pDensity[index(_width + 1, _height + 1)] =
      0.5f * (_pDensity[index(_width, _height + 1)] +
              _pDensity[index(_width + 1, _height)]);
}

void Fluid::advectDensity() {
  int xVectorIndex, yVectorIndex, rightOfXVectorIndex, belowYVectorIndex;
  float xVector, yVector, decimalBetweenXVectors, decimalBetweenYVectors,
      betweenXVectors, betweenYVectors;
  int linearIndex;

  const float xdt = _dt * _width;
  const float ydt = _dt * _height;

  for (int y = 1; y <= _height; y++) {
    for (int x = 1; x <= _width; x++) {
      linearIndex = index(x, y);
      xVector = x - xdt * _pOldVelocity[linearIndex].x;
      yVector = y - ydt * _pOldVelocity[linearIndex].y;

      xVector = clamp(xVector, 0.5f, _width + 0.5);

      xVectorIndex = (int)xVector;
      rightOfXVectorIndex = xVectorIndex + 1;

      yVector = clamp(yVector, 0.5f, _height + 0.5);

      yVectorIndex = (int)yVector;
      belowYVectorIndex = yVectorIndex + 1;

      betweenXVectors = xVector - xVectorIndex;
      decimalBetweenXVectors = 1 - betweenXVectors;
      betweenYVectors = yVector - yVectorIndex;
      decimalBetweenYVectors = 1 - betweenYVectors;

      _pDensity[linearIndex] =
          decimalBetweenXVectors *
              (decimalBetweenYVectors *
                   _pOldDensity[index(xVectorIndex, yVectorIndex)] +
               betweenYVectors *
                   _pOldDensity[index(xVectorIndex, belowYVectorIndex)]) +
          betweenXVectors *
              (decimalBetweenYVectors *
                   _pOldDensity[index(rightOfXVectorIndex, yVectorIndex)] +
               betweenYVectors *
                   _pOldDensity[index(rightOfXVectorIndex, belowYVectorIndex)]);
    }
  }
  setDensityBoundary();
}

void Fluid::diffuseDensity() {
  double a = _dt * _diffusion * _width * _height;
  int linearIndex;
  double c = 1. / (1. + 4. * a);
  for (int iteration = 0; iteration < _solverIterations; iteration++) {
    for (int y = 1; y <= _height; y++) {
      linearIndex = index(1, y);
      for (int x = 1; x <= _width; x++) {
        _pDensity[linearIndex] =
            ((_pDensity[linearIndex - 1] + _pDensity[linearIndex + 1] +
              _pDensity[linearIndex - _trueWidth] +
              _pDensity[linearIndex + _trueWidth]) *
                 a +
             _pOldDensity[linearIndex]) *
            c;
        linearIndex++; // index(x,y)
      }
    }
    setDensityBoundary();
  }
}

void Fluid::fadeDensity() {
  for (int i = 0; i < _numCells; i++) {
    // _pDensity[i] *= _densityFade;
    _pOldDensity[i] = 0;
    _pOldVelocity[i] = sf::Vector2f();
  }
}

void Fluid::project() {
  // change divergence to curl I think
  int linearIndex;
  float h = -0.5f / _width;
  for (int y = 1; y <= _height; y++) {
    linearIndex = index(1, y);
    for (int x = 1; x <= _width; x++) {
      _pOldVelocity[linearIndex].x =
          h * (_pVelocity[linearIndex + 1].x - _pVelocity[linearIndex - 1].x +
               _pVelocity[linearIndex + _trueWidth].y -
               _pVelocity[linearIndex - _trueWidth].y);
      _pOldVelocity[linearIndex].y = 0;
      linearIndex++;
    }
  }
}

void Fluid::addGravity() {
  for (int i = 1; i < _width; i++) {
    addVelocity(sf::Vector2f(0, -1), i, 1);
  }
}

void Fluid::addDensitySource() {
  int halfCenterSizeX = _drawScale * _width / 30;
  int halfCenterSizeY = _drawScale * _height / 30;
  int centerWinX = _width * _drawScale / 2;
  int centerWinY = _height * _drawScale / 2;
  for (int i = -halfCenterSizeY; i < halfCenterSizeY + 1; i++) {
    for (int j = -halfCenterSizeX; j < halfCenterSizeX + 1; j++) {
      addDensity(1, j + centerWinX, i + centerWinY);
    }
  }
}

void Fluid::diffuseVelocity() {}

void Fluid::advectVelocity() {
  int xVectorIndex, yVectorIndex, rightOfXVectorIndex, belowYVectorIndex;
  float xVector, yVector, decimalBetweenXVectors, decimalBetweenYVectors,
      betweenXVectors, betweenYVectors;
  int linearIndex;

  const float dx = _dt * _width;
  const float dy = _dt * _height;

  for (int y = 1; y <= _height; y++) {
    for (int x = 1; x <= _width; x++) {
      linearIndex = index(x, y);
      xVector = x - dx * _pOldVelocity[linearIndex].x;
      yVector = y - dy * _pOldVelocity[linearIndex].y;

      xVector = clamp(xVector, 0.5f, _width + 0.5);

      xVectorIndex = (int)xVector;
      rightOfXVectorIndex = xVectorIndex + 1;

      yVector = clamp(yVector, 0.5f, _height + 0.5);

      yVectorIndex = (int)yVector;
      belowYVectorIndex = yVectorIndex + 1;

      betweenXVectors = xVector - xVectorIndex;
      decimalBetweenXVectors = 1 - betweenXVectors;
      betweenYVectors = yVector - yVectorIndex;
      decimalBetweenYVectors = 1 - betweenYVectors;

      _pVelocity[linearIndex].x =
          decimalBetweenXVectors *
              (decimalBetweenYVectors *
                   _pOldVelocity[index(xVectorIndex, yVectorIndex)].x +
               betweenYVectors *
                   _pOldVelocity[index(xVectorIndex, belowYVectorIndex)].x) +
          betweenXVectors *
              (decimalBetweenYVectors *
                   _pOldVelocity[index(rightOfXVectorIndex, yVectorIndex)].x +
               betweenYVectors *
                   _pOldVelocity[index(rightOfXVectorIndex, belowYVectorIndex)]
                       .x);
      _pVelocity[linearIndex].y =
          decimalBetweenXVectors *
              (decimalBetweenYVectors *
                   _pOldVelocity[index(xVectorIndex, yVectorIndex)].y +
               betweenYVectors *
                   _pOldVelocity[index(xVectorIndex, belowYVectorIndex)].y) +
          betweenXVectors *
              (decimalBetweenYVectors *
                   _pOldVelocity[index(rightOfXVectorIndex, yVectorIndex)].y +
               betweenYVectors *
                   _pOldVelocity[index(rightOfXVectorIndex, belowYVectorIndex)]
                       .y);
    }
  }
  //   setBoundary2d(1, _pVelocity);
  //   setBoundary2d(2, _pVelocity);
}

void Fluid::update(const float dt) {
  _dt = dt;

  addDensitySource();

  addGravity();

  addSource(_pVelocity, _pOldVelocity);

  swap(&_pVelocity, &_pOldVelocity);

  // diffuseUV(_viscocity);

  //   project();

  swap(&_pVelocity, &_pOldVelocity);

  // advect2d(_pVelocity, _pOldVelocity);

  //   project();

  addSource(_pDensity, _pOldDensity);

  swap(&_pDensity, &_pOldDensity);

  diffuseDensity();

  swap(&_pDensity, &_pOldDensity);

  advectDensity();

  fadeDensity();

  color();
}