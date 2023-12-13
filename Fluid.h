#ifndef FLUID_H
#define FLUID_H

#include "Globals.h"
#include <thread>

class Fluid {
public:
  Fluid();
  ~Fluid(){};

  Fluid &setSize(int xSize, int ySize, int scale);
  void reinitialize();
  void clear();
  void draw(sf::RenderWindow &window);
  void addVelocity(const sf::Vector2f velocity, const uint16_t x,
                   const uint16_t y);
  void addDensity(const double density, const uint16_t x, const uint16_t y);
  void update(const float dt);
  void saveImage();

private:
  void curl(const uint16_t x, const uint16_t y);
  void advectDensity();
  void boundDensity();
  void color();
  void diffuseDensity();
  void project();
  void addGravity();
  void addDensitySource();
  void diffuseVelocity();
  void advectVelocity();
  void boundVelocity();
  void boundDivergence();
  void pressurize();
  void boundPressure();
  template <typename T> void swap(T **from, T **to) {
    T *temp = *from;
    *from = *to;
    *to = temp;
  }
  template <typename T> void addSource(T *to, T *from) {
    for (int i = 0; i < _numCells; i++) {
      to[i] += from[i] * (float)_dt;
    }
  }
  template <typename T> void copy(T *from, T *to) {
    for (int i = 0; i < _numCells; i++) {
      to[i] = from[i];
    }
  }
  void fadeDensity();
  uint16_t index(const uint16_t X, const uint16_t Y) const;
  sf::Image _image;
  sf::Texture _tex;
  sf::Sprite _sprite;

  sf::VertexArray *_pVectors;
  sf::Color *_pColor;
  sf::Vector2f *_pVelocity, *_pOldVelocity, *_pVelocityStorage;
  double *_pDensity, *_pOldDensity, *_pDensityStorage, *_pCurl, _invTrueWidth,
      _invTrueHeight, _invWidth, _invHeight, _invNumCells, _dt, _diffusion,
      _densityFade, _viscosity;
  int _solverIterations, _width, _height, _numCells, _trueWidth, _trueHeight,
      _drawNumCells, _drawScale;
};

#endif