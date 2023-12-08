#ifndef GAME_H
#define GAME_H

#include "Fluid.h"
#include "Globals.h"

class Game {
public:
  Game();
  ~Game();

  void create(int fluidSizeX, int fluidSizeY, int windowSizeScale);
  void run();
  void processEvents();
  void render();

private:
  sf::RenderWindow _window;
  sf::Event _gameEvent = sf::Event();
  int _winX, _winY, _scale;
  sf::Clock _clock;

  Fluid _fluid;
};

#endif