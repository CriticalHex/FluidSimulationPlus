#include "Game.h"
#include <iostream>

using namespace std;

Game::Game() { _winX = _winY = _scale = 0; }
Game::~Game() {
  // TODO
}

void Game::create(int fluidSizeX, int fluidSizeY, int windowSizeScale = 1) {
  _scale = windowSizeScale;
  _winX = fluidSizeX * _scale, _winY = fluidSizeY * _scale;
  _window.create(sf::VideoMode(_winX, _winY), "Fluid Simulation",
                 sf::Style::None);
  _window.setVerticalSyncEnabled(true);
  _window.setFramerateLimit(144);
  _fluid.setSize(fluidSizeX, fluidSizeY, _scale);
}

void Game::run() {
  while (_window.isOpen()) {
    processEvents();
    _fluid.update(_clock.getElapsedTime().asSeconds());
    _clock.restart();
    render();
  }
}

void Game::processEvents() {
  while (_window.pollEvent(_gameEvent)) {
    switch (_gameEvent.type) {
    case sf::Event::Closed:
      _window.close();
      break;

    case sf::Event::KeyPressed:
      switch (_gameEvent.key.code) {
      case sf::Keyboard::Q:
        _window.close();
        break;
      }
      break;
    case sf::Event::MouseMoved:
    case sf::Event::MouseButtonPressed:
      sf::Vector2i mousePos = sf::Mouse::getPosition(_window);
      if (sf::Mouse::isButtonPressed(sf::Mouse::Left))
        for (int i = 0; i < (_winY / 15) + 1; i++) {
          for (int j = 0; j < (_winX / 15) + 1; j++) {
            _fluid.addDensity(.1, mousePos.x + j, mousePos.y + i);
          }
        }
      break;
    }
  }
}

void Game::render() {
  _window.clear();
  _fluid.draw(_window);
  _window.display();
}