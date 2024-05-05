#ifndef WINDOW_H
#define WINDOW_H

#include <vector>
#include <SDL2/SDL.h>

#include <algorithm>
#include <string>

void RenderWindow_2plots(SDL_Renderer* renderer, const std::vector<double>& new_sol_R_call, const std::vector<double>& sol_pde_C_call, double xScale, double yScale, int windowHeight);
void RenderWindow1plot(SDL_Renderer* renderer, const std::vector<double>& diff_call, double xScale, double yScale, int windowHeight);

#endif // WINDOW_H
