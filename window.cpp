#include "window.h"


void RenderWindow_2plots(SDL_Renderer* renderer, const std::vector<double>& new_sol_R_call, const std::vector<double>& sol_pde_C_call, double xScale, double yScale,int windowHeight) {



    // Draw PDE Call Solution
    SDL_SetRenderDrawColor(renderer, 255, 0, 0, 255); // Red for Call Solution
    for (size_t i = 0; i < sol_pde_C_call.size() - 1; ++i) {
        int x1 = static_cast<int>(i * xScale) + 10;
        int y1 = windowHeight - 10 - static_cast<int>(sol_pde_C_call[i] * yScale);
        int x2 = static_cast<int>((i + 1) * xScale) + 10;
        int y2 = windowHeight - 10 - static_cast<int>(sol_pde_C_call[i + 1] * yScale);
        SDL_RenderDrawLine(renderer, x1, y1, x2, y2);
    }

    // Draw new PDE Call Solution
    SDL_SetRenderDrawColor(renderer, 0, 255, 0, 255); // Green for new Call Solution
    for (size_t i = 0; i < new_sol_R_call.size() - 1; ++i) {
        int x1 = static_cast<int>(i * xScale) + 10;
        int y1 = windowHeight - 10 - static_cast<int>(new_sol_R_call[i] * yScale);
        int x2 = static_cast<int>((i + 1) * xScale) + 10;
        int y2 = windowHeight - 10 - static_cast<int>(new_sol_R_call[i + 1] * yScale);
        SDL_RenderDrawLine(renderer, x1, y1, x2, y2);
    }


    

}



void RenderWindow1plot(SDL_Renderer* renderer, const std::vector<double>& diff_call, double xScale, double yScale,int windowHeight) {




    SDL_SetRenderDrawColor(renderer, 0, 255, 255, 255); // Cyan for error curve
    for (size_t i = 0; i < diff_call.size() - 1; ++i) {
        int x1 = 10 + static_cast<int>(i * xScale);
        int y1 = windowHeight - static_cast<int>(diff_call[i] * yScale) + 10;
        int x2 = 10 + static_cast<int>((i + 1) * xScale);
        int y2 = windowHeight - static_cast<int>(diff_call[i + 1] * yScale) + 10;
        SDL_RenderDrawLine(renderer, x1, y1, x2, y2);
    }


    

}