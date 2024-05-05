#include "Payoff.h"
#include "PDE.h"
#include "PDESolverOneStep.h"
#include "Tridiag.h"
#include "iostream"
#include "window.h"
#include <cmath>
#include <fstream>
#include "Discretisation.h"



int main(){

    //PDE COMPLETE

    Diffusion diffusion;
    diffusion.r=0.1;
    diffusion.sigma=0.1;

    PDEparam pdeparam;
    pdeparam.L=300.;
    pdeparam.l=0.;
    pdeparam.T=1.;
    pdeparam.M=1000;
    pdeparam.N=1000;
    pdeparam.type=PDEType::C_N;
    double K=100.;
    auto call =std::make_unique<Call>(K);
    auto put =std::make_unique<Put>(K);

    auto pde_c_call = std::make_unique<PDE_C>(std::move(call),diffusion,pdeparam);
    auto pde_c_put = std::make_unique<PDE_C>(std::move(put),diffusion,pdeparam);

    PDESolver pde_C_solver_call(std::move(pde_c_call));
    PDESolver pde_C_solver_put(std::move(pde_c_put));

    std::vector<double> sol_pde_C_call=pde_C_solver_call.solve();
    std::vector<double> sol_pde_C_put=pde_C_solver_put.solve();

    
    //PDE_REDUITE

    PDEparam pdeparam2;
    pdeparam2.L=log(300.);
    pdeparam2.l=log(0.5);
    pdeparam2.T=1.;
    pdeparam2.M=1000;
    pdeparam2.N=1000;
    pdeparam2.type=PDEType::IMP;
    
    auto call2 =std::make_unique<Call>(K);
    auto put2 =std::make_unique<Put>(K);

    auto pde_r_call = std::make_unique<PDE_R>(std::move(call2),diffusion,pdeparam2);
    auto pde_r_put = std::make_unique<PDE_R>(std::move(put2),diffusion,pdeparam2);

    PDESolver pde_R_solver_call(std::move(pde_r_call));
    PDESolver pde_R_solver_put(std::move(pde_r_put));

    std::vector<double> sol_pde_R_call=pde_R_solver_call.solve();
    std::vector<double> sol_pde_R_put=pde_R_solver_put.solve();

    //DIFFERENCE ENTRE LES DEUX 
    Discretisation d_C(pdeparam.l,pdeparam.L,pdeparam.N,STYPE::STANDARD);
    Discretisation d_R(0.5,300,1000,STYPE::LOG);
    std::vector<double> grid_C =d_C.getGrid();
    std::vector<double> grid_R =d_R.getGrid();
    std::vector<double> new_sol_R_call, new_sol_R_put;
    d_R.computeSolutionOnNewGrid(sol_pde_R_call,grid_C,new_sol_R_call);
    d_R.computeSolutionOnNewGrid(sol_pde_R_put,grid_C,new_sol_R_put);
    
    //CALCUL DE L'erreur
    std::vector<double> diff_call(new_sol_R_call.size(),0.0), diff_put(new_sol_R_put.size(),0.0);
    for (size_t i=0;i<new_sol_R_call.size();++i){
        diff_call[i]=std::abs(sol_pde_C_call[i]-new_sol_R_call[i]);
        diff_put[i]=std::abs(sol_pde_C_put[i]-new_sol_R_put[i]);
    }


    if (SDL_Init(SDL_INIT_VIDEO) < 0) {
            SDL_Log("SDL could not initialize! SDL_Error: %s\n", SDL_GetError());
            return 1;
        }



// Main loop flag
    bool quit = false;

        // Event handler
    SDL_Event e1;
    int windowWidth = 800;
    int windowHeight = 600;
    int usableWindowWidth = windowWidth - 20; // 10 pixels margin on each side
    int usableWindowHeight = windowHeight - 20; // 10 pixels margin on top and bottom
    double xScale = static_cast<double>(usableWindowWidth) / (pdeparam.M);
    double yScale = static_cast<double>(usableWindowHeight) / (210.); //yscale du CALL
    double yScale2 = static_cast<double>(usableWindowHeight) / (90.); // yscale du put 
    double yScale_1plot = static_cast<double>(usableWindowHeight) / (20.0); //yscale de diff_call
    double yScale_1plot2 = static_cast<double>(usableWindowHeight) / (4.0); //yscale de diff_put

    const int NUM_WINDOWS = 4;
    SDL_Window* windows[NUM_WINDOWS];
    SDL_Renderer* renderers[NUM_WINDOWS];
    const char* windowTitles[NUM_WINDOWS] = {"Graphique Call", "Graphique Différence Call", "Graphique Put", "Graphique Différence Put"};

    for (int i = 0; i < NUM_WINDOWS; ++i) {
        windows[i] = SDL_CreateWindow(windowTitles[i], SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, 800, 600, SDL_WINDOW_SHOWN);
        if (windows[i] == nullptr) {
            std::cerr << "Window " << i << " could not be created! SDL_Error: " << SDL_GetError() << std::endl;
            // Clean up previously created windows
            for (int j = 0; j < i; ++j) {
                SDL_DestroyWindow(windows[j]);
            }
            SDL_Quit();
            return 1;
        }

        renderers[i] = SDL_CreateRenderer(windows[i], -1, SDL_RENDERER_ACCELERATED);
        if (renderers[i] == nullptr) {
            std::cerr << "Renderer for window " << i << " could not be created! SDL_Error: " << SDL_GetError() << std::endl;
            // Clean up all windows and renderers created so far
            SDL_DestroyWindow(windows[i]);
            for (int j = 0; j <= i; ++j) {
                if (renderers[j]) {
                    SDL_DestroyRenderer(renderers[j]);
                }
                SDL_DestroyWindow(windows[j]);
            }
            SDL_Quit();
            return 1;
        }
    

    }
    while (!quit) {
        // Handle events on queue
        while (SDL_PollEvent(&e1) != 0) {
            // User requests quit
            if (e1.type == SDL_QUIT) {
                quit = true;
            }
        }

        // Render each window with its own dedicated function
        SDL_SetRenderDrawColor(renderers[0], 0, 0, 0, 255); // Black background
        SDL_RenderClear(renderers[0]);
        RenderWindow_2plots(renderers[0],new_sol_R_call ,sol_pde_C_call, xScale, yScale, windowHeight);
        SDL_RenderPresent(renderers[0]);

        SDL_SetRenderDrawColor(renderers[1], 0, 0, 0, 255); // Black background
        SDL_RenderClear(renderers[1]);
        RenderWindow1plot(renderers[1], diff_call, xScale, yScale_1plot, windowHeight);
        SDL_RenderPresent(renderers[1]);

        SDL_SetRenderDrawColor(renderers[2], 0, 0, 0, 255); // Black background
        SDL_RenderClear(renderers[2]);
        RenderWindow_2plots(renderers[2], new_sol_R_put,sol_pde_C_put, xScale, yScale2,windowHeight);
        SDL_RenderPresent(renderers[2]);

        SDL_SetRenderDrawColor(renderers[3], 0, 0, 0, 255); // Black background
        SDL_RenderClear(renderers[3]);
        RenderWindow1plot(renderers[3], diff_put, xScale, yScale_1plot2, windowHeight);
        SDL_RenderPresent(renderers[3]);

        }


    for (int i = 0; i < NUM_WINDOWS; ++i) {
        SDL_DestroyRenderer(renderers[i]);
        SDL_DestroyWindow(windows[i]);
    }

    SDL_Quit();




    return 0;
}
