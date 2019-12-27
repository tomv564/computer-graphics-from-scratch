#include <iostream>

#include <SDL2/SDL.h>
#include <utility>

#define WINDOW_WIDTH 600
#define HALF_CANVAS 300
// #define CANVAS_WIDTH 600
#define VIEWPORT_SIZE 1
#define PROJECTION_PLANE_Z 1
#define INF 100000

struct Vec3 {
  float x, y, z;
};

struct Color {
  int r, g, b;
};

struct Sphere {
  Vec3 center;
  int radius;
  Color color;
};

Sphere spheres[3] = {{{0, -1, 3}, 1, {255, 0, 0}},
                     {{2, 0, 4}, 1, {0, 0, 255}},
                     {{-2, 0, 4}, 1, {0, 255, 0}}};

Vec3 CanvasToViewport(float x, float y) {
  return {x * VIEWPORT_SIZE / WINDOW_WIDTH, y * VIEWPORT_SIZE / WINDOW_WIDTH,
          PROJECTION_PLANE_Z};
}

float DotProduct(Vec3 v1, Vec3 v2) {
  return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

Vec3 Subtract(Vec3 v1, Vec3 v2) {
  return {v1.x - v2.x, v1.y - v2.y, v1.z - v2.z};
}

std::pair<float, float> IntersectRaySphere(Vec3 origin, Vec3 direction,
                                       Sphere sphere) {
  Vec3 oc = Subtract(origin, sphere.center);
  float k1 = DotProduct(direction, direction);
  float k2 = 2 * DotProduct(oc, direction);
  float k3 = DotProduct(oc, oc) - sphere.radius * sphere.radius;

  float discriminant = k2 * k2 - 4 * k1 * k3;
  if (discriminant < 0)
    return std::make_pair(INF, INF);

  float t1 = (-k2 + sqrt(discriminant) / (2 * k1));
  float t2 = (-k2 - sqrt(discriminant) / (2 * k1));

  return std::make_pair(t1, t2);
}

Color TraceRay(Vec3 origin, Vec3 direction, int min_t, int max_t) {
  int closest_t = INF;
  const Sphere* closest_sphere = nullptr;
  for (auto const &sphere : spheres) {
    std::pair<float, float> ts = IntersectRaySphere(origin, direction, sphere);
    if (ts.first < closest_t && min_t < ts.first && ts.first < max_t) {
      closest_t = ts.first;
      closest_sphere = &sphere;
    }
    if (ts.second < closest_t && min_t < ts.second && ts.second < max_t) {
      closest_t = ts.second;
      closest_sphere = &sphere;
    }
  }

  if (closest_sphere == nullptr)
    return {255, 255, 255};

  return closest_sphere->color;
}

int main() {
  // std::cout << "hello" << std::endl;
  SDL_Window *window;
  SDL_Renderer *renderer;
  SDL_Event event;

  Vec3 camera_position = {0, 0, 0};

  SDL_Init(SDL_INIT_VIDEO);
  SDL_CreateWindowAndRenderer(WINDOW_WIDTH, WINDOW_WIDTH, 0, &window,
                              &renderer);
  SDL_RenderSetScale(renderer, 2, 2);
  SDL_SetRenderDrawColor(renderer, 255, 255, 255, 255);
  SDL_RenderClear(renderer);

  for (int x = -HALF_CANVAS; x < HALF_CANVAS; ++x) {
    for (int y = -HALF_CANVAS; y < HALF_CANVAS; ++y) {
      auto direction = CanvasToViewport(x, y);
      auto color = TraceRay(camera_position, direction, 1, INF /* Inf ? */);
      SDL_SetRenderDrawColor(renderer, color.r, color.g, color.b, 255);

      int screen_x = x + HALF_CANVAS;
      int screen_y = HALF_CANVAS - y - 1;
      SDL_RenderDrawPoint(renderer, screen_x, screen_y);
      // std::cout << screen_x << screen_y << std::endl;
    }
  }

  SDL_RenderPresent(renderer);

  while (1) {
    if (SDL_PollEvent(&event) && event.type == SDL_QUIT)
      break;
  }

  SDL_DestroyRenderer(renderer);
  SDL_DestroyWindow(window);
  SDL_Quit();
  return 0;
}
