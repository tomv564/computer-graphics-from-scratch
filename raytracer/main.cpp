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

enum LightType {
	ambient,
	point,
	directional
};

struct Light {
	LightType type;
	float intensity;
	Vec3 position;
	Vec3 direction;
};

Light lights[3] = {
	{LightType::ambient, 0.2},
	{LightType::point, 0.6, {2, 1, 0}},
	{LightType::directional, 0.2, {0, 0, 0}, {1, 4, 4}}
};

Sphere spheres[4] = {{{0, -1, 3}, 1, {255, 0, 0}}, // red
                    {{2, 0, 4}, 1, {0, 0, 255}}, // blue
   									{{-2, 0, 4}, 1, {0, 255, 0}}, // green
                    {{0, -5001, 0}, 5000, {255, 255, 0}}};

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

Vec3 Add(Vec3 v1, Vec3 v2) {
  return {v1.x + v2.x, v1.y + v2.y, v1.z + v2.z};
}

float Length(Vec3 vec) {
	return sqrt(DotProduct(vec, vec));
}

Vec3 operator*(Vec3 vec, float factor) {
  return {vec.x * factor, vec.y * factor, vec.z * factor};
}

Color operator*(Color color, float factor) {
  return {static_cast<int>(color.r * factor), static_cast<int>(color.g * factor), static_cast<int>(color.b * factor)};
}

Vec3 operator*(float factor, Vec3 vec) {
  return {vec.x * factor, vec.y * factor, vec.z * factor};
}

Vec3 operator/(Vec3 vec, float divider) {
	return {vec.x / divider, vec.y / divider, vec.z / divider};
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

  float t1 = (-k2 + sqrt(discriminant)) / (2 * k1);
  float t2 = (-k2 - sqrt(discriminant)) / (2 * k1);

  return std::make_pair(t1, t2);
}

float ComputeLightning(Vec3 point, Vec3 normal)
{
	float intensity = 0.0;
    float normal_length = Length(normal);

	for (auto const &light : lights) {
		if (light.type == LightType::ambient)
		{
			 intensity += light.intensity;
		}
		else
		{
			Vec3 direction;

			if (light.type == LightType::directional)
				direction = light.direction;
			else
				direction = Subtract(light.position, point);

				  float angle = DotProduct(normal, direction);
				  if (angle > 0)
				  	intensity += (light.intensity*angle) / (normal_length*Length(direction));
		}
	}
	return intensity;
}

Color TraceRay(Vec3 origin, Vec3 direction, int min_t, int max_t) {
  float closest_t = INF;
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
    return {0, 0, 0};

  // surface point
  Vec3 point = Add(origin, direction * closest_t);
  // normal direction
  Vec3 normal = Subtract(point, closest_sphere->center);
  // normal vector should be based on unit of 1.
  normal = normal * (1 / Length(normal));

  float intensity = ComputeLightning(point, normal);

  return closest_sphere->color * intensity;
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
