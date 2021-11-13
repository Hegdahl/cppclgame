#pragma GCC optimize("O3,unroll-loops")
#pragma GCC target("avx2")
#include <array>
#include <chrono>
#include <random>
#include "../src/canvas.cpp"
#include "../src/geometry.cpp"

std::mt19937 rng(1337);

using quad = std::array<V3, 4>;

int main() {
  quad q_back {{
    {-1, -1, -1},
    {1, -1, -1},
    {1, 1, -1},
    {-1, 1, -1},
  }};
  quad q_front {{
    {-1, -1, 1},
    {1, -1, 1},
    {1, 1, 1},
    {-1, 1, 1},
  }};
  quad q_left {{
    {-1, -1, -1},
    {-1, 1, -1},
    {-1, 1, 1},
    {-1, -1, 1},
  }};
  quad q_right {{
    {1, -1, -1},
    {1, 1, -1},
    {1, 1, 1},
    {1, -1, 1},
  }};
  quad q_top {{
    {-1, -1, -1},
    {1, -1, -1},
    {1, -1, 1},
    {-1, -1, 1},
  }};
  quad q_bottom {{
    {-1, 1, -1},
    {1, 1, -1},
    {1, 1, 1},
    {-1, 1, 1},
  }};

  const int AA = 4;
  const int H = 144*AA;
  const int W = 240*AA;
  Canvas c(H, W, AA);
  V3 center {W/2, H/2, 0};

  double ax = 0;
  double ay = 0;
  double az = 0;

  const int TARGET_FPS = 60;

  for (int i = 0; i < TARGET_FPS*60; ++i) {
    auto t0 = std::chrono::steady_clock::now();

    auto rotation = identity<double, 3>() * (30.0 * AA);

    ax += 0.01;
    ay += 0.02;
    az += 0.03;

    rotation *= rot_x(ax);
    rotation *= rot_y(ay);
    rotation *= rot_z(az);

    Paint::convex(c, q_back * rotation + center, {255, 255, 255});
    Paint::convex(c, q_front * rotation + center, {255, 255, 0});
    Paint::convex(c, q_left * rotation + center, {255, 0, 0});
    Paint::convex(c, q_right * rotation + center, {255, 127, 0});
    Paint::convex(c, q_top * rotation + center, {0, 0, 255});
    Paint::convex(c, q_bottom * rotation + center, {0, 255, 0});

    c.gen_frame();
    c.put_frame();
    c.clear();

    auto tf = std::chrono::steady_clock::now();
    long time_left = long(1e6/TARGET_FPS) -
      std::chrono::duration_cast<std::chrono::microseconds>(tf-t0).count();
    if (time_left > 0)
      usleep((int)time_left);
  }
}
