#pragma once
#include <algorithm>
#include <array>
#include <cassert>
#include <cstdio>
#include <cstdint>
#include <cstring>
#include <random>
#include <sys/ioctl.h>
#include <termios.h>
#include <unistd.h>

#include "geometry.cpp"

struct Color {
  uint32_t v_;

  constexpr Color() : v_(0) {}
  constexpr Color(uint32_t v) : v_(v) {}
  constexpr Color(uint8_t r, uint8_t g, uint8_t b, uint8_t a = 255U) : v_(0) {
    v_ |= r;
    v_ <<= 8;
    v_ |= g;
    v_ <<= 8;
    v_ |= b;
    v_ <<= 8;
    v_ |= a;
  }

  constexpr uint8_t r() const { return uint8_t(v_>>(8*3)); }
  constexpr uint8_t g() const { return uint8_t(v_>>(8*2)); }
  constexpr uint8_t b() const { return uint8_t(v_>>(8*1)); }
  constexpr uint8_t a() const { return uint8_t(v_>>(8*0)); }

  operator bool() const { return v_; }
};

struct Canvas {
  static std::array<unsigned short, 2> window_size() {
    struct winsize ws;
    ioctl (STDOUT_FILENO, TIOCGWINSZ, &ws);
    return {ws.ws_row, ws.ws_col};
  }

  const int H, W, AA;
  Color background;

  static constexpr int buffer_sz = 1e8;
  char char_buffer[buffer_sz];
  std::size_t char_buffer_pos;
  Color *pixel_buffer;
  double *depth_buffer;

  Canvas(int H_, int W_, int AA_)
    : H(H_), W(W_), AA(AA_), background(0, 0, 0),
      pixel_buffer(new Color[(H+1)/2*2*W]()), 
      depth_buffer(new double[(H+1)/2*2*W]()) {
    assert(H%AA == 0);
    assert(W%AA == 0);

    auto [wh, ww] = window_size();
    for (int i = 0; i < wh-1; ++i)
      char_buffer[i] = '\n';
    char_buffer[wh-1] = 0;
    clear();
    fputs("\033[?25l", stdout);
    fputs(char_buffer, stdout);
  }
  Canvas(Canvas&) = delete;

  ~Canvas() {
    fputs("\033[?25h\033[H\033[2J", stdout);
    delete pixel_buffer;
  }

  void clear() {
    std::fill(pixel_buffer, pixel_buffer+(H+1)/2*2*W, background);
    std::fill(depth_buffer, depth_buffer+(H+1)/2*2*W,
        std::numeric_limits<double>::infinity());
  }

  void set_pixel(int i, int j, double depth, Color c) {
    int I = i*W + j;
    if (depth_buffer[I] > depth) {
      depth_buffer[I] = depth;
      pixel_buffer[I] = c;
    }
  }

  void gen_frame() {
    char_buffer_pos = 0;
    append("\033[H");

    auto [wh, ww] = window_size();
    int yextra = wh-(H/AA+1)/2;
    int xextra = ww-W/AA;

    if (yextra < 0 || xextra < 0) {
      append("\033[32;1;mTOO SMALL WINDOW\033[m");
      return;
    }

    int extra_top = yextra/2;
    int extra_left = xextra/2;

    for (int i = 0; i < extra_top; ++i)
      char_buffer[char_buffer_pos++] = '\n';
    for (int i = 0; i < H/AA; i += 2) {
      for (int j = 0; j < extra_left; ++j)
        char_buffer[char_buffer_pos++] = ' ';
      for (int j = 0; j < W/AA; ++j) {

        uint32_t fr = 0;
        uint32_t fg = 0;
        uint32_t fb = 0;
        uint32_t fa = 0;

        uint32_t br = 0;
        uint32_t bg = 0;
        uint32_t bb = 0;
        uint32_t ba = 0;

        for (int di = 0; di < AA; ++di) {
          for (int dj = 0; dj < AA; ++dj) {
            Color c = pixel_buffer[(i*AA+di)*W+j*AA+dj];
            fr += (uint32_t)c.r()*c.r();
            fg += (uint32_t)c.g()*c.g();
            fb += (uint32_t)c.b()*c.b();
            fa += (uint32_t)c.a()*c.a();
          }
        }

        for (int di = 0; di < AA; ++di) {
          for (int dj = 0; dj < AA; ++dj) {
            Color c = pixel_buffer[((i+1)*AA+di)*W+j*AA+dj];
            br += (uint32_t)c.r()*c.r();
            bg += (uint32_t)c.g()*c.g();
            bb += (uint32_t)c.b()*c.b();
            ba += (uint32_t)c.a()*c.a();
          }
        }

        fr /= AA*AA;
        fg /= AA*AA;
        fb /= AA*AA;
        fa /= AA*AA;
        br /= AA*AA;
        bg /= AA*AA;
        bb /= AA*AA;
        ba /= AA*AA;

        fr = (uint32_t)std::sqrt(fr);
        fg = (uint32_t)std::sqrt(fg);
        fb = (uint32_t)std::sqrt(fb);
        fa = (uint32_t)std::sqrt(fa);
        br = (uint32_t)std::sqrt(br);
        bg = (uint32_t)std::sqrt(bg);
        bb = (uint32_t)std::sqrt(bb);
        ba = (uint32_t)std::sqrt(ba);

        Color fc{uint8_t(fr), uint8_t(fg), uint8_t(fb), uint8_t(fa)};
        Color bc{uint8_t(br), uint8_t(bg), uint8_t(bb), uint8_t(ba)};
        
        append_pixel(fc, bc);
      }
      append("\033[0m");
      if (i+2 < H)
        append("\n");
    }
  }

  template<int N>
  void append(const char (&s)[N]) {
    assert(char_buffer_pos+N < sizeof(char_buffer));
    for (int i = 0; i < N-1; ++i)
      char_buffer[char_buffer_pos++] = s[i];
  }

  void append_int(int x) {
    int w = 0;
    int tmp = x;
    do { ++w; } while (tmp /= 10);
    char_buffer_pos += w;
    assert(char_buffer_pos <= sizeof(char_buffer));
    w = 0;
    do { char_buffer[char_buffer_pos + --w] = char(x%10 + '0');
    } while (x /= 10);
  }

  void append_pixel(Color fg, Color bg) {
    if (!fg && !bg) {
      append("\033[m ");
      return;
    }

    if (!fg) {
      append("\033[m");
      append("\033[38;2;");
      append_int(bg.r());
      append(";");
      append_int(bg.g());
      append(";");
      append_int(bg.b());
      append("m▄");
      return;
    }

    if (!bg) {
      append("\033[m");
      append("\033[38;2;");
      append_int(fg.r());
      append(";");
      append_int(fg.g());
      append(";");
      append_int(fg.b());
      append("m▀");
      return;
    }

    append("\033[38;2;");
    append_int(fg.r());
    append(";");
    append_int(fg.g());
    append(";");
    append_int(fg.b());
    append("m\033[48;2;");
    append_int(bg.r());
    append(";");
    append_int(bg.g());
    append(";");
    append_int(bg.b());
    append("m▀");
  }

  void put_frame() {
    char_buffer[char_buffer_pos] = 0;
    fputs(char_buffer, stdout);
    fflush(stdout);
  }

};

struct Paint {
  static constexpr double EPS = 1e-9;

  static void triangle(Canvas &canvas,
      double x0, double y0, double z0,
      double x1, double y1, double z1,
      double x2, double y2, double z2,
      Color color) {

    double px[] = {x0, x1, x2};
    double py[] = {y0, y1, y2};
    double pz[] = {z0, z1, z2};

    double vx[3], vy[3], vz[3], vy_inv[3];

    for (int I = 0; I < 3; ++I) {
      vx[I] = px[(I+1)%3] - px[I];
      vy[I] = py[(I+1)%3] - py[I];
      vz[I] = pz[(I+1)%3] - pz[I];
      vy_inv[I] = 1/vy[I];
    }

    double area = vy[0]*vx[1] - vy[1]*vx[0];
    if (area < 0) return triangle(canvas, x0, y0, z0, x2, y2, z2, x1, y1, z1, color);

    double r = vy[0]*vz[1] - vz[0]*vy[1];
    double s = vz[0]*vx[1] - vx[0]*vz[1];
    double t = vx[0]*vy[1] - vy[0]*vx[1];

    double t_inv = 1/t;

    double min_y = std::min({y0, y1, y2});
    double max_y = std::max({y0, y1, y2});
    double all_min_x = std::min({x0, x1, x2});
    double all_max_x = std::max({x0, x1, x2});

    for (int y = std::max((int)ceil(min_y), 0); y <= std::min((int)max_y, canvas.H-1); ++y) {

      double min_x = all_min_x, max_x = all_max_x;

      for (int I = 0; I < 3; ++I) {
        if (std::abs(vy[I]) < EPS) continue;
        // vy * x >= vy * px - vx * py + vx * y
        double rhs = vy[I] * px[I] - vx[I] * py[I] + vx[I] * y;
        if (vy[I] > 0) min_x = std::max(min_x, rhs*vy_inv[I]);
        else max_x = std::min(max_x, rhs*vy_inv[I]);
      }

      for (int x = std::max((int)ceil(min_x), 0); x <= std::min((int)max_x, canvas.W-1); ++x) {
        double z = t_inv * (r*x0 + s*y0 + t*z0 - r*x - s*y);
        canvas.set_pixel(y, x, z, color);
      }
    }
  }

  template<class Container>
  static void convex(Canvas &canvas,
      const Container &vertices,
      Color color) {
    for (int i = 1; i+1 < (int)vertices.size(); ++i)
      triangle(canvas,
        vertices.begin()[0  ].x, vertices.begin()[0  ].y, vertices.begin()[0  ].z,
        vertices.begin()[i  ].x, vertices.begin()[i  ].y, vertices.begin()[i  ].z,
        vertices.begin()[i+1].x, vertices.begin()[i+1].y, vertices.begin()[i+1].z,
        color);
  }
};
