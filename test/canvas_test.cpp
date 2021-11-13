#include <unistd.h>
#include "../src/canvas.cpp"

int main() {

  const int AA = 2;
  const int H = 144*AA, W = 256*AA;
  Canvas c(H, W, AA);

  Paint::convex(c,
    {{0, 0, 0}, {W, 0, 0}, {W, H, 0}},
    {255, 0, 0});
  c.gen_frame();
  c.put_frame();

  usleep(1e6);
}
