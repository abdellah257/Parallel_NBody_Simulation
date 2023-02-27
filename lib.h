
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define DEFAULT_N 1000          // Number of particles
#define DEFAULT_TIME 1000       // Number of iterations
#define G 6.67300e-11 
#define X_WIDTH 1.0e6           
#define Y_HEIGHT 1.0e6          
#define RBOUND 10      
#define DELTAT 0.01             // Time delta
#define THETA 1.0               // Theta parametre     
#define M_UNKNOWN 1.899e12   


struct Vec2 {

  double x, y;

  Vec2(double x = 0, double y = 0) : x(x), y(y) {}

  Vec2 operator+(const Vec2& other) const {
    return Vec2(x + other.x, y + other.y);
  }

  Vec2 operator-(const Vec2& other) const {
    return Vec2(x - other.x, y - other.y);
  }

  Vec2 operator*(double scale) const {
    return Vec2(x * scale, y * scale);
  }

  Vec2 operator/(double scale) const {
    return Vec2(x / scale, y / scale);
  }

  Vec2 &operator+= (const Vec2 &v){
    x += v.x;
    y += v.y;
    return *this;
  }

  Vec2 &operator=(const Vec2 & v) {
    x = v.x;
    y = v.y;
    return *this;
   }

  double norm() const {
    return sqrt(x * x + y * y);
  }

  Vec2 normalized() const {
    double n = norm();
    return Vec2(x / n, y / n);
  }

};


/* Cubic cell representing tree node in Barnes-Hut algo. */
struct Cell  {

   int index;                    // Index into arrays to identify particle's 
   int no_childs;                // Indicate whether cell is leaf or not
   double mass;                  // Mass of particle of total mass of subtree
   double x, y;                  // Location of cell(cube) in space
   double cx, cy;                // Location of center of mass of cell
   double width, height;         // Width, Height, and Depth of cell
   struct Cell* children[4];     // Pointers to child nodes

};


// #endif