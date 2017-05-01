#include <iostream>
#include <glm/glm.hpp>
#include <SDL/SDL.h>
#include "SDLauxiliary.h"
#include "TestModel.h"

using namespace std;
using glm::vec3;
using glm::mat3;
using glm::ivec2;
using glm::vec2;

/* ----------------------------------------------------------------------------*/
/* GLOBAL VARIABLES                                                            */

const int SCREEN_WIDTH = 500;
const int SCREEN_HEIGHT = 500;
SDL_Surface* screen;
int t;
float f = SCREEN_HEIGHT;
std::vector<Triangle> triangles;
vec3 cameraPos( 0, 0, -3.001 );
mat3 R;
float yaw = 0; // Yaw angle controlling camera rotation around y-axis

/* ----------------------------------------------------------------------------*/
/* FUNCTIONS                                                                   */

void Update();
void Draw();
void Interpolate( vec3 a, vec3 b, vector<vec3>& result );
void VertexShader(const vec3& v, glm::ivec2& p);
void DrawPolygonEdges( const vector<vec3>& vertices );
void DrawLineSDL( SDL_Surface* surface, ivec2 a, ivec2 b, vec3 color );

int main( int argc, char* argv[] )
{
    screen = InitializeSDL( SCREEN_WIDTH, SCREEN_HEIGHT );
    t = SDL_GetTicks(); // Set start value for timer.

    while ( NoQuitMessageSDL() )
    {
        Update();
        Draw();
    }

    SDL_SaveBMP( screen, "screenshot.bmp" );
    return 0;
}

void Update()
{
    // Compute frame time:
    int t2 = SDL_GetTicks();
    float dt = float(t2 - t);
    t = t2;
    cout << "Render time: " << dt << " ms." << endl;
}

void Draw()
{
    SDL_FillRect( screen, 0, 0 );
    if ( SDL_MUSTLOCK(screen) )
        SDL_LockSurface(screen);
    LoadTestModel(triangles);
    for (size_t i = 0; i < triangles.size(); ++i)
    {
        vector<vec3> vertices(3);
        vertices[0] = triangles[i].v0;
        vertices[1] = triangles[i].v1;
        vertices[2] = triangles[i].v2;
        for (int v = 0; v < 3; ++v)
        {
            DrawPolygonEdges( vertices );
        }
    }
    if ( SDL_MUSTLOCK(screen) )
        SDL_UnlockSurface(screen);
    SDL_UpdateRect( screen, 0, 0, 0, 0 );
}

void VertexShader( const vec3& v, ivec2& p )
{
    float x_ = v[0] / (v[2] - cameraPos[2]) * f + SCREEN_WIDTH / 2;
    float y_ = v[1] / (v[2] - cameraPos[2]) * f + SCREEN_HEIGHT / 2;
    p.x = x_;
    p.y = y_;
}

void Interpolate( ivec2 a, ivec2 b, vector<ivec2>& result )
{
    int N = result.size();
    vec2 step = vec2(b - a) / float(max(N - 1, 1));
    vec2 current( a );
    for ( int i = 0; i < N; ++i )
    {
        result[i] = current;
        current += step;
    }
}

void DrawLineSDL( SDL_Surface* surface, ivec2 a, ivec2 b, vec3 color )
{
    ivec2 delta = glm::abs( a - b );
    int pixels = glm::max( delta.x, delta.y ) + 1;

    vector<ivec2> line( pixels );
    Interpolate( a, b, line );
    for (int i = 0; i < pixels; i++)
        PutPixelSDL( surface, line[i].x, line[i].y, color );
}

void DrawPolygonEdges( const vector<vec3>& vertices )
{
    int V = vertices.size();
    // Transform each vertex from 3D world position to 2D image position:
    vector<ivec2> projectedVertices( V );
    for ( int i = 0; i < V; ++i )
    {
        VertexShader( vertices[i], projectedVertices[i] );
    }
    // Loop over all vertices and draw the edge from it to the next vertex:
    for ( int i = 0; i < V; ++i )
    {
        int j = (i + 1) % V; // The next vertex
        vec3 color( 1, 1, 1 );
        DrawLineSDL( screen, projectedVertices[i], projectedVertices[j], color );
    }
}

// void ComputePolygonRows( const vector<ivec2>& vertexPixels, vector<ivec2>& leftPixels, vector<ivec2>& rightPixels )
// {
//     int min = numeric_limits<int>::max();
//     int max = 0;
//     for (int i = 0; i < vertexPixels.size(); i++)
//     {
//         if (min > vertexPixels[i].y)
//             min = vertexPixels[i].y;
//         if (max < vertexPixels[i].y)
//             max = vertexPixels[i].y;
//     }
//     int ROWS = max - min + 1;

//     leftPixels = vector<ivec2>(ROWS);
//     rightPixels = vector<ivec2>(ROWS);
//     leftPixels[i].x = +numeric_limits<int>::max();
//     rightPixels[i].x = -numeric_limits<int>::max();
// }