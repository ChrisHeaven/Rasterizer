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
const int SCREEN_WIDTH = 800;
const int SCREEN_HEIGHT = 800;
SDL_Surface* screen;
int t;
float f = SCREEN_HEIGHT;
std::vector<Triangle> triangles;
vec3 cameraPos(0, 0, -3.001);
mat3 R;
float yaw = 0; // Yaw angle controlling camera rotation around y-axis
static vec3 anti_aliasing[SCREEN_WIDTH / 2][SCREEN_HEIGHT / 2];
static vec3 original_img[SCREEN_WIDTH][SCREEN_HEIGHT];
float depthBuffer[SCREEN_HEIGHT][SCREEN_WIDTH];
vec3 lightPos(0, -0.5, -0.7);
vec3 lightPower = 14.f * vec3( 1, 1, 1 );
vec3 indirectLightPowerPerArea = 0.5f * vec3( 1, 1, 1 );

struct Pixel
{
    int x;
    int y;
    float zinv;
    // vec3 illumination;
    vec3 pos3d;
    int triangle_index;
};

/* ----------------------------------------------------------------------------*/
/* FUNCTIONS                                                                   */
void Update();
void Draw();
void Interpolate(vec3 a, vec3 b, vector<vec3>& result);
void VertexShader(const vec3& v, Pixel& p, int triangle_index);
void DrawPolygon(const vector<vec3>& vertices, vec3 color, int triangle_index);
void DrawLineSDL(SDL_Surface* surface, Pixel a, Pixel b);
void ComputePolygonRows(const vector<Pixel>& vertexPixels, vector<Pixel>& leftPixels, vector<Pixel>& rightPixels, vec3 color);
void DrawPolygonRows(vector<Pixel>& leftPixels, vector<Pixel>& rightPixels, vec3 color);

int main(int argc, char* argv[])
{
    screen = InitializeSDL(SCREEN_WIDTH / 2, SCREEN_HEIGHT / 2);
    t = SDL_GetTicks(); // Set start value for timer.

    while (NoQuitMessageSDL())
    {
        Update();
        Draw();
    }

    SDL_SaveBMP(screen, "screenshot.bmp");
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
    SDL_FillRect(screen, 0, 0);
    if (SDL_MUSTLOCK(screen))
        SDL_LockSurface(screen);
    LoadTestModel(triangles);
    vec3 currentColor;

    for ( int y = 0; y < SCREEN_HEIGHT; y++ )
    {
        for ( int x = 0; x < SCREEN_WIDTH; x++ )
            depthBuffer[y][x] = 0;
    }

    for (size_t i = 0; i < triangles.size(); i++)
    {
        vector<vec3> vertices(3);
        vertices[0] = triangles[i].v0;
        vertices[1] = triangles[i].v1;
        vertices[2] = triangles[i].v2;
        currentColor = triangles[i].color;
        for (int v = 0; v < 3; ++v)
            DrawPolygon(vertices , currentColor, i);
    }
    if (SDL_MUSTLOCK(screen))
        SDL_UnlockSurface(screen);
    SDL_UpdateRect(screen, 0, 0, 0, 0);
}

void VertexShader(const vec3& v, Pixel& p, int triangle_index)
{
    float x_ = v[0] / (v[2] - cameraPos[2]) * f + SCREEN_WIDTH / 2;
    float y_ = v[1] / (v[2] - cameraPos[2]) * f + SCREEN_HEIGHT / 2;
    float z_ = 1 / (v[2] - cameraPos[2]);
    p.x = x_;
    p.y = y_;
    p.zinv = z_;
    p.pos3d = v;
    p.triangle_index = triangle_index;
}

void Interpolate(Pixel a, Pixel b, vector<Pixel>& result)
{
    int N = result.size();
    // vec2 step = vec2(b - a) / float(max(N - 1, 1));
    float step_x = float(b.x - a.x) / float(max(N - 1, 1));
    float step_y = float(b.y - a.y) / float(max(N - 1, 1));
    float step_zinv = float(b.zinv - a.zinv) / float(max(N - 1, 1));
    // float step_pos3d0 = float(b.pos3d[0] * b.zinv - a.pos3d[0] * a.zinv) / float(max(N - 1, 1));
    // float step_pos3d1 = float(b.pos3d[1] * b.zinv - a.pos3d[1] * a.zinv) / float(max(N - 1, 1));
    // float step_pos3d2 = float(1.0f / (b.pos3d[2] - cameraPos[2]) - 1.0f / (a.pos3d[2] - cameraPos[2])) / float(max(N - 1, 1));

    for (int i = 0; i < N; i++)
    {
        result[i].x = a.x + step_x * i;
        result[i].y = a.y + step_y * i;
        result[i].zinv = a.zinv + step_zinv * i;

        result[i].pos3d[2] = 1 / result[i].zinv + cameraPos[2];
        result[i].pos3d[1] = (result[i].y - SCREEN_HEIGHT / 2) / f * (result[i].pos3d[2] - cameraPos[2]);
        result[i].pos3d[0] = (result[i].x - SCREEN_WIDTH / 2) / f * (result[i].pos3d[2] - cameraPos[2]);
        // result[i].pos3d = vec3(a.pos3d[0] + step_pos3d0 * i, a.pos3d[1] + step_pos3d1 * i, a.pos3d[2] + step_pos3d2 * i);
        result[i].triangle_index = a.triangle_index;
        // result[i] = current;
        // current += step;
    }

    // for (int j = 0; j < N; j++)
    // {
    //     result[j].pos3d[0] = result[j].pos3d[0] / a.zinv;
    //     result[j].pos3d[1] = result[j].pos3d[1] / a.zinv;
    //     result[j].pos3d[2] = 1.0f / result[j].pos3d[2];
    // }
}

void DrawLineSDL(SDL_Surface * surface, Pixel a, Pixel b, vec3 color)
{
    int delta_x = glm::abs(a.x - b.x);
    int delta_y = glm::abs(a.y - b.y);
    int pixels = glm::max(delta_x, delta_y) + 1;
    vector<Pixel> line(pixels);
    vec3 light_area;
    // line = vector<Pixel> (pixels);
    Interpolate(a, b, line);
    for (int i = 0; i < pixels; i++)
    {
        // PutPixelSDL(surface, line[i].x, line[i].y, color);
        if ( line[i].zinv > depthBuffer[line[i].x][line[i].y] )
        {
            depthBuffer[line[i].x][line[i].y] = line[i].zinv;
            vec3 dis = lightPos - line[i].pos3d;
            float r = glm::length(dis);
            float result = dis[0] * triangles[line[i].triangle_index].normal[0] + dis[1] * triangles[line[i].triangle_index].normal[1] + dis[2] * triangles[line[i].triangle_index].normal[2];
            float camera_pos = 4.0 * 3.1415926 * r * r;
            if (result > 0.0)
                light_area = result / camera_pos * lightPower;
            else
                light_area = vec3(0.0, 0.0, 0.0);

            light_area = 0.5f * (indirectLightPowerPerArea + light_area);
            original_img[line[i].x][line[i].y] = color * light_area;
        }
    }
}

void DrawPolygon(const vector<vec3>& vertices, vec3 color, int triangle_index)
{
    int V = vertices.size();
    // Transform each vertex from 3D world position to 2D image position:
    vector<Pixel> vertexPixels(V);
    for (int i = 0; i < V; i++)
        VertexShader(vertices[i], vertexPixels[i], triangle_index);

    vector<Pixel> leftPixels;
    vector<Pixel> rightPixels;

    ComputePolygonRows(vertexPixels, leftPixels, rightPixels, color);
    DrawPolygonRows(leftPixels, rightPixels, color);

    int height = 0;
    for (int i = 0; i < SCREEN_HEIGHT; i = i + 2)
    {
        int width = 0;
        for (int j = 0; j < SCREEN_WIDTH; j = j + 2)
        {
            anti_aliasing[width][height] = (original_img[j][i] + original_img[j + 1][i] + original_img[j][i + 1] + original_img[j + 1][i + 1]) / vec3(4.0f, 4.0f, 4.0f);
            PutPixelSDL(screen, width, height, anti_aliasing[width][height]);
            width++;
        }
        height++;
    }
}

void ComputePolygonRows(const vector<Pixel>& vertexPixels, vector<Pixel>& leftPixels, vector<Pixel>& rightPixels, vec3 color)
{
    int min = numeric_limits<int>::max();
    int max = 0;
    vector<Pixel> edge;

    for (size_t i = 0; i < vertexPixels.size(); i++)
    {
        if (min > vertexPixels[i].y)
            min = vertexPixels[i].y;
        if (max < vertexPixels[i].y)
            max = vertexPixels[i].y;
    }
    int ROWS = max - min + 1;

    leftPixels = vector<Pixel>(ROWS);
    rightPixels = vector<Pixel>(ROWS);

    for (int i = 0; i < ROWS; i++)
    {
        leftPixels[i].x = +numeric_limits<int>::max();
        rightPixels[i].x = -numeric_limits<int>::max();
        leftPixels[i].y = min + i;
        rightPixels[i].y = min + i;
    }

    for (size_t i = 0; i < vertexPixels.size(); i++)
    {
        int j = (i + 1) % vertexPixels.size(); // The next vertex
        // DrawLineSDL(screen, vertexPixels[i], vertexPixels[j], color, edge);
        int delta_x = glm::abs(vertexPixels[i].x - vertexPixels[j].x);
        int delta_y = glm::abs(vertexPixels[i].y - vertexPixels[j].y);
        int pixels = glm::max(delta_x, delta_y) + 1;
        edge = vector<Pixel> (pixels);
        Interpolate(vertexPixels[i], vertexPixels[j], edge);

        for (int a = 0; a < ROWS; a++)
        {
            for (size_t b = 0; b < edge.size(); b++)
            {
                if (edge[b].y == min + a)
                {
                    if (edge[b].x < leftPixels[a].x)
                    {
                        leftPixels[a].x = edge[b].x;
                        leftPixels[a].zinv = edge[b].zinv;
                        leftPixels[a].pos3d = edge[b].pos3d;
                        leftPixels[a].triangle_index = edge[b].triangle_index;
                    }
                    if (edge[b].x > rightPixels[a].x)
                    {
                        rightPixels[a].x = edge[b].x;
                        rightPixels[a].zinv = edge[b].zinv;
                        rightPixels[a].pos3d = edge[b].pos3d;
                        rightPixels[a].triangle_index = edge[b].triangle_index;
                    }
                }
            }
        }
    }
}

void DrawPolygonRows(vector<Pixel>& leftPixels, vector<Pixel>& rightPixels, vec3 color)
{
    for (size_t i = 0; i < leftPixels.size(); i++)
        DrawLineSDL(screen, leftPixels[i], rightPixels[i], color);
}