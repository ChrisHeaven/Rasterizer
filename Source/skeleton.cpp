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
vec3 camera_pos(0, 0, -3.001);
float yaw = 0.0f * 3.1415926 / 180; // Yaw angle controlling camera rotation around y-axis
mat3 R;
static vec3 anti_aliasing[SCREEN_WIDTH / 2][SCREEN_HEIGHT / 2];
static vec3 original_img[SCREEN_WIDTH][SCREEN_HEIGHT];
float depthBuffer[SCREEN_HEIGHT][SCREEN_WIDTH];
vec3 light_pos(0, -0.5, -0.7);
vec3 lightPower = 14.f * vec3(1, 1, 1);
vec3 indirectLightPowerPerArea = 0.5f * vec3(1, 1, 1);

struct Pixel
{
    int x;
    int y;
    float zinv;
    // vec3 illumination;
    vec3 pos3d;
    int triangle_index;
};

struct Intersection
{
    vec3 position;
    float distance;
    int triangle_index;
    vec3 colour;
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
bool closest_intersection(vec3 start, vec3 dir, const vector<Triangle>& triangles, Intersection& cloestIntersection);

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

    Uint8* keystate = SDL_GetKeyState(0);
    if (keystate[SDLK_UP])
    {
        // Move camera forward
        camera_pos.z += 0.1f;
    }
    if (keystate[SDLK_DOWN])
    {
        // Move camera backward
        camera_pos.z -= 0.1f;
    }
    if (keystate[SDLK_LEFT])
    {
        // Rotate the camera along the y axis
        yaw += 0.1f;
    }
    if (keystate[SDLK_RIGHT])
    {
        // Rotate the camera along the y axis
        yaw -= 0.1f;
    }
    if (keystate[SDLK_w])
    {
        // Move lightsource forward
        light_pos.z += 0.1f;
    }
    if (keystate[SDLK_s])
    {
        // Move lightsource backward
        light_pos.z -= 0.1f;
    }
    if (keystate[SDLK_a])
    {
        // Move lightsource left
        light_pos.x -= 0.1f;
    }
    if (keystate[SDLK_d])
    {
        // Move lightsource right
        light_pos.x += 0.1f;
    }
    if (keystate[SDLK_q])
    {
        // Move lightsource up
        light_pos.y -= 0.1f;
    }
    if (keystate[SDLK_e])
    {
        // Move lightsource down
        light_pos.y += 0.1f;
    }
}

void Draw()
{
    SDL_FillRect(screen, 0, 0);
    if (SDL_MUSTLOCK(screen))
        SDL_LockSurface(screen);
    LoadTestModel(triangles);
    vec3 currentColor;
    R = mat3(cos(yaw), 0, sin(yaw), 0, 1, 0, -sin(yaw), 0, cos(yaw));
    for (int y = 0; y < SCREEN_HEIGHT; y++)
    {
        for (int x = 0; x < SCREEN_WIDTH; x++)
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
    vec3 nv = (v - camera_pos) * R + camera_pos;

    float x_ = (nv[0] - camera_pos[0]) / (nv[2] - camera_pos[2]) * f + SCREEN_WIDTH / 2;
    float y_ = (nv[1] - camera_pos[1]) / (nv[2] - camera_pos[2]) * f + SCREEN_HEIGHT / 2;
    float z_ = 1 / (nv[2] - camera_pos[2]);
    p.x = x_;
    p.y = y_;
    p.zinv = z_;
    p.pos3d = nv;
    p.triangle_index = triangle_index;
}

void Interpolate(Pixel a, Pixel b, vector<Pixel>& result)
{
    int N = result.size();

    float step_x = float(b.x - a.x) / float(max(N - 1, 1));
    float step_y = float(b.y - a.y) / float(max(N - 1, 1));
    float step_zinv = float(b.zinv - a.zinv) / float(max(N - 1, 1));

    for (int i = 0; i < N; i++)
    {
        result[i].x = a.x + step_x * i;
        result[i].y = a.y + step_y * i;
        result[i].zinv = a.zinv + step_zinv * i;

        result[i].pos3d[2] = 1 / result[i].zinv + camera_pos[2];
        result[i].pos3d[1] = (result[i].y - SCREEN_HEIGHT / 2) / f * (result[i].pos3d[2] - camera_pos[2]) + camera_pos[1];
        result[i].pos3d[0] = (result[i].x - SCREEN_WIDTH / 2) / f * (result[i].pos3d[2] - camera_pos[2]) + camera_pos[0];
        result[i].triangle_index = a.triangle_index;
        // result[i] = current;
        // current += step;
    }
}

void DrawLineSDL(SDL_Surface * surface, Pixel a, Pixel b, vec3 color)
{
    int delta_x = glm::abs(a.x - b.x);
    int delta_y = glm::abs(a.y - b.y);
    int pixels = glm::max(delta_x, delta_y) + 1;
    vector<Pixel> line(pixels);
    vec3 light_area;
    Intersection inter, shadow_inter;
    // line = vector<Pixel> (pixels);
    Interpolate(a, b, line);
    for (int i = 0; i < pixels; i++)
    {
        // PutPixelSDL(surface, line[i].x, line[i].y, color);
        if (line[i].zinv > depthBuffer[line[i].x][line[i].y])
        {
            depthBuffer[line[i].x][line[i].y] = line[i].zinv;
            vec3 dis = light_pos - line[i].pos3d;
            float r = glm::length(dis);
            float result = dis[0] * triangles[line[i].triangle_index].normal[0] + dis[1] * triangles[line[i].triangle_index].normal[1] + dis[2] * triangles[line[i].triangle_index].normal[2];
            float camera_pos = 4.0 * 3.1415926 * r * r;
            if (result > 0.0)
                light_area = result / camera_pos * lightPower;
            else
                light_area = vec3(0.0, 0.0, 0.0);

            if (light_pos[1] < -0.20f)
            {
                if (line[i].pos3d[1] >= -0.20f)
                {
                    if (closest_intersection(line[i].pos3d, dis, triangles, inter))
                    {
                        vec3 dis_ = inter.position - line[i].pos3d;
                        if (r > glm::length(dis_) && result > 0.0 && line[i].triangle_index != inter.triangle_index)
                            light_area = vec3(0.0, 0.0, 0.0);
                    }

                    light_area = 0.5f * (indirectLightPowerPerArea + light_area);
                    original_img[line[i].x][line[i].y] = color * light_area;
                }
                else
                {
                    light_area = 0.5f * (indirectLightPowerPerArea + light_area);
                    original_img[line[i].x][line[i].y] = color * light_area;
                }
            }
            else
            {
                if (closest_intersection(line[i].pos3d, dis, triangles, inter))
                {
                    vec3 dis_ = inter.position - line[i].pos3d;
                    if (r > glm::length(dis_) && result > 0.0 && line[i].triangle_index != inter.triangle_index)
                        light_area = vec3(0.0, 0.0, 0.0);
                }

                light_area = 0.5f * (indirectLightPowerPerArea + light_area);
                original_img[line[i].x][line[i].y] = color * light_area;
            }

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

bool closest_intersection(vec3 start, vec3 dir, const vector<Triangle>& triangles, Intersection& cloestIntersection)
{
    // printf("aaa\n");
    bool flag = false;
    float min = 0.0;
    int triangle_index;
    vec3 v0, v1, v2, e1, e2, b, x;
    vec3 reflect_position, reflect_dir;
    mat3 A;

    for (size_t i = 0; i < triangles.size(); i++)
    {
        v0 = triangles[i].v0;
        v1 = triangles[i].v1;
        v2 = triangles[i].v2;

        e1 = v1 - v0;
        e2 = v2 - v0;
        b = start - v0;

        A = mat3(-dir, e1, e2);
        x = glm::inverse(A) * b;
        if (x[1] >= 0 && x[2] >= 0 && (x[1] + x[2]) <= 1 && x[0] > 0.03f)
        {
            if (!flag)
            {
                min = x[0];
                triangle_index = i;
            }
            flag = true;
            if (min > x[0])
            {
                min = x[0];
                triangle_index = i;
            }
        }
    }

    if (flag)
    {
        cloestIntersection.position = start + min * dir;
        cloestIntersection.distance = min;
        cloestIntersection.triangle_index = triangle_index;
        cloestIntersection.colour = triangles[triangle_index].color;
        return true;
    }
    else
        return false;
}