//task list
//
//TODO: operator overloads for vector math
//TODO: better/more precise performance measurements
//TODO: build options
//TODO: put library files into external directory and symlink into projects
//
//TODO: perspective correctness
//TODO: global material database
//TODO: switch to view space for lighting calculations
//TODO: more physically accurate fresnel
//TODO: per-material reflectivity
//TODO: fresnel calculation based on IOR?
//TODO: deferred rendering?
//
//TODO: Windows build
//TODO: packaged Mac build
//
//TODO: basic player control
//TODO: fixed update tick
//TODO: fixed camera
//TODO: basic terrain pieces
//TODO: random terrain gen?
//
//TODO: SIMD rendering
//TODO: multithreaded rendering?

//NOTE(miles): I'm not sure if this is right! <SDL.h> alone doesn't compile for me on OSX,
//but I think that's what you're supposed to use, at least on other platforms.
//let me know if you have problems with this include
#ifdef __APPLE__
# include <SDL2/SDL.h>
#else
# include <SDL.h>
#endif

#include "system.hpp"
#include "std.hpp"
#include "math.hpp"
#include "ArrayList.hpp"

template <typename TYPE>
inline void swap(TYPE & a, TYPE & b) {
    TYPE temp = a;
    a = b;
    b = temp;
}

////////////////////////////////////////////////////////////////////////////////
///
////////////////////////////////////////////////////////////////////////////////

u64 applicationStartupTimeValue;

double get_time() {
    u64 currentTimeValue = SDL_GetPerformanceCounter();
    u64 diffTimeValue = currentTimeValue - applicationStartupTimeValue;
    double elapsedSeconds = (double)diffTimeValue / (double)SDL_GetPerformanceFrequency();
    return elapsedSeconds;
}

////////////////////////////////////////////////////////////////////////////////
///
////////////////////////////////////////////////////////////////////////////////

//debug
inline void print_matrix(Mat4 m) {
    printf("%8.3f%8.3f%8.3f%8.3f\n", m.m00, m.m01, m.m02, m.m03);
    printf("%8.3f%8.3f%8.3f%8.3f\n", m.m10, m.m11, m.m12, m.m13);
    printf("%8.3f%8.3f%8.3f%8.3f\n", m.m20, m.m21, m.m22, m.m23);
    printf("%8.3f%8.3f%8.3f%8.3f\n", m.m30, m.m31, m.m32, m.m33);
}

////////////////////////////////////////////////////////////////////////////////
///
////////////////////////////////////////////////////////////////////////////////

struct Vert {
    float x, y, z; //fully transformed position
    float vx, vy, vz; //position in world space (for lighting)
    float lx, ly, lz; //position in light space (for shadow mapping)
};

struct VertPair {
    u16 p, n; //position index, normal index
};

struct Triangle {
    VertPair pairs[3];
    u16 mat; //material index
};

struct Model {
    size_t vertexCount;
    Vec3 * vertices;
    Vert * verts; //post-transform coordinates

    size_t normalCount;
    Vec3 * normals;
    Vec3 * norms; //post-transform normals

    size_t triangleCount;
    Triangle * triangles;
};

////////////////////////////////////////////////////////////////////////////////
///
////////////////////////////////////////////////////////////////////////////////

Model * load_model(const char * path) {
    Model * model = (Model *) read_entire_file(path);

    //DEBUG PRINT
    // printf("%s\n", path);
    // printf("    vertexCount: %lu\n", model->vertexCount);
    // printf("    vertices: %lu\n", (size_t) model->vertices);
    // printf("    normalCount: %lu\n", model->normalCount);
    // printf("    normals: %lu\n", (size_t) model->normals);
    // printf("    triangleCount: %lu\n", model->triangleCount);
    // printf("    triangles: %lu\n", (size_t) model->triangles);

    //patch relative pointers
    size_t base = (size_t) model;
    model->vertices = (Vec3 *) ((size_t) model->vertices + base);
    model->normals = (Vec3 *) ((size_t) model->normals + base);
    model->triangles = (Triangle *) ((size_t) model->triangles + base);

    //allocate post-transform arrays
    model->verts = (Vert *) malloc(model->vertexCount * sizeof(Vert));
    model->norms = (Vec3 *) malloc(model->normalCount * sizeof(Vec3));

    return model;
}

void free_model(Model * model) {
    free(model->verts);
    free(model->norms);
    free(model);
}

////////////////////////////////////////////////////////////////////////////////
///
////////////////////////////////////////////////////////////////////////////////

typedef struct Pixel {
    u8 r, g, b, a;
} Color;

struct Material {
    float dr, dg, db; //diffuse color
};

// const int WIDTH = 256*2;
// const int HEIGHT = 192*2;
const int WIDTH = 192*2;
const int HEIGHT = 120*2;
const int SCALE = 3;

const int SHADOW_WIDTH = 256*2;
const int SHADOW_HEIGHT = 256*2;

//DEBUG GLOBALS
int totalTris;
int drawnTris;
int shadowTotalTris;
int shadowDrawnTris;
int globalFrame;

//TODO: reduce the massive amount of code duplication going on here!
void draw_static_model_shadow(Model * model, Mat4 modelMatrix, Mat4 shadowMatrix,
    float * shadowBuffer) {

    Mat4 modelshadow = shadowMatrix * modelMatrix;

    //transform coordinates to light space
    for (int i = 0; i < model->vertexCount; ++i) {
        Vec3 a = model->vertices[i];

        //matrix transforms ("vertex shader")
        Vec4 b = modelshadow * vec4(a, 1);

        //transform from clip space to pixel space
        //TODO: incorporate this into the light-space matrix
        //      or should we calculate this during the rasterization phase?
        b.x = (b.x + 1) *  0.5f * SHADOW_WIDTH;
        b.y = (b.y - 1) * -0.5f * SHADOW_HEIGHT;

        //OPTIMIZE: should we worry about storing this less sparsely?
        model->verts[i] = { b.x, b.y, b.z, 0, 0, 0, 0, 0, 0 };
    }

    //draw triangles
    for (int t = 0; t < model->triangleCount; ++t) {
        Triangle tri = model->triangles[t];
        Vert triangle[3] = {
            model->verts[tri.pairs[0].p],
            model->verts[tri.pairs[1].p],
            model->verts[tri.pairs[2].p],
        };

        ++shadowTotalTris;

        //alternative method: https://www.geeksforgeeks.org/orientation-3-ordered-points/
        if ((triangle[1].y - triangle[0].y) * (triangle[2].x - triangle[1].x) -
            (triangle[1].x - triangle[0].x) * (triangle[2].y - triangle[1].y) >= 0) {
            continue;
        }

        if (triangle[0].x < 0 && triangle[1].x < 0 && triangle[2].x < 0) {
            continue;
        }
        if (triangle[0].x > SHADOW_WIDTH && triangle[1].x > SHADOW_WIDTH &&
            triangle[2].x > SHADOW_WIDTH) {
            continue;
        }

        if (triangle[0].y < 0 && triangle[1].y < 0 && triangle[2].y < 0) {
            continue;
        }
        if (triangle[0].y > SHADOW_HEIGHT && triangle[1].y > SHADOW_HEIGHT &&
            triangle[2].y > SHADOW_HEIGHT) {
            continue;
        }

        //we have no concept of near/far plane for shadow buffer,
        //so we don't discard based on z

        ++shadowDrawnTris;

        struct Edge {
            Vert v1, v2; //sorted by y
            float yfactor; //factor for interpolating vertex attributes vertically
        };

        Edge edges[3];

        float miny = SHADOW_HEIGHT;
        float maxy = 0;

        //convert triangles to edges
        for (int i = 0; i < 3; ++i) {
            Vert v1 = triangle[i];
            Vert v2 = triangle[(i + 1) % 3];

            //update the triangle's vertical extent
            miny = v1.y < miny? v1.y : miny;
            maxy = v1.y > maxy? v1.y : maxy;

            //sort vertices by y
            edges[i].v1 = v1.y < v2.y? v1 : v2;
            edges[i].v2 = v1.y > v2.y? v1 : v2;
            edges[i].yfactor = 1.0f/(edges[i].v2.y - edges[i].v1.y);
        }

        //convert the triangle's vertical extent to pixels
        int firstLine = miny + 1;
        int lastLine = maxy;
        //clamp vertical extent of triangle to within the screen for rasterization
        if (firstLine < 0) firstLine = 0;
        if (lastLine > SHADOW_HEIGHT - 1) lastLine = SHADOW_HEIGHT - 1;

        for (int y = firstLine; y <= lastLine; ++y) {
            float * srow = shadowBuffer + y * SHADOW_WIDTH;

            //the current pixel row will be within the vertical extend of only two
            //of the three edges at any time, so find those two and discard the third
            Edge e1, e2;
            if (y < edges[0].v1.y || y > edges[0].v2.y) {
                e1 = edges[1];
                e2 = edges[2];
            } else if (y < edges[1].v1.y || y > edges[1].v2.y) {
                e1 = edges[0];
                e2 = edges[2];
            } else {
                e1 = edges[0];
                e2 = edges[1];
            }

            //calculate vertical blend amounts for this scanline
            float f1a = (e1.v2.y - y) * e1.yfactor;
            float f2a = (e2.v2.y - y) * e2.yfactor;
            float f1b = 1 - f1a;
            float f2b = 1 - f2a;

            //find intersection with each edge by interpolating x along the edge
            float x1 = f1a * e1.v1.x + f1b * e1.v2.x;
            float x2 = f2a * e2.v1.x + f2b * e2.v2.x;

            //sort edges based on intersections
            float minx = x1 < x2? x1 : x2;
            float maxx = x1 > x2? x1 : x2;
            if (x1 > x2) {
                swap(e1, e2);
                swap(f1a, f2a);
                swap(f1b, f2b);
            }

            //interpolate vertex attributes at intersection points
            float z1 = f1a * e1.v1.z + f1b * e1.v2.z;
            float z2 = f2a * e2.v1.z + f2b * e2.v2.z;
            //factor for interpolating vertex attributes horizontally
            float xfactor = 1.0f/(maxx - minx);

            //convert horizontal extent to pixels
            int first = minx + 1;
            int last = maxx;
            //clamp horizontal extent of scanline to pixels
            if (first < 0) first = 0;
            if (last > SHADOW_WIDTH - 1) last = SHADOW_WIDTH - 1;

            for (int x = first; x <= last; ++x) {
                //calculate horizontal interpolation factor for this pixel
                float fa = (maxx - x) * xfactor;
                float fb = 1 - fa;

                //interpolate vertex attributes for this pixel
                float z = fa * z1 + fb * z2;

                srow[x] = fmin(srow[x], z);
            }
        }
    }
}

void draw_static_model(Model * model,
    Mat4 modelMatrix, Mat4 viewMatrix, Mat4 projectionMatrix, Mat4 shadowMatrix,
    SDL_Surface * canvas, float * depthBuffer, float * shadowBuffer,
    Vec3 lightPos, Vec3 lightDir, Vec3 cameraPos, Material * materials) {

    Mat4 viewproj = projectionMatrix * viewMatrix;
    Mat4 modelviewproj = viewproj * modelMatrix;
    Mat4 modelshadow = shadowMatrix * modelMatrix;
    Mat3 normalMatrix = inverse_transpose(mat3(modelMatrix));

    //transform positions
    for (int i = 0; i < model->vertexCount; ++i) {
        Vec4 a = vec4(model->vertices[i], 1);

        //matrix transforms ("vertex shader")
        Vec4 b = modelviewproj * a;
        Vec4 v = modelMatrix * a;
        Vec4 l = modelshadow * a;

        //I don't know why we need to knock out the sign of w.
        //I just know that if we don't, horrble things happen.
        b.w = fabs(b.w);

        //perspective divide
        b.x /= b.w;
        b.y /= b.w;
        b.z /= b.w;

        //transform from clip space to screen space
        //TODO: incorporate this into the modelviewproj matrix
        //      (or should we calculate this during the rasterization phase?)
        b.x = (b.x + 1) *  0.5f * WIDTH;
        b.y = (b.y - 1) * -0.5f * HEIGHT;

        l.x = (l.x + 1) *  0.5f * SHADOW_WIDTH;
        l.y = (l.y - 1) * -0.5f * SHADOW_HEIGHT;

        model->verts[i] = { b.x, b.y, b.z, v.x, v.y, v.z, l.x, l.y, l.z };
    }

    //transform normals
    for (int i = 0; i < model->normalCount; ++i) {
        model->norms[i] = normalize(normalMatrix * model->normals[i]);
    }

    //draw triangles
    for (int t = 0; t < model->triangleCount; ++t) {
        Triangle tri = model->triangles[t];
        Vert triangle[3] = {
            model->verts[tri.pairs[0].p],
            model->verts[tri.pairs[1].p],
            model->verts[tri.pairs[2].p],
        };
        Vec3 normals[3] = {
            model->norms[tri.pairs[0].n],
            model->norms[tri.pairs[1].n],
            model->norms[tri.pairs[2].n],
        };

        ++totalTris;

        //alternative method: https://www.geeksforgeeks.org/orientation-3-ordered-points/
        if ((triangle[1].y - triangle[0].y) * (triangle[2].x - triangle[1].x) -
            (triangle[1].x - triangle[0].x) * (triangle[2].y - triangle[1].y) <= 0) {
            continue;
        }

        if (triangle[0].x < 0 && triangle[1].x < 0 && triangle[2].x < 0) {
            continue;
        }
        if (triangle[0].x > WIDTH && triangle[1].x > WIDTH && triangle[2].x > WIDTH) {
            continue;
        }

        if (triangle[0].y < 0 && triangle[1].y < 0 && triangle[2].y < 0) {
            continue;
        }
        if (triangle[0].y > HEIGHT && triangle[1].y > HEIGHT && triangle[2].y > HEIGHT) {
            continue;
        }

        if (triangle[0].z < -1 && triangle[1].z < -1 && triangle[2].z < -1) {
            continue;
        }
        if (triangle[0].z > 1 && triangle[1].z > 1 && triangle[2].z > 1) {
            continue;
        }

        ++drawnTris;

        struct Edge {
            Vert v1, v2; //sorted by y
            Vec3 n1, n2; //sorted by y
            float yfactor; //factor for interpolating vertex attributes vertically
        };

        //URGENT TODO: keep normals attached to vertices!

        Edge edges[3];

        float miny = HEIGHT;
        float maxy = 0;

        //convert triangles to edges
        for (int i = 0; i < 3; ++i) {
            Vert v1 = triangle[i];
            Vert v2 = triangle[(i + 1) % 3];

            Vec3 n1 = normals[i];
            Vec3 n2 = normals[(i + 1) % 3];

            //update the triangle's vertical extent
            miny = v1.y < miny? v1.y : miny;
            maxy = v1.y > maxy? v1.y : maxy;

            //sort vertices by y
            edges[i].v1 = v1.y < v2.y? v1 : v2;
            edges[i].v2 = v1.y < v2.y? v2 : v1;
            //and the normals, too
            edges[i].n1 = v1.y < v2.y? n1 : n2;
            edges[i].n2 = v1.y < v2.y? n2 : n1;

            edges[i].yfactor = 1.0f/(edges[i].v2.y - edges[i].v1.y);
        }

        //convert the triangle's vertical extent to pixels
        int firstLine = miny + 1;
        int lastLine = maxy;
        //clamp vertical extent of triangle to within the screen for rasterization
        if (firstLine < 0) firstLine = 0;
        if (lastLine > HEIGHT - 1) lastLine = HEIGHT - 1;

        for (int y = firstLine; y <= lastLine; ++y) {
            Pixel * row = (Pixel *) ((u8 *)canvas->pixels + y * canvas->pitch);
            float * zrow = depthBuffer + y * WIDTH;

            //the current pixel row will be within the vertical extend of only two
            //of the three edges at any time, so find those two and discard the third
            Edge e1, e2;
            if (y < edges[0].v1.y || y > edges[0].v2.y) {
                e1 = edges[1];
                e2 = edges[2];
            } else if (y < edges[1].v1.y || y > edges[1].v2.y) {
                e1 = edges[0];
                e2 = edges[2];
            } else {
                e1 = edges[0];
                e2 = edges[1];
            }

            //calculate vertical blend amounts for this scanline
            float f1a = (e1.v2.y - y) * e1.yfactor;
            float f2a = (e2.v2.y - y) * e2.yfactor;
            float f1b = 1 - f1a;
            float f2b = 1 - f2a;

            //find intersection with each edge by interpolating x along the edge
            float x1 = f1a * e1.v1.x + f1b * e1.v2.x;
            float x2 = f2a * e2.v1.x + f2b * e2.v2.x;

            //sort edges based on intersections
            float minx = x1 < x2? x1 : x2;
            float maxx = x1 > x2? x1 : x2;
            if (x1 > x2) {
                swap(e1, e2);
                swap(f1a, f2a);
                swap(f1b, f2b);
            }

            //interpolate vertex attributes at intersection points
            float z1 = f1a * e1.v1.z + f1b * e1.v2.z;
            float z2 = f2a * e2.v1.z + f2b * e2.v2.z;

            float nx1 = f1a * e1.n1.x + f1b * e1.n2.x;
            float nx2 = f2a * e2.n1.x + f2b * e2.n2.x;
            float ny1 = f1a * e1.n1.y + f1b * e1.n2.y;
            float ny2 = f2a * e2.n1.y + f2b * e2.n2.y;
            float nz1 = f1a * e1.n1.z + f1b * e1.n2.z;
            float nz2 = f2a * e2.n1.z + f2b * e2.n2.z;

            float vx1 = f1a * e1.v1.vx + f1b * e1.v2.vx;
            float vx2 = f2a * e2.v1.vx + f2b * e2.v2.vx;
            float vy1 = f1a * e1.v1.vy + f1b * e1.v2.vy;
            float vy2 = f2a * e2.v1.vy + f2b * e2.v2.vy;
            float vz1 = f1a * e1.v1.vz + f1b * e1.v2.vz;
            float vz2 = f2a * e2.v1.vz + f2b * e2.v2.vz;

            float lx1 = f1a * e1.v1.lx + f1b * e1.v2.lx;
            float lx2 = f2a * e2.v1.lx + f2b * e2.v2.lx;
            float ly1 = f1a * e1.v1.ly + f1b * e1.v2.ly;
            float ly2 = f2a * e2.v1.ly + f2b * e2.v2.ly;
            float lz1 = f1a * e1.v1.lz + f1b * e1.v2.lz;
            float lz2 = f2a * e2.v1.lz + f2b * e2.v2.lz;
            //factor for interpolating vertex attributes horizontally
            float xfactor = 1.0f/(maxx - minx);

            //convert horizontal extent to pixels
            int first = minx + 1;
            int last = maxx;
            //clamp horizontal extent of scanline to pixels
            if (first < 0) first = 0;
            if (last > WIDTH - 1) last = WIDTH - 1;

            //fetch material data for this triangle
            Material material = materials[tri.mat];
            float dr = material.dr;
            float dg = material.dg;
            float db = material.db;

            for (int x = first; x <= last; ++x) {
                //calculate horizontal interpolation factor for this pixel
                float fa = (maxx - x) * xfactor;
                float fb = 1 - fa;

                //interpolate vertex attributes for this pixel
                float z = fa * z1 + fb * z2;

                //depth test early-out
                if (z > 1.0f || z < -1.0f || z >= zrow[x]) {
                    continue;
                }

                float nx = fa * nx1 + fb * nx2;
                float ny = fa * ny1 + fb * ny2;
                float nz = fa * nz1 + fb * nz2;

                float vx = fa * vx1 + fb * vx2;
                float vy = fa * vy1 + fb * vy2;
                float vz = fa * vz1 + fb * vz2;

                float lx = fa * lx1 + fb * lx2;
                float ly = fa * ly1 + fb * ly2;
                float lz = fa * lz1 + fb * lz2;



                //"fragment shader"
                float r, g, b;

                //shadow mapping
                int ix = lx + 0.5f;
                int iy = ly + 0.5f;
                float shadow = 1;
                if (ix >= 0 && ix < SHADOW_WIDTH && iy > 0 && iy < SHADOW_HEIGHT) {
                    int idx = iy * SHADOW_WIDTH + ix;
                    if (lz > shadowBuffer[idx] + 0.2f) {
                        shadow = 0;
                    }
                }

                //apply illumination from directional light
                Vec3 normal = normalize(vec3(nx, ny, nz));
                float diffuse = shadow;
                diffuse *= fmax(dot(normal, lightDir), 0);

                //add illumination from point light
                Vec3 pointDir = lightPos - vec3(vx, vy, vz);
                float invFalloff = len2(pointDir);
                float falloff = 1 / invFalloff;
                pointDir = normalize(pointDir);
                diffuse += 2 * fmax(dot(normal, pointDir), 0) * falloff;

                //add ambient light
                diffuse += 0.01f;

                //apply material diffuse color
                r = g = b = diffuse;
                r *= dr;
                g *= dg;
                b *= db;

                //compute specular lighting value for point light
                Vec3 viewVector = normalize(vec3(vx, vy, vz) - cameraPos);
                // float r = d - 2 * dot(d, n) * n;
                Vec3 reflected = viewVector - 2 * dot(viewVector, normal) * normal;
                float pointSpec = fmax(dot(reflected, pointDir), 0);
                pointSpec = 2.0f * pow(pointSpec, 32) / sqrtf(invFalloff);

                //compute specular lighting value for directional light
                float sunSpec = fmax(dot(reflected, lightDir), 0);
                sunSpec = 2.0f * shadow * pow(sunSpec, 32);

                //apply a loose approximation of the fresnel effect
                float fresnel = sq(sq(1 - dot(-viewVector, normal)));
                float minusFresnel = fmax(0, 1 - fresnel);
                float specular = (pointSpec + sunSpec) * fresnel;
                // float specular = sunSpec * fresnel;

                //add specular lighting
                // r += specular;
                // g += specular;
                // b += specular;

                //use specular lighting in a (hopefully) physically correct way
                r = r * minusFresnel + specular;
                g = g * minusFresnel + specular;
                b = b * minusFresnel + specular;

                // r = fresnel;
                // g = fresnel;
                // b = fresnel;

                //apply gamma correction
                r = sqrtf(r);
                g = sqrtf(g);
                b = sqrtf(b);

                //clamp
                r = fmin(r, 1);
                g = fmin(g, 1);
                b = fmin(b, 1);

                row[x] = { (u8)(r * 255), (u8)(g * 255), (u8)(b * 255), row[x].a };
                zrow[x] = z;
            }
        }
    }
}

////////////////////////////////////////////////////////////////////////////////
///
////////////////////////////////////////////////////////////////////////////////

void debug_print_model_info(Model * model) {
    Triangle * faces = model->triangles;
    printf("\n\n\nfaces:\n");
    for (int i = 0; i < model->triangleCount; ++i) {
        printf("    { %d, %d }, { %d, %d }, { %d, %d }, %d\n",
            faces[i].pairs[0].p, faces[i].pairs[0].n,
            faces[i].pairs[1].p, faces[i].pairs[1].n,
            faces[i].pairs[2].p, faces[i].pairs[2].n,
            faces[i].mat);
    }
}

//PRECONDITION: x1 >= 0, y1 >= 0, x1 < WIDTH, y1 < HEIGHT
void draw_bresenham(SDL_Surface * canvas, int x1, int y1, int x2, int y2, Color color) {
    //avoid drawing invalid lines
    if (x1 < 0 || x1 >= WIDTH  || y2 < 0 || x2 >= WIDTH ||
        y1 < 0 || y1 >= HEIGHT || y2 < 0 || y2 >= HEIGHT) {
        // printf("BRESENHAM: %f, %f | %f, %f\n");
        return;
    }

    int dx =  abs(x2 - x1);
    int dy = -abs(y2 - y1);

    int sx = x1 < x2? 1 : -1;
    int sy = y1 < y2? 1 : -1;

    int err = dx + dy;

    while (true) {
        if (x1 < 0 || x1 >= WIDTH) {
            printf("BAD x1: %d, %d\n", x1, y1);
            exit(1);
        } else if (y1 < 0 || y1 >= HEIGHT) {
            printf("BAD y1: %d, %d\n", x1, y1);
        }

        Pixel * row = (Pixel *) ((u8 *)canvas->pixels + y1 * canvas->pitch);
        row[x1] = color;

        if (x1 == x2 && y1 == y2) {
            return;
        }

        int e2 = 2 * err;
        if (e2 > dy) {
            err += dy;
            x1 += sx;
        }
        if (e2 < dx) {
            err += dx;
            y1 += sy;
        }
    }
}

void draw_line(SDL_Surface * canvas, float x1, float y1, float x2, float y2, Color color) {
    //convert to slope-intercept form
    float dx = x2 - x1;
    float dy = y2 - y1;
    float m1 = dy/dx; //slope
    float m2 = dx/dy; //inverse (reciprocal) slope
    float b1 = y1 - m1 * x1; //intercept
    float b2 = x1 - m2 * y1; //inverse (reciprocal) intercept

    //early outs
    if (x1 < 0 && x2 < 0) { return; }
    if (y1 < 0 && y2 < 0) { return; }
    if (x1 >= WIDTH && x2 >= WIDTH) { return; }
    if (y1 >= HEIGHT && y2 >= HEIGHT) { return; }

    //crop y
    if (y1 < 0) {
        y1 = 0;
        x1 = m2 * y1 + b2;
    } else if (y1 > HEIGHT - 1) {
        y1 = HEIGHT - 1;
        x1 = m2 * y1 + b2;
    }

    if (y2 < 0) {
        y2 = 0;
        x2 = m2 * y2 + b2;
    } else if (y2 > HEIGHT - 1) {
        y2 = HEIGHT - 1;
        x2 = m2 * y2 + b2;
    }

    //crop x
    if (x1 < 0) {
        x1 = 0;
        y1 = m1 * x1 + b1;
    } else if (x1 > WIDTH - 1) {
        x1 = WIDTH - 1;
        y1 = m1 * x1 + b1;
    }

    if (x2 < 0) {
        x2 = 0;
        y2 = m1 * x2 + b1;
    } else if (x2 > WIDTH - 1) {
        x2 = WIDTH - 1;
        y2 = m1 * x2 + b1;
    }

    draw_bresenham(canvas, x1, y1, x2, y2, color);
}

void draw_line(SDL_Surface * canvas, Mat4 viewproj,
               float x1, float y1, float z1, float x2, float y2, float z2, Color color) {
    Vec4 start = viewproj * vec4(x1, y1, z1, 1);
    Vec4 end   = viewproj * vec4(x2, y2, z2, 1);
    start.x /= start.w;
    start.y /= start.w;
    end.x /= end.w;
    end.y /= end.w;
    start.x = (start.x + 1) *  0.5f * WIDTH;
    start.y = (start.y - 1) * -0.5f * HEIGHT;
    end.x = (end.x + 1) *  0.5f * WIDTH;
    end.y = (end.y - 1) * -0.5f * HEIGHT;
    draw_line(canvas, start.x, start.y, end.x, end.y, color);
}

void draw_point(SDL_Surface * canvas, float x0, float y0, int radius, Color color) {
    int x = x0;
    int y = y0;

    int minx = imax(x - radius, 0);
    int miny = imax(y - radius, 0);
    int maxx = imin(x + radius, WIDTH - 1);
    int maxy = imin(y + radius, HEIGHT - 1);

    for (int y = miny; y <= maxy; ++y) {
        Pixel * row = (Pixel *) ((u8 *)canvas->pixels + y * canvas->pitch);
        for (int x = minx; x <= maxx; ++x) {
            row[x] = color;
        }
    }
}

void draw_point(SDL_Surface * canvas, Mat4 viewproj,
                float x0, float y0, float z0, int radius, Color color) {
    Vec4 point = viewproj * vec4(x0, y0, z0, 1);
    point.x /= point.w;
    point.y /= point.w;
    point.x = (point.x + 1) *  0.5f * WIDTH;
    point.y = (point.y - 1) * -0.5f * HEIGHT;

    draw_point(canvas, point.x, point.y, radius, color);
}

struct DrawCall {
    Model * model;
    Mat4 matrix; //model-to-world transform
    bool shadow; //we will probably replace this with a flags mask later
};

int main() {
    //initialize timer
    applicationStartupTimeValue = SDL_GetPerformanceCounter();

    if (SDL_Init(SDL_INIT_VIDEO | SDL_INIT_JOYSTICK | SDL_INIT_GAMECONTROLLER)) {
        printf("SDL FAILED TO INIT: %s\n", SDL_GetError());
        return 1;
    }

    printf("SDL init: %f seconds\n", get_time());

    SDL_Window * window = SDL_CreateWindow("Test Window",
        SDL_WINDOWPOS_CENTERED_DISPLAY(1), SDL_WINDOWPOS_CENTERED_DISPLAY(1),
        WIDTH * SCALE, HEIGHT * SCALE,
        SDL_WINDOW_SHOWN);
    if (window == nullptr) {
        printf("SDL FAILED TO CREATE WINDOW: %s\n", SDL_GetError());
        return 1;
    }

    printf("SDL create window: %f seconds\n", get_time());

    SDL_Surface * canvas = SDL_CreateRGBSurface(0, WIDTH, HEIGHT, 32,
        0x000000FF,
        0x0000FF00,
        0x00FF0000,
        0xFF000000);
    if (canvas == nullptr) {
        printf("SDL FAILED TO CREATE SURFACE: %s\n", SDL_GetError());
        return 1;
    }

    //DEBUG: surface for shadow map visualization
    SDL_Surface * debugShadow = SDL_CreateRGBSurface(0, SHADOW_WIDTH, SHADOW_HEIGHT, 32,
        0x000000FF,
        0x0000FF00,
        0x00FF0000,
        0xFF000000);
    if (debugShadow == nullptr) {
        printf("SDL FAILED TO CREATE DEBUG SURFACE: %s\n", SDL_GetError());
        return 1;
    }



    struct Gamepad {
        SDL_Joystick * joy;
        SDL_GameController * con;
        SDL_JoystickID id;
        //maps to enumeration SDL_CONTROLLER_AXIS_*
        float axis[SDL_CONTROLLER_AXIS_MAX]; //6
        //maps to enumeration SDL_CONTROLLER_BUTTON_*
        bool button[SDL_CONTROLLER_BUTTON_MAX]; //15
    };

    printf("sizeof SDL_Event: %lu\n", sizeof(SDL_Event));
    printf("sizeof Gamepad: %lu\n", sizeof(Gamepad));
    printf("sizeof DrawCall: %lu\n", sizeof(DrawCall));
    printf("sizeof Model: %lu\n", sizeof(Model));
    printf("SDL_PRESSED: %d\n", SDL_PRESSED);
    printf("SDL_RELEASED: %d\n", SDL_RELEASED);

    printf("SDL full init: %f seconds\n", get_time());



    Model * carModel = load_model("res/ae86.73");
    Model * tableModel = load_model("res/table.73");
    Model * droidModel = load_model("res/table.73");
    Model * floorModel = load_model("res/test_ground.73");

    //allocate buffers
    float * depthBuffer = (float *) malloc(WIDTH * HEIGHT * sizeof(float));
    float * shadowBuffer = (float *) malloc(SHADOW_WIDTH * SHADOW_HEIGHT * sizeof(float));



    //timestep and framerate info
    float frameTimes[10] = {};
    float time = get_time();
    float lastTime = 0;
    int frame = 0;



    //DEBUG: randomly scattered tables
    Vec4 * tablePos = (Vec4 *) malloc(100 * sizeof(Vec4));
    for (int i = 0; i < 100; ++i) {
        float variance = 20; //variance
        tablePos[i] = {
            random_float(-variance, variance),
            -1,
            random_float(-variance, variance),
            random_float(TWO_PI), //rotation around y axis
        };
    }

    //DEBUG: randomly assigned material colors
    Material * materials = (Material *) malloc(100 * sizeof(Material));
    for (int i = 0; i < 100; ++i) {
        materials[i].dr = sq(random_float());
        materials[i].dg = sq(random_float());
        materials[i].db = sq(random_float());
    }



    //initialize render list
    ArrayList<DrawCall> drawList = create_array_list<DrawCall>(100);



    //keymap
    //TODO: make maps for tracking key events which are reset every frame
    //      (so that we can check for keydown and keyup events outside the event loop)
    bool isDown[256] = {};

    //input mode state
    bool captured = false;
    if (captured) {
        SDL_SetRelativeMouseMode(SDL_TRUE);
    }

    const int MAX_PADS = 4;
    int padCount = 0;
    Gamepad gamepads[MAX_PADS] = {};

    //TODO: make this a Vec3 the whole way?
    Vec3 camPos = vec3(0, 0, 10);
    float cameraYaw = 0;
    float cameraPitch = 0;

    Vec3 carPos = vec3(0, 0, 0);
    float carYaw = 0;

    printf("done initializing: %f seconds\n", get_time());







    bool shouldExit = false;
    while (!shouldExit) {
        float dmx = 0;
        float dmy = 0;

        SDL_Event event;
        while (SDL_PollEvent(&event)) {
            if (event.type == SDL_QUIT) {
                shouldExit = true;
            } else if (event.type == SDL_KEYDOWN) {
                isDown[event.key.keysym.scancode] = true;
            } else if (event.type == SDL_KEYUP) {
                isDown[event.key.keysym.scancode] = false;
            } else if (event.type == SDL_MOUSEMOTION) {
                dmx += event.motion.xrel;
                dmy += event.motion.yrel;
            } else if (event.type == SDL_CONTROLLERDEVICEADDED) {
                printf("CONTROLLER ADDED: %d, %d\n",
                    event.cdevice.timestamp, event.cdevice.which);

                if (padCount < MAX_PADS) {
                    gamepads[padCount].joy = SDL_JoystickOpen(event.cdevice.which);
                    gamepads[padCount].con = SDL_GameControllerOpen(event.cdevice.which);
                    gamepads[padCount].id = SDL_JoystickInstanceID(gamepads[padCount].joy);
                    padCount += 1;
                }
            } else if (event.type == SDL_CONTROLLERDEVICEREMOVED) {
                printf("CONTROLLER REMOVED: %d, %d\n",
                    event.cdevice.timestamp, event.cdevice.which);

                for (int i = 0; i < MAX_PADS; ++i) {
                    if (gamepads[i].id == event.cdevice.which) {
                        printf("Found controller to disconnect (idx %d)\n", i);
                        SDL_GameControllerClose(gamepads[i].con);
                        SDL_JoystickClose(gamepads[i].joy);

                        //remove controller from list
                        for (int j = i + 1; j < MAX_PADS; ++j) {
                            gamepads[j - 1] = gamepads[j];
                        }
                        padCount -= 1;

                        break;
                    }
                }
            } else if (event.type == SDL_CONTROLLERAXISMOTION) {
                DEBUG_ASSERT(event.caxis.axis < SDL_CONTROLLER_AXIS_MAX);
                for (int i = 0; i < padCount; ++i) {
                    if (gamepads[i].id == event.caxis.which) {
                        gamepads[i].axis[event.caxis.axis] =
                            fmax(event.caxis.value / 32767.0f, -1);
                        break;
                    }
                }
            } else if (event.type == SDL_CONTROLLERBUTTONDOWN ||
                       event.type == SDL_CONTROLLERBUTTONUP) {
                DEBUG_ASSERT(event.cbutton.button < SDL_CONTROLLER_BUTTON_MAX);
                for (int i = 0; i < padCount; ++i) {
                    if (gamepads[i].id == event.cbutton.which) {
                        gamepads[i].button[event.cbutton.button] = event.cbutton.state;
                        break;
                    }
                }
            }
        }

        //update timestep
        lastTime = time;
        time = get_time();
        float dt = time - lastTime;



        if (isDown[SDL_SCANCODE_ESCAPE]) {
            captured = false;
            SDL_SetRelativeMouseMode(SDL_FALSE);
        } else if (isDown[SDL_SCANCODE_GRAVE]) {
            if (!captured) {
                //set the camera closer to the car
                //so I don't have to fly all the way over there every time!
                camPos = carPos + vec3(5, 5, 5);
            }
            captured = true;
            SDL_SetRelativeMouseMode(SDL_TRUE);
        }



        //camera controls
        if (captured) {
            //first-person debug cam
            cameraYaw += dmx * 0.01f;
            cameraPitch += dmy * 0.01f;

            //generate 3D unit vector representing camera facing direction
            Vec3 facing = normalize({ sinf(cameraYaw), -sinf(cameraPitch), -cosf(cameraYaw) });

            float sp = dt * 10;
            if (isDown[SDL_SCANCODE_LSHIFT]) {
                sp *= 0.25f;
            }

            if (isDown[SDL_SCANCODE_A]) {
                camPos.x += facing.z * sp;
                camPos.z -= facing.x * sp;
            }
            if (isDown[SDL_SCANCODE_D]) {
                camPos.x -= facing.z * sp;
                camPos.z += facing.x * sp;
            }
            if (isDown[SDL_SCANCODE_W]) {
                camPos.x += facing.x * sp;
                camPos.y += facing.y * sp;
                camPos.z += facing.z * sp;
            }
            if (isDown[SDL_SCANCODE_S]) {
                camPos.x -= facing.x * sp;
                camPos.y -= facing.y * sp;
                camPos.z -= facing.z * sp;
            }
            if (isDown[SDL_SCANCODE_UP]) { //if up is down, something is very wrong... or not
                camPos.y += sp;
            }
            if (isDown[SDL_SCANCODE_DOWN]) { //BUG: this should always return true, but doesn't
                camPos.y -= sp;
            }
        } else {
            //follow cam
            float dist = 30;
            camPos = carPos + vec3(dist, dist * 1.0f, dist);
        }

        //lock camera pitch to within reasonable bounds
        if (cameraPitch > radians(80)) {
            cameraPitch = radians(80);
        } else if (cameraPitch < radians(-80)) {
            cameraPitch = radians(-80);
        }



        SDL_LockSurface(canvas);

        //clear the screen
        for (int y = 0; y < HEIGHT; ++y) {
            Pixel * row = (Pixel *) ((u8 *)canvas->pixels + y * canvas->pitch);
            float * zrow = depthBuffer + y * WIDTH;
            for (int x = 0; x < WIDTH; ++x) {
                row[x] = { 63, 191, 191, 255 };
                // row[x] = { 0, 0, 0, 255 };
                zrow[x] = 1.0f;
                // row[x] = { (u8)x, (u8)y, (u8)(x * y), 255 }; //test pattern
            }
        }

        //clear the shadow buffer
        for (int y = 0; y < SHADOW_HEIGHT; ++y) {
            float * row = shadowBuffer + y * SHADOW_WIDTH;
            for (int x = 0; x < SHADOW_WIDTH; ++x) {
                //arbitrarily high number (we don't clip the shadow buffer to [0, 1])
                row[x] = 1000000.0f;
            }
        }



        //clear the render queue
        drawList.len = 0;

        Vec3 lightPos = { -sinf(get_time() * 0.9f) * 5, 3, -cosf(get_time() * 0.9f) * 2 };
        // Vec3 lightPos = { 100, 100, 100 }; //out of the way

        //place things in the render queue
        {
            Mat4 modelMatrix = IDENTITY_4;
            // modelMatrix = translate(modelMatrix, Vec4{ sinf(get_time() * 0.1f) * 10, 0, 0 });
            // modelMatrix = rotateX(modelMatrix, get_time());
            // modelMatrix = rotateY(modelMatrix, get_time() * 0.5f);
            // modelMatrix = rotateZ(modelMatrix, get_time() * 0.25f);

            float deadZone = 0.3f;
            float xaxis = gamepads[0].axis[SDL_CONTROLLER_AXIS_LEFTX];
            float yaxis = gamepads[0].axis[SDL_CONTROLLER_AXIS_LEFTY];
            if (len(xaxis, yaxis) > deadZone) {
                carYaw = atan2f(xaxis, yaxis) + radians(45);
            }

            modelMatrix = rotateY(modelMatrix, carYaw);

            drawList.add({ carModel, modelMatrix, true });
        }

        {
            drawList.add({ floorModel, IDENTITY_4, true });
        }

        for (int i = 0; i < 100; ++i) {
            Mat4 matrix = IDENTITY_4;
            matrix = scale(matrix, 0.5f);
            matrix = rotateY(matrix, tablePos[i].w);
            matrix = translate(matrix, tablePos[i]);

            drawList.add({ tableModel, matrix, true });
        }

        {
            Mat4 ground = IDENTITY_4;
            ground = scale(ground, 1);
            ground = translate(ground, vec3(0, -5, 0));

            drawList.add({ tableModel, ground, false });
        }

        {
            Mat4 lightMatrix = IDENTITY_4;
            lightMatrix = scale(lightMatrix, 0.25f);
            lightMatrix = translate(lightMatrix, lightPos);

            drawList.add({ droidModel, lightMatrix, false });
        }



        //determine shadow direction
        float shadowPitch = radians(-45);
        float shadowYaw = get_time() * 0.1f;
        // float shadowYaw = radians(45);

        //setup shadow matrix
        Mat4 shadowMatrix = IDENTITY_4;
        shadowMatrix = rotateY(shadowMatrix, shadowYaw);
        shadowMatrix = rotateX(shadowMatrix, shadowPitch);
        float shadowScale = 1.0f / 10;
        shadowMatrix = scale(shadowMatrix, shadowScale, shadowScale, 1);

        //render to shadow buffer
        for (int i = 0; i < drawList.len; ++i) {
            if (drawList[i].shadow) {
                draw_static_model_shadow(drawList[i].model, drawList[i].matrix, shadowMatrix,
                    shadowBuffer);
            }
        }

        //apply view transforms
        Mat4 view = IDENTITY_4;
        Mat4 projection = IDENTITY_4;
        if (captured) {
            view = translate(view, vec3(-camPos.x, -camPos.y, -camPos.z));
            view = rotateY(view, cameraYaw);
            view = rotateX(view, cameraPitch);

            projection = perspective(radians(60), (float)WIDTH / (float)HEIGHT, 1, 100);
        } else {
            view = look_at(camPos, carPos, vec3(0, 1, 0));

            projection = perspective(radians(20), (float)WIDTH / (float)HEIGHT, 2, 200);
        }

        //setup lights
        Vec3 lightDir = {  sinf(shadowYaw),
                          -sinf(shadowPitch),
                          -cosf(shadowYaw) };
        lightDir = normalize(lightDir);

        //render to frame buffer
        for (int i = 0; i < drawList.len; ++i) {
            draw_static_model(drawList[i].model,
                drawList[i].matrix, view, projection, shadowMatrix,
                canvas, depthBuffer, shadowBuffer,
                lightPos, lightDir, camPos, materials);
        }



        //draw debug GUI
        if (!captured) {
            int mx, my;
            /*u32 mbflags =*/ SDL_GetMouseState(&mx, &my);
            // draw_line(canvas, WIDTH/2, HEIGHT/2, mx / SCALE, my / SCALE, { 255, 0, 255, 255 });
            // draw_point(canvas, mx / SCALE, my / SCALE, 5, { 255, 0, 255, 255 });
        }

        Mat4 viewproj = projection * view;
        draw_line(canvas, viewproj, 0, 0, 0, 1, 0, 0, { 255, 0, 0, 255 });
        draw_line(canvas, viewproj, 0, 0, 0, 0, 1, 0, { 0, 255, 0, 255 });
        draw_line(canvas, viewproj, 0, 0, 0, 0, 0, 1, { 0, 0, 255, 255 });
        draw_point(canvas, viewproj, 1, 0, 0, 1, { 255, 0, 0, 255 });
        draw_point(canvas, viewproj, 0, 1, 0, 1, { 0, 255, 0, 255 });
        draw_point(canvas, viewproj, 0, 0, 1, 1, { 0, 0, 255, 255 });

        SDL_UnlockSurface(canvas);



        //update sliding window filter for framerate
        float timeSum = 0;
        for (int i = 1; i < ARR_SIZE(frameTimes); ++i) {
            frameTimes[i - 1] = frameTimes[i];
            timeSum += frameTimes[i - 1];
        }
        frameTimes[ARR_SIZE(frameTimes) - 1] = dt;
        timeSum += dt;

        float framerate = ARR_SIZE(frameTimes) / timeSum;
        if (frame % 100 == 99) {
            // printf("%d fps\n", (int)(framerate + 0.5f));
            printf("total: %6d   drawn: %6d   shadow total: %6d   shadow: %6d"
                   "   draw calls: %6d   fps: %3d\n",
                   totalTris, drawnTris, shadowTotalTris, shadowDrawnTris,
                   (int) drawList.len, (int)(framerate + 0.5f));
        }

        totalTris = drawnTris = shadowTotalTris = shadowDrawnTris = 0;



        SDL_Surface * surface = SDL_GetWindowSurface(window);
        SDL_BlitScaled(canvas, nullptr, surface, nullptr);



        //DEBUG: draw shadow map
        for (int y = 0; y < SHADOW_HEIGHT; ++y) {
            Pixel * row = (Pixel *) ((u8 *)debugShadow->pixels + y * debugShadow->pitch);
            float * srow = shadowBuffer + y * SHADOW_WIDTH;
            for (int x = 0; x < SHADOW_WIDTH; ++x) {
                u8 s = srow[x] * 256;
                row[x] = { s, s, s, 255 };
            }
        }
        SDL_Rect blitDest = { 0, 0, 256, 256 };
        SDL_BlitScaled(debugShadow, nullptr, surface, &blitDest);



        SDL_UpdateWindowSurface(window);

        fflush(stdout);
        fflush(stderr);

        frame += 1;

        globalFrame = frame;

        //uncomment this to make the game exit immediately (good for testing compile+load times)
        // shouldExit = true;
    }







    printf("time alive: %f seconds\n", get_time());

    //free resources (not strictly necessary)
    drawList.finalize();

    free_model(carModel);
    free_model(tableModel);
    free_model(droidModel);

    SDL_FreeSurface(debugShadow);
    SDL_FreeSurface(canvas);
    SDL_DestroyWindow(window);
    SDL_Quit();

    printf("shutdown: %f seconds\n", get_time());

    return 0;
}
