#include <stdint.h>  /* Declarations of uint_32 and the like */
#include <stdlib.h>
#include <stdio.h>
#include <pic32mx.h> /* Declarations of system-specific addresses etc */
#include "header.h" /* Declatations for these labs */

#define TABLE_SIZE 720  // One entry per degree
#define PI 3.14159265358979323846
#define GRID_WIDTH 128
#define GRID_HEIGHT 32
#define MAX_TRIANGLES 90
float sine_table[TABLE_SIZE];
float cosine_table[TABLE_SIZE];

float fTheta=1;
float fYaw;
struct vec3d vLookDir = {0.0f,0.0f,0.0f,1.0f};
struct vec3d vCamera = {0.0f,0.0f,0.0f,1.0f};

int maxTriangles = 0;



uint8_t icon1[GRID_HEIGHT*GRID_WIDTH/8];

/**
 * @brief Initializes the icon array with a value of 0xFF.
 * 
 * This function iterates over the icon array and sets each element to 0xFF to turn screen off.
 * The icon array is the size of GRID_WIDTH*4 since the screen is split in 4 rows of 8bit numbers where each bit represent a pixel.
 * 
 * @param icons The icon array to be initialized.
 */
void initializeIcon(uint8_t icons[])
{
    int i;
    for (i = 0; i < GRID_WIDTH*4; i++)
    {
        icons[i] = 0xFF;
    }
}


/**
 * @brief Sets a pixel at the specified coordinates on the icon.
 * 
 * This function calculates the index for the byte in the array based on the coordinates (x, y).
 * It then finds the correct bit in the byte to change and clears that bit.
 * 
 * @param icon1 Pointer to the icon array.
 * @param x The x-coordinate of the pixel.
 * @param y The y-coordinate of the pixel.
 */
void setPixel(uint8_t *icon1, int x, int y)
{
    // Calculate index for the byte in the array
    int index = (y / 8) * GRID_WIDTH + x;

    // Find the correct bit in the byte to change
    int bit = y % 8;

    // Clear the bit
    icon1[index] &= ~(1 << bit);
}

/**
 * @brief Draws a line on the icon using Bresenham's line algorithm.
 * 
 * This function draws a line from the coordinates (x0, y0) to (x1, y1) on the icon.
 * It uses Bresenham's line algorithm to determine the pixels to set along the line.
 * 
 * @param icon Pointer to the icon array.
 * @param x0 The x-coordinate of the starting point of the line.
 * @param y0 The y-coordinate of the starting point of the line.
 * @param x1 The x-coordinate of the ending point of the line.
 * @param y1 The y-coordinate of the ending point of the line.
 */
void drawLine(uint8_t* icon, int x0, int y0, int x1, int y1) {
    int dx = abs(x1 - x0), sx = x0 < x1 ? 1 : -1;
    int dy = -abs(y1 - y0), sy = y0 < y1 ? 1 : -1; 
    int err = dx + dy, e2; /* error value e_xy */

    while (1) {
        setPixel(icon, x0, y0); // Set the pixel at (x0, y0) on the icon
        if (x0 == x1 && y0 == y1) break; // Exit the loop if (x0, y0) reaches (x1, y1)
        e2 = 2 * err;
        if (e2 >= dy) { /* e_xy+e_x > 0 */
            err += dy; x0 += sx; // Increment x0 and update the error value
        }
        if (e2 <= dx) { /* e_xy+e_y < 0 */
            err += dx; y0 += sy; // Increment y0 and update the error value
        }
    }
}

/**
 * @brief Draws a triangle on the icon.
 * 
 * This function draws a triangle on the icon by drawing three lines connecting the specified vertices.
 * 
 * @param icon Pointer to the icon array.
 * @param x0 The x-coordinate of the first vertex.
 * @param y0 The y-coordinate of the first vertex.
 * @param x1 The x-coordinate of the second vertex.
 * @param y1 The y-coordinate of the second vertex.
 * @param x2 The x-coordinate of the third vertex.
 * @param y2 The y-coordinate of the third vertex.
 */
void drawTriangle(uint8_t* icon, int x0, int y0, int x1, int y1, int x2, int y2){
    drawLine(icon, x0, y0, x1, y1);
    drawLine(icon, x0, y0, x2, y2);
    drawLine(icon, x1, y1, x2, y2);
}

/**
 * Initializes a vector with the given coordinates.
 *
 * @param vec Pointer to the vector to be initialized.
 * @param x The x-coordinate of the vector.
 * @param y The y-coordinate of the vector.
 * @param z The z-coordinate of the vector.
 * @param w The w-coordinate of the vector.
 */
void InitializeVector(struct vec3d* vec, float x, float y, float z, float w) {
    vec->x = x;
    vec->y = y;
    vec->z = z;
    vec->w = w;
}

/**
 * Initializes a triangle with the given coordinates for its vertices.
 *
 * @param tri Pointer to the triangle to be initialized.
 * @param x1 The x-coordinate of the first vertex.
 * @param y1 The y-coordinate of the first vertex.
 * @param z1 The z-coordinate of the first vertex.
 * @param w1 The w-coordinate of the first vertex.
 * @param x2 The x-coordinate of the second vertex.
 * @param y2 The y-coordinate of the second vertex.
 * @param z2 The z-coordinate of the second vertex.
 * @param w2 The w-coordinate of the second vertex.
 * @param x3 The x-coordinate of the third vertex.
 * @param y3 The y-coordinate of the third vertex.
 * @param z3 The z-coordinate of the third vertex.
 * @param w3 The w-coordinate of the third vertex.
 */
void InitializeTriangle(struct triangle* tri, float x1, float y1, float z1, float w1,
                        float x2, float y2, float z2, float w2,
                        float x3, float y3, float z3, float w3) {
    InitializeVector(&tri->p[0], x1, y1, z1, w1);
    InitializeVector(&tri->p[1], x2, y2, z2, w2);
    InitializeVector(&tri->p[2], x3, y3, z3, w3);
}

/**
 * Initializes a triangle with zero coordinates for its vertices.
 *
 * @param tri Pointer to the triangle to be initialized.
 */
void InitializeTriangle2(struct triangle* tri) {
    InitializeVector(&tri->p[0], 0.0f, 0.0f, 0.0f, 1.0f);
    InitializeVector(&tri->p[1], 0.0f, 0.0f, 0.0f, 1.0f);
    InitializeVector(&tri->p[2], 0.0f, 0.0f, 0.0f, 1.0f);
}

/**
 * Pushes a triangle into an array of triangles.
 *
 * @param array The array of triangles.
 * @param size The current size of the array.
 * @param newTriangle The triangle to be pushed into the array.
 */
void push_triangle(struct triangle* array, int* size, struct triangle newTriangle) {
    if (*size < MAX_TRIANGLES) {
        array[*size] = newTriangle;
        (*size)++;
    } else {
        // Handle the case when the array is full
    }
}

/**
 * Pops a triangle from an array of triangles.
 *
 * @param array The array of triangles.
 * @param size The current size of the array.
 * @param poppedTriangle A pointer to store the popped triangle.
 */
void pop_triangle(struct triangle* array, int* size, struct triangle* poppedTriangle) {
    if (*size > 0) {
        // Shift all elements one position to the left
        *poppedTriangle = array[0];
        int i;
        for (i = 0; i < *size - 1; i++) {
            array[i] = array[i + 1];
        }
        (*size)--; // Decrease the size of the array
    } else {
        // Handle the case when the array is empty
    }
}

/**
 * Compares two triangles based on their average z-coordinate.
 *
 * @param a A pointer to the first triangle.
 * @param b A pointer to the second triangle.
 * @return -1 if the average z-coordinate of the first triangle is greater than the second triangle,
 *          1 if the average z-coordinate of the first triangle is less than the second triangle,
 *          0 if the average z-coordinate of both triangles are equal.
 */
int compareTriangles(const void* a, const void* b) {
    const struct triangle* t1 = (const struct triangle*)a;
    const struct triangle* t2 = (const struct triangle*)b;
    float z1 = (t1->p[0].z + t1->p[1].z + t1->p[2].z) / 3.0f;
    float z2 = (t2->p[0].z + t2->p[1].z + t2->p[2].z) / 3.0f;
    if (z1 > z2) return -1;
    if (z1 < z2) return 1;
    return 0;
}

struct mat4x4 matProj; // Use the mat4x4 type directly
struct mesh meshCube;

/**
 * Initializes the rendering for the application.
 * 
 * This function sets up the meshCube and initializes the projection matrix for a 1x1 cube.
 * It also manually copies the initial triangles to the allocated memory.
 * Finally, it creates a perspective projection matrix using the given parameters.
 * 
 * @return 1 if the initialization is successful, 0 otherwise.
 */
int OnUserCreate() {
    meshCube.size = 12;
    InitializeMatrix(&matProj, 0);

    // Initialize triangles (manually copy each triangle)
    struct triangle init_triangles[12] = {
        // SOUTH
        {
            .p[0] = {0.0f, 0.0f, 0.0f},
            .p[1] = {0.0f, 1.0f, 0.0f},
            .p[2] = {1.0f, 1.0f, 0.0f}
        },
        {
            .p[0] = {0.0f, 0.0f, 0.0f},
            .p[1] = {1.0f, 1.0f, 0.0f},
            .p[2] = {1.0f, 0.0f, 0.0f}
        },

        // EAST
        {
            .p[0] = {1.0f, 0.0f, 0.0f},
            .p[1] = {1.0f, 1.0f, 0.0f},
            .p[2] = {1.0f, 1.0f, 1.0f}
        },
        {
            .p[0] = {1.0f, 0.0f, 0.0f},
            .p[1] = {1.0f, 1.0f, 1.0f},
            .p[2] = {1.0f, 0.0f, 1.0f}
        },

        // NORTH
        {
            .p[0] = {1.0f, 0.0f, 1.0f},
            .p[1] = {1.0f, 1.0f, 1.0f},
            .p[2] = {0.0f, 1.0f, 1.0f}
        },
        {
            .p[0] = {1.0f, 0.0f, 1.0f},
            .p[1] = {0.0f, 1.0f, 1.0f},
            .p[2] = {0.0f, 0.0f, 1.0f}
        },

        // WEST
        {
            .p[0] = {0.0f, 0.0f, 1.0f},
            .p[1] = {0.0f, 1.0f, 1.0f},
            .p[2] = {0.0f, 1.0f, 0.0f}
        },
        {
            .p[0] = {0.0f, 0.0f, 1.0f},
            .p[1] = {0.0f, 1.0f, 0.0f},
            .p[2] = {0.0f, 0.0f, 0.0f}
        },

        // TOP
        {
            .p[0] = {0.0f, 1.0f, 0.0f},
            .p[1] = {0.0f, 1.0f, 1.0f},
            .p[2] = {1.0f, 1.0f, 1.0f}
        },
        {
            .p[0] = {0.0f, 1.0f, 0.0f},
            .p[1] = {1.0f, 1.0f, 1.0f},
            .p[2] = {1.0f, 1.0f, 0.0f}
        },

        // BOTTOM
        {
            .p[0] = {1.0f, 0.0f, 1.0f},
            .p[1] = {0.0f, 0.0f, 1.0f},
            .p[2] = {0.0f, 0.0f, 0.0f}
        },
        {
            .p[0] = {1.0f, 0.0f, 1.0f},
            .p[1] = {0.0f, 0.0f, 0.0f},
            .p[2] = {1.0f, 0.0f, 0.0f}
        }
        
    };



    // Copy the initialized triangles to the allocated memory
    int i;
    for (i = 0; i < meshCube.size; ++i) {
        meshCube.p[i] = init_triangles[i];
    }

    // Creates a perspective projection matrix
    float fNear = 0.1f;
    float fFar = 1000.0f;
    float fFov = 45.0f;
    float fAspectRatio = GRID_HEIGHT / (float)GRID_WIDTH;

    InitializeMatrix(&matProj,0);
    Matrix_MakeProjection(&matProj,fFov, fAspectRatio, fNear, fFar);
    
    return 1;
}

/**
 * 
 * This function creates a pyramid mesh and sets up the perspective projection matrix.
 * 
 * @return 1 if the initialization is successful, otherwise 0.
 */
int OnUserCreate2() {
    // Create pyramid
    meshCube.size = 6;
    
    // Initialize triangles for pyramid
    struct triangle init_triangles2[6] = {
        // Base
        {
            .p[0] = {2.0f, 0.0f, 2.0f},
            .p[1] = {0.0f, 0.0f, 2.0f},
            .p[2] = {0.0f, 0.0f, 0.0f}
        },
        {
            .p[0] = {2.0f, 0.0f, 2.0f},
            .p[1] = {0.0f, 0.0f, 0.0f},
            .p[2] = {2.0f, 0.0f, 0.0f}
        },

        // Sides
        {
            .p[0] = {0.0f, 0.0f, 0.0f},
            .p[1] = {1.0f, 1.0f, 1.0f},
            .p[2] = {2.0f, 0.0f, 0.0f}
        },
        {
            .p[0] = {2.0f, 0.0f, 0.0f},
            .p[1] = {1.0f, 1.0f, 1.0f},
            .p[2] = {2.0f, 0.0f, 2.0f}
        },
        {
            .p[0] = {2.0f, 0.0f, 2.0f},
            .p[1] = {1.0f, 1.0f, 1.0f},
            .p[2] = {0.0f, 0.0f, 2.0f}
        },
        {
            .p[0] = {0.0f, 0.0f, 2.0f},
            .p[1] = {1.0f, 1.0f, 1.0f},
            .p[2] = {0.0f, 0.0f, 0.0f}
        }
    };

    int b;
    for (b = 0; b < meshCube.size; ++b) {
        meshCube.p[b] = init_triangles2[b];
    }

    // Creates a perspective projection matrix
    float fNear = 0.1f;
    float fFar = 1000.0f;
    float fFov = 45.0f;
    float fAspectRatio = GRID_HEIGHT / (float)GRID_WIDTH;

    InitializeMatrix(&matProj, 0);
    Matrix_MakeProjection(&matProj, fFov, fAspectRatio, fNear, fFar);
    
    return 1;
}
/**
 * Initializes the cube mesh for a block and creates a perspective projection matrix.
 * 
 * @return 1 if successful, otherwise 0.
 */
int OnUserCreate3() {
    meshCube.size = 12;
    InitializeMatrix(&matProj, 0);

    // Initialize triangles (manually copy each triangle)
    struct triangle init_triangles[12] = {
        // SOUTH
        {
            .p[0] = {0.0f, 0.0f, 0.0f},
            .p[1] = {0.0f, 2.0f, 0.0f},
            .p[2] = {1.0f, 2.0f, 0.0f}
        },
        {
            .p[0] = {0.0f, 0.0f, 0.0f},
            .p[1] = {1.0f, 2.0f, 0.0f},
            .p[2] = {1.0f, 0.0f, 0.0f}
        },

        // EAST
        {
            .p[0] = {1.0f, 0.0f, 0.0f},
            .p[1] = {1.0f, 2.0f, 0.0f},
            .p[2] = {1.0f, 2.0f, 1.0f}
        },
        {
            .p[0] = {1.0f, 0.0f, 0.0f},
            .p[1] = {1.0f, 2.0f, 1.0f},
            .p[2] = {1.0f, 0.0f, 1.0f}
        },

        // NORTH
        {
            .p[0] = {1.0f, 0.0f, 1.0f},
            .p[1] = {1.0f, 2.0f, 1.0f},
            .p[2] = {0.0f, 2.0f, 1.0f}
        },
        {
            .p[0] = {1.0f, 0.0f, 1.0f},
            .p[1] = {0.0f, 2.0f, 1.0f},
            .p[2] = {0.0f, 0.0f, 1.0f}
        },

        // WEST
        {
            .p[0] = {0.0f, 0.0f, 1.0f},
            .p[1] = {0.0f, 2.0f, 1.0f},
            .p[2] = {0.0f, 2.0f, 0.0f}
        },
        {
            .p[0] = {0.0f, 0.0f, 1.0f},
            .p[1] = {0.0f, 2.0f, 0.0f},
            .p[2] = {0.0f, 0.0f, 0.0f}
        },

        // TOP
        {
            .p[0] = {0.0f, 2.0f, 0.0f},
            .p[1] = {0.0f, 2.0f, 1.0f},
            .p[2] = {1.0f, 2.0f, 1.0f}
        },
        {
            .p[0] = {0.0f, 2.0f, 0.0f},
            .p[1] = {1.0f, 2.0f, 1.0f},
            .p[2] = {1.0f, 2.0f, 0.0f}
        },

        // BOTTOM
        {
            .p[0] = {1.0f, 0.0f, 1.0f},
            .p[1] = {0.0f, 0.0f, 1.0f},
            .p[2] = {0.0f, 0.0f, 0.0f}
        },
        {
            .p[0] = {1.0f, 0.0f, 1.0f},
            .p[1] = {0.0f, 0.0f, 0.0f},
            .p[2] = {1.0f, 0.0f, 0.0f}
        }
    };
        // Copy the initialized triangles to the allocated memory
    int i;
    for (i = 0; i < meshCube.size; ++i) {
        meshCube.p[i] = init_triangles[i];
    }

    // Creates a perspective projection matrix
    float fNear = 0.1f;
    float fFar = 1000.0f;
    float fFov = 45.0f;
    float fAspectRatio = GRID_HEIGHT / (float)GRID_WIDTH;

    InitializeMatrix(&matProj,0);
    Matrix_MakeProjection(&matProj,fFov, fAspectRatio, fNear, fFar);
    
    return 1;
}
/**
 * Initializes the cube mesh for a triangle on top of another triangle and creates a perspective projection matrix.
 * 
 * @return 1 if successful, otherwise 0.
 */
int OnUserCreate4() {
    // Create pyramid
    meshCube.size = 12;
    // Initialize triangles for pyramid
    struct triangle init_triangles2[12] = {
        //Base
        {
            .p[0] = {1.0f, 0.0f, 1.0f},
            .p[1] = {0.0f, 0.0f, 1.0f},
            .p[2] = {0.0f, 0.0f, 0.0f}
        },
        {
            .p[0] = {1.0f, 0.0f, 1.0f},
            .p[1] = {0.0f, 0.0f, 0.0f},
            .p[2] = {1.0f, 0.0f, 0.0f}
        },

        // Sides
        {
            .p[0] = {0.0f, 0.0f, 0.0f},
            .p[1] = {1.0f, 1.0f, 1.0f},
            .p[2] = {1.0f, 0.0f, 0.0f}
        },
        {
            .p[0] = {1.0f, 0.0f, 0.0f},
            .p[1] = {1.0f, 1.0f, 1.0f},
            .p[2] = {1.0f, 0.0f, 1.0f}
        },
        {
            .p[0] = {1.0f, 0.0f, 1.0f},
            .p[1] = {1.0f, 1.0f, 1.0f},
            .p[2] = {0.0f, 0.0f, 1.0f}
        },
        {
            .p[0] = {0.0f, 0.0f, 1.0f},
            .p[1] = {1.0f, 1.0f, 1.0f},
            .p[2] = {0.0f, 0.0f, 0.0f}
        },

	//top
	    {
            .p[0] = {1.0f, 2.0f, 1.0f},
            .p[1] = {0.0f, 2.0f, 1.0f},
            .p[2] = {0.0f, 2.0f, 0.0f}
        },
        {
            .p[0] = {1.0f, 2.0f, 1.0f},
            .p[1] = {0.0f, 2.0f, 0.0f},
            .p[2] = {1.0f, 2.0f, 0.0f}
        },
	//sides
	    {
            .p[0] = {1.0f, 1.0f, 1.0f},
            .p[1] = {0.0f, 2.0f, 0.0f},
            .p[2] = {1.0f, 2.0f, 0.0f}
        },
        {
            .p[0] = {1.0f, 1.0f, 1.0f},
            .p[1] = {1.0f, 2.0f, 0.0f},
            .p[2] = {1.0f, 2.0f, 1.0f}
        },
        {
            .p[0] = {1.0f, 1.0f, 1.0f},
            .p[1] = {1.0f, 2.0f, 1.0f},
            .p[2] = {0.0f, 2.0f, 1.0f}
        },
        {
            .p[0] = {1.0f, 1.0f, 1.0f},
            .p[1] = {0.0f, 2.0f, 1.0f},
            .p[2] = {0.0f, 2.0f, 0.0f}
        }

    };

    int b;
        for (b = 0; b < meshCube.size; ++b) {
            meshCube.p[b] = init_triangles2[b];
    }

    // Creates a perspective projection matrix
    float fNear = 0.1f;
    float fFar = 1000.0f;
    float fFov = 45.0f;
    float fAspectRatio = GRID_HEIGHT / (float)GRID_WIDTH;

    InitializeMatrix(&matProj,0);
    Matrix_MakeProjection(&matProj,fFov, fAspectRatio, fNear, fFar);
    
    return 1;
}



// This function initializes the meshes based on the user's choice
// It calls the appropriate initialization functions based on the choice parameter
void labinit(int choice) 
{
    initializeIcon(icon1);
    initialize_tables();
    
    // Call the appropriate OnUserCreate function based on the choice parameter
    if (choice == 0) {
        OnUserCreate();
    }
    if (choice == 1) {
        OnUserCreate2();
    }
    if (choice == 2) {
        //----
        OnUserCreate3();
    }
    if (choice == 3) {
        OnUserCreate4();
    }
    
    inputInit();
    return;
}
/**
 * 
 * This function is called for each frame update and handles user input and rendering.
 *
 * Parameters: 
 *     fElapsedTime: The elapsed time since the last frame update.
 */
void OnUserUpdate(float fElapsedTime){
// Set up "World Tranmsform" though not updating theta 
		// makes this a bit redundant
        volatile unsigned int *pointerToPORTE;
        pointerToPORTE = (volatile unsigned int *)0xbf886110;
        unsigned int currBtns = getbtns();
        unsigned int currSw = getsw();
        if (currBtns ) {
            if (currBtns & 0x4)  // If bit number four is set, increment fSigma
            {
                vCamera.y += 0.2f;
            }
            if (currBtns & 0x2)  // If bit number two is set, increment fTheta
            {
                vCamera.y -= 0.2f;
            }
            if (currBtns & 0x1)  // If bit number one is set, decrement fSigma
            {
                vCamera.x -= 0.2f;
            }
            if (currBtns & 0x8)  // If bit number eight is set, decrement fTheta
            {
                vCamera.x += 0.2f;
            }
        }

        struct vec3d vForward; 
        InitializeVector(&vForward, 0.0f, 0.0f, 0.0f, 1.0f);
        Vector_Mul(&vForward, &vLookDir, 1.0f * fElapsedTime);
        if(fYaw > 2*PI){
            fYaw = 0;
        }

        if (currSw) {
            if (currSw & 0x4)  // If bit number four is set, increment ty
            {
                Vector_Add(&vCamera, &vCamera, &vForward);
            }
            if (currSw & 0x2)  // If bit number two is set, increment tx
            {
                Vector_Sub(&vCamera, &vCamera, &vForward);
            }
            if (currSw & 0x1)  // If bit number one is set, increment tz
            {
                fYaw -= 0.3f * fElapsedTime;
                if (fYaw < 0) {
                    fYaw += PI * 2.0f;
                }            }
            if (currSw & 0x8)  // If bit number eight is set, reset to default values
            {
                fYaw += 0.3f * fElapsedTime;
            }
        }

        struct mat4x4 matWorld;
        InitializeMatrix(&matWorld, 0);
        Matrix_MakeIdentity(&matWorld);
        SetWorldMatrix(&matWorld, fElapsedTime);
        // printmatrices(&matWorld, "matWorld1");
        

		// Create "Point At" Matrix for camera
        struct mat4x4 matView;
        InitializeMatrix(&matView, 0);

        struct vec3d vUp = { 0,1,0,1 };
        struct vec3d vTarget = { 0,0,1,1 };
        struct mat4x4 matCameraRot;
        InitializeMatrix(&matCameraRot,0);
        Matrix_MakeRotationY(&matCameraRot, fYaw);
        Matrix_MultiplyVector(&vTarget, &vLookDir,&matCameraRot);//man kan inte ändra const

        Vector_Add(&vTarget, &vCamera, &vLookDir);
        struct mat4x4 matCamera;
        InitializeMatrix(&matCamera, 0);
        Matrix_PointAt(&matCamera, &vCamera, &vTarget, &vUp);

        // Make view matrix from camera
        Matrix_QuickInverse(&matView, &matCamera);

        // printmatrices(&matView, "matView");
		// Store triagles for rastering later
        struct triangle vecTrianglesToRaster[meshCube.size];
        int ii;
        for (ii = 0; ii < meshCube.size; ++ii) {
            InitializeTriangle2(&vecTrianglesToRaster[ii]);
        }
        int sizeVecRaster = 0;
		// Draw Triangles

    unsigned int i;
    for (i=0; i < meshCube.size; ++i) {
            struct triangle triProjected, triTransformed, triViewed;
            InitializeTriangle2(&triProjected);
            InitializeTriangle2(&triTransformed);
            InitializeTriangle2(&triViewed);
			// World Matrix Transform
            // printvec(&meshCube.p[i].p[0], "meshCube.p[i].p[0]");
            // printvec(&meshCube.p[i].p[1], "meshCube.p[i].p[1]");
            // printvec(&meshCube.p[i].p[2], "meshCube.p[i].p[2]");
            // printmatrices(&matWorld, "matWorldMultiply");
			Matrix_MultiplyVector( &meshCube.p[i].p[0], &triTransformed.p[0], &matWorld);
			Matrix_MultiplyVector( &meshCube.p[i].p[1], &triTransformed.p[1], &matWorld);
			Matrix_MultiplyVector( &meshCube.p[i].p[2], &triTransformed.p[2], &matWorld);
            // printvec(&triTransformed.p[0], "triTransformed.p[0]");
            // printvec(&triTransformed.p[1], "triTransformed.p[1]");
            // printvec(&triTransformed.p[2], "triTransformed.p[2]");

			// Calculate triangle Normal
			struct vec3d normal= {0.0f,0.0f,0.0f,1.0f}, line1= {0.0f,0.0f,0.0f,1.0f}, line2= {0.0f,0.0f,0.0f,1.0f};

			// Get lines either side of triangle
			Vector_Sub(&line1, &triTransformed.p[1], &triTransformed.p[0]);
			Vector_Sub(&line2, &triTransformed.p[2], &triTransformed.p[0]);
            // printvec(&line1, "line1");
            // printvec(&line2, "line2");
			// Take cross product of lines to get normal to triangle surface
			Vector_CrossProduct(&normal, &line1, &line2);
            // printvec(&normal, "normalCross");

			// You normally need to normalise a normal!
			Vector_Normalise(&normal, &normal);
            // printvec(&normal, "normalNorm");
			// Get Ray from triangle to camera
            struct vec3d vCameraRay = {0.0f,0.0f,0.0f,1.0f};
			Vector_Sub(&vCameraRay, &triTransformed.p[0], &vCamera);
            // printvec(&vCameraRay, "vCameraRay");
            
            // printf("DotProduct: %f\n", Vector_DotProduct(&normal, &vCameraRay));
            if (Vector_DotProduct(&normal, &vCameraRay) < 0.0f)
			{
				// Convert World Space --> View Space
				Matrix_MultiplyVector( &triTransformed.p[0], &triViewed.p[0], &matView);
				Matrix_MultiplyVector( &triTransformed.p[1], &triViewed.p[1], &matView);
				Matrix_MultiplyVector( &triTransformed.p[2], &triViewed.p[2], &matView);
                // printvec(&triViewed.p[0], "triViewed.p[0]");
                // printvec(&triViewed.p[1], "triViewed.p[1]");
                // printvec(&triViewed.p[2], "triViewed.p[2]");


				int nClippedTriangles = 0;
                struct triangle clipped[2];
                //initialize clipped
                InitializeTriangle2(&clipped[0]);
                InitializeTriangle2(&clipped[1]);


                struct vec3d plane_p = { 0.0f, 0.0f, 0.1f, 1.0f };
                struct vec3d plane_n = { 0.0f, 0.0f, 1.0f, 1.0f };
                
				nClippedTriangles = Triangle_ClipAgainstPlane(&plane_p, &plane_n, &triViewed, &clipped[0], &clipped[1]);
                // printf("nClippedTriangles: %d\n", nClippedTriangles);

                int n;
				for (n = 0; n < nClippedTriangles; n++)
				{
					// Project triangles from 3D --> 2D
					Matrix_MultiplyVector(&clipped[n].p[0], &triProjected.p[0], &matProj);
                    Matrix_MultiplyVector(&clipped[n].p[1], &triProjected.p[1], &matProj);
                    Matrix_MultiplyVector(&clipped[n].p[2], &triProjected.p[2], &matProj);
                    // printvec(&triProjected.p[0], "CliptriProjected.p[0]");
                    // printvec(&triProjected.p[1], "triProjected.p[1]");
                    // printvec(&triProjected.p[2], "triProjected.p[2]");



					Vector_Div(&triProjected.p[0], &triProjected.p[0], triProjected.p[0].w);
					Vector_Div(&triProjected.p[1], &triProjected.p[1], triProjected.p[1].w);
					Vector_Div(&triProjected.p[2], &triProjected.p[2], triProjected.p[2].w);
                    // printvec(&triProjected.p[0], "DivtriProjected.p[0]");
                    // printvec(&triProjected.p[1], "DivtriProjected.p[1]");
                    // printvec(&triProjected.p[2], "DivtriProjected.p[2]");


					triProjected.p[0].x *= -1.0f;
					triProjected.p[1].x *= -1.0f;
					triProjected.p[2].x *= -1.0f;
					triProjected.p[0].y *= -1.0f;
					triProjected.p[1].y *= -1.0f;
					triProjected.p[2].y *= -1.0f;


					struct vec3d vOffsetView = { 1,1,0,1 };
					Vector_Add(&triProjected.p[0], &triProjected.p[0], &vOffsetView);
					Vector_Add(&triProjected.p[1], &triProjected.p[1], &vOffsetView);
					Vector_Add(&triProjected.p[2], &triProjected.p[2], &vOffsetView);
                    // printvec(&triProjected.p[0], "OffsettriProjected.p[0]");
                    // printvec(&triProjected.p[1], "OffsettriProjected.p[1]");
                    // printvec(&triProjected.p[2], "OffsettriProjected.p[2]");
					triProjected.p[0].x *= 0.5f * (float)GRID_WIDTH;
					triProjected.p[0].y *= 0.5f * (float)GRID_HEIGHT;
					triProjected.p[1].x *= 0.5f * (float)GRID_WIDTH;
					triProjected.p[1].y *= 0.5f * (float)GRID_HEIGHT;
					triProjected.p[2].x *= 0.5f * (float)GRID_WIDTH;
					triProjected.p[2].y *= 0.5f * (float)GRID_HEIGHT;
                    // printvec(&triProjected.p[0], "ScaletriProjected.p[0]");
                    // printvec(&triProjected.p[1], "ScaletriProjected.p[1]");
                    // printvec(&triProjected.p[2], "ScaletriProjected.p[2]");
                    push_triangle(triangles, &sizeVecRaster, triProjected);
				}			
                //print pushed triangles in triangles
                struct triangle triToRaster3;
                InitializeTriangle2(&triToRaster3);
                int m;
                for (m = 0; m < sizeVecRaster; m++)
                {
                    triToRaster3 = triangles[m];
                    // printf("thisisitbefore\n  ");
                    // printf("Triangle %d:\n", m + 1);
                    // printTriangle(&triToRaster3);
                }
			}
        }
        //print triangles
		// Loop through all transformed, viewed, projected, and sorted triangles
        struct triangle triToRaster;
        InitializeTriangle2(&triToRaster);
        qsort(triangles, sizeVecRaster, sizeof(struct triangle), compareTriangles);
        //print qsort triangles in triangles
        int m;
        for (m = 0; m < sizeVecRaster; m++)
        {
            triToRaster = triangles[m];
            // printf("thisisitqsort\n  ");
            // printf("Triangle %d:\n", m + 1);
            // printTriangle(&triToRaster);
        }

        int sizelistTriangles = 0;
        int j;
		for (j=0; j< sizeVecRaster; ++j)//12 för size av vecTriangleToRaster
		{
            

			struct triangle clipped[2];
            //initialize clipped
            InitializeTriangle2(&clipped[0]);
            InitializeTriangle2(&clipped[1]);
			
            //initialize listTriangles
            struct triangle listTriangles[sizeVecRaster];
            //initialize listTriangles
            size_t i;
            for (i = 0; i < sizeVecRaster; i++)
            {
                InitializeTriangle2(&listTriangles[i]);
            }
            
            //från o me här
			// Add initial triangle through push
            push_triangle(listTriangles, &sizelistTriangles, triangles[j]);
            //print pushed triangles in triangles

            struct triangle triToRaster2;
            InitializeTriangle2(&triToRaster2);
            int m;
            for (m = 0; m < sizelistTriangles; m++)
            {
                triToRaster2 = listTriangles[m];
                // printf("thisisitafter\n  ");
                // printf("Triangle %d:\n", m + 1);
                // printTriangle(&triToRaster2);
            }
            //till o med här

			int nNewTriangles = 1;
            int k=0;
            int p;
			for (p = 0; p < 4; p++)
			{
				int nTrisToAdd = 0;
				while (nNewTriangles > 0)
				{
					// Take triangle from front of queue
                    struct triangle test;
                    InitializeTriangle2(&test);
                    pop_triangle(listTriangles, &sizelistTriangles, &test);
                    k++;
                    nNewTriangles--;
                    //print what pop gives (&test)
                    // printf("pop_triangle %d:\n", k);
                    // printTriangle(&test);


	
                    // printf("p: %d\n", p);
                    switch (p)
                    {
                        case 0:
                            {
                                struct vec3d planeNormal = { 0.0f, 1.0f, 0.0f, 1.0f };
                                struct vec3d planePoint = { 0.0f, 0.0f, 0.0f, 1.0f };
                                nTrisToAdd = Triangle_ClipAgainstPlane(&planePoint, &planeNormal, &test, &clipped[0], &clipped[1]);
                                // printf("nTrisToAdd2: %d\n", nTrisToAdd);
                                break;
                            }
                        case 1:
                            {
                                struct vec3d planeNormal = { 0.0f, -1.0f, 0.0f, 1.0f };
                                struct vec3d planePoint = { 0.0f, (float)GRID_HEIGHT - 1.0f, 0.0f, 1.0f };
                                nTrisToAdd = Triangle_ClipAgainstPlane(&planePoint, &planeNormal, &test, &clipped[0], &clipped[1]);
                                // printf("nTrisToAdd21: %d\n", nTrisToAdd);
                                break;
                            }
                        case 2:
                            {
                                struct vec3d planeNormal = { 1.0f, 0.0f, 0.0f, 1.0f };
                                struct vec3d planePoint = { 0.0f, 0.0f, 0.0f, 1.0f };
                                nTrisToAdd = Triangle_ClipAgainstPlane(&planePoint, &planeNormal, &test, &clipped[0], &clipped[1]);
                                break;
                            }
                        case 3:
                            {
                                struct vec3d planeNormal = { -1.0f, 0.0f, 0.0f, 1.0f };
                                struct vec3d planePoint = { (float)GRID_WIDTH-1.0f, 0.0f, 0.0f, 1.0f };
                                nTrisToAdd = Triangle_ClipAgainstPlane(&planePoint, &planeNormal, &test, &clipped[0], &clipped[1]);
                                break;
                            }
                    }


                    int m;
                    for (m = 0; m < nTrisToAdd; m++)
                    {
                        push_triangle(listTriangles, &sizelistTriangles, clipped[m]);
                        // printf("size: %d\n", sizelistTriangles);
                        // //print what has been pushed
                        // printf("push_triangle %d:\n", m);
                        // printTriangle(&clipped[m]);
                    }
    
				}
                nNewTriangles = sizelistTriangles;
            }
            for (i = 0; i < sizelistTriangles; i++)
            {
                triToRaster = listTriangles[i];
                // printf("thisisit\n  ");
                // printf("Triangle %d:\n", i + 1);
                // printTriangle(&triToRaster);
                drawTriangle(icon1, (int)triToRaster.p[0].x, (int)triToRaster.p[0].y, (int)triToRaster.p[1].x, (int)triToRaster.p[1].y, (int)triToRaster.p[2].x, (int)triToRaster.p[2].y);
            }
	    }
}

/**
 * Converts a floating-point number to a string representation with a specified number of decimal places.
 *
 * @param buffer The buffer to store the resulting string.
 * @param value The floating-point number to convert.
 * @param decimalPlaces The number of decimal places to include in the string.
 */
void floatToString(char *buffer, float value, int decimalPlaces) {
    if (decimalPlaces < 0) decimalPlaces = 0;

    // Handle negative numbers
    if (value < 0.0) {
        *buffer++ = '-';
        value = -value;
    }

    // Round the value based on the decimal places
    float rounding = 0.5;
    int i;
    for (i = 0; i < decimalPlaces; ++i) {
        rounding /= 10.0;
    }
    value += rounding;

    // Extract integer part
    long intPart = (long)value;
    value -= intPart;

    // Convert integer part to string
    char* intStr = buffer;
    if (intPart == 0) {
        *intStr++ = '0'; // Handle zero integer part
    } else {
        long temp = intPart;
        // Convert integer part to string in reverse order
        while (temp > 0) {
            *intStr++ = (char)(temp % 10 + '0');
            temp /= 10;
        }
    }

    // Reverse integer part string
    *intStr = '\0';
    char* end = intStr - 1;
    while (buffer < end) {
        char tempChar = *buffer;
        *buffer++ = *end;
        *end-- = tempChar;
    }

    // Add decimal point
    if (decimalPlaces > 0) {
        *intStr++ = '.';
    }

    // Convert fractional part
    while (decimalPlaces-- > 0) {
        value *= 10.0;
        int digit = (int)value;
        *intStr++ = (char)(digit + '0');
        value -= digit;
    }

    // Null terminate the string
    *intStr = '\0';
}




/* Interrupt Service Routine */
void user_isr( void )
{
  return;
}

float time=0.0f;
float fYaw=0.0f;
float counterTime2 = 0;
float timeInSeconds = 0;
float timebefore;
float elapsedTime;
void newlabwork(void)
{
    if( IFS(0) & (1<<8)){
        counterTime2 = 0.1;
        IFSCLR(0) = (1 << 8);  // Clear the Timer 2 interrupt flag
    }
    OnUserUpdate(counterTime2);
    display_update();
    display_image(0, icon1);
    initializeIcon(icon1);
}