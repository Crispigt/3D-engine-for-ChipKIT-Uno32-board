#include <stdint.h>  /* Declarations of uint_32 and the like */
#include <stdlib.h>
#include <stdio.h>
#include "header.h" /* Declatations for these labs */


// Function to calculate the absolute value of an integer
int abs(int x){
    if (x < 0){
       return x = -x;
    }
    else{
        return x;
    }
}

// Function to initialize a matrix with a given initial value
void InitializeMatrix(struct mat4x4 *matrix, float initialValue) {
    int i, j;
    for (i = 0; i < 4; ++i) {
        for (j = 0; j < 4; ++j) {
            matrix->m[i][j] = initialValue;
        }
    }
}

// Function to calculate the factorial of a number
float factorial(int n) {
    float result = 1.0f;
    int i;
    for (i = 2; i <= n; ++i) {
        result *= i;
    }
    return result;
}

// Function to calculate the power of a number
double pow1(double base, int exponent) {
    if (exponent == 0) return 1;  // Any number to the power of 0 is 1
    double result = 1;
    int i;
    for (i = 0; i < exponent; i++) {
        result *= base;
    }
    return result;
}
/**
 * Calculates the absolute value of a double number.
 *
 * @param x The number to calculate the absolute value of.
 * @return The absolute value of the input number.
 */
double my_fabs(double x) {
    if (x < 0) {
        return -x; // If x is negative, return its negation (making it positive)
    }
    return x; // If x is positive or zero, return it as is
}

/**
 * Calculates the inverse square root of a float number using the fast inverse square root algorithm and then inverts it.
 *
 * @param x The number to calculate the inverse square root of.
 * @return The inverse square root of the input number.
 */
float longSqrt2D(float x) {
    float halfx = 0.5f * x;
    uint32_t i;
    memcpy(&i, &x, sizeof(i)); // Proper type punning in C
    i = 0x5f3759df - (i >> 1); // Magic number for 32-bit floats
    float y;
    memcpy(&y, &i, sizeof(y)); // Convert bits back to float

    y = y * (1.5f - halfx * y * y); // First iteration
    y = y * (1.5f - halfx * y * y); // Second iteration (repeat as needed for precision)

    return 1/y;
}

/**
 * Calculates the sine of a given angle using the Taylor series approximation.
 *
 * @param x The angle in radians.
 * @return The sine of the angle.
 */
float taylor_sin(float x) {
    float result = 0.0f;
    int terms = 20; // More terms for higher precision
    int i;
    for (i = 0; i < terms; ++i) {
        float term = (i % 2 == 0 ? 1 : -1) * pow1(x, 2 * i + 1) / factorial(2 * i + 1);
        result += term;
    }

    return result;
}
/**
 * Calculates the cosine of a given angle using the Taylor series approximation.
 *
 * @param x The angle in radians.
 * @return The cosine of the angle.
 */
float taylor_cos(float x) {
    float result = 0.0f;
    int terms = 20; // More terms for higher precision
    int i;
    for (i = 0; i < terms; ++i) {
        float term = (i % 2 == 0 ? 1 : -1) * pow1(x, 2 * i) / factorial(2 * i);
        result += term;
    }

    return result;
}
/**
 * Initializes the sine and cosine lookup tables.
 * The tables store pre-calculated values for angles in the range [0, 359].
 */
void initialize_tables() {
    int i;
    for (i = 0; i < TABLE_SIZE; ++i) {
        float radians = (float)i * PI / 180.0f;
        sine_table[i] = taylor_sin(radians);
        cosine_table[i] = taylor_cos(radians);
    }
}
/**
 * Approximates the sine of a given angle using the pre-calculated lookup table.
 *
 * @param x The angle in radians.
 * @return The approximate sine of the angle.
 */
float sin_approx(float x) {
    // Convert x to an index in the range [0, 359]
    int index = (int)(x * 180.0f / PI) % TABLE_SIZE;
    return sine_table[index];
}

/**
 * Approximates the cosine of a given angle using the pre-calculated lookup table.
 *
 * @param x The angle in radians.
 * @return The approximate cosine of the angle.
 */
float cos_approx(float x) {
    // Convert x to an index in the range [0, 359]
    int index = (int)(x * 180.0f / PI) % TABLE_SIZE;
    return cosine_table[index];
}
/**
 * Approximates the tangent of a given angle using the pre-calculated lookup table.
 * Returns NaN (Not a Number) when the tangent is undefined.
 *
 * @param x The angle in radians.
 * @return The approximate tangent of the angle.
 */
float tan_approx(float x) {
    float cosine = cos_approx(x);
    // Check if cosine is very close to zero to avoid division by zero
    if (my_fabs(cosine) < 1e-6) {
        return (float)1.0; // Return NaN (Not a Number) when tangent is undefined
    }
    return sin_approx(x) / cosine;
}

/**
 * Multiplies a vector by a 4x4 matrix.
 *
 * @param i The input vector.
 * @param o The output vector.
 * @param m The 4x4 matrix.
 */
void Matrix_MultiplyVector(struct vec3d* i, struct vec3d* o, struct mat4x4* m ){
    float ox = i->x * m->m[0][0] + i->y * m->m[1][0] + i->z * m->m[2][0] + m->m[3][0];
    float oy = i->x * m->m[0][1] + i->y * m->m[1][1] + i->z * m->m[2][1] + m->m[3][1];
    float oz = i->x * m->m[0][2] + i->y * m->m[1][2] + i->z * m->m[2][2] + m->m[3][2];
    float w = i->x * m->m[0][3] + i->y * m->m[1][3] + i->z * m->m[2][3] + m->m[3][3];
    
    o->x = ox;
    o->y = oy;
    o->z = oz;
    o->w = w;
}

/**
 * Initializes an array of triangle pointers with the given triangle array.
 *
 * @param triangles The array of triangles.
 * @param n The number of triangles in the array.
 * @param triArr The array of triangle pointers.
 */
void meshInit(struct triangle triangles[], int n, struct triangle *triArr[])
{
    int i;

    for (i = 0; i < n; ++i)
    {
        triArr[i] = &triangles[i];
    }
}
/**
 * @brief Makes the given matrix an identity matrix.
 * 
 * @param matrix Pointer to the matrix structure.
 */
void Matrix_MakeIdentity(struct mat4x4* matrix)
{
    matrix->m[0][0] = 1.0f;
    matrix->m[1][1] = 1.0f;
    matrix->m[2][2] = 1.0f;
    matrix->m[3][3] = 1.0f;
}

/**
 * @brief Makes the given matrix a rotation matrix around the X-axis.
 * 
 * @param matrix Pointer to the matrix structure.
 * @param fAngleRad The rotation angle in radians.
 */
void Matrix_MakeRotationX(struct mat4x4* matrix, float fAngleRad)
{
    matrix->m[0][0] = 1.0f;
    matrix->m[1][1] = cos_approx(fAngleRad);
    matrix->m[1][2] = sin_approx(fAngleRad);
    matrix->m[2][1] = -sin_approx(fAngleRad);
    matrix->m[2][2] = cos_approx(fAngleRad);
    matrix->m[3][3] = 1.0f;
}

/**
 * @brief Makes the given matrix a rotation matrix around the Y-axis.
 * 
 * @param matrix Pointer to the matrix structure.
 * @param fAngleRad The rotation angle in radians.
 */
void Matrix_MakeRotationY(struct mat4x4* matrix, float fAngleRad)
{
    matrix->m[0][0] = cos_approx(fAngleRad);
    matrix->m[0][2] = sin_approx(fAngleRad);
    matrix->m[2][0] = -sin_approx(fAngleRad);
    matrix->m[1][1] = 1.0f;
    matrix->m[2][2] = cos_approx(fAngleRad);
    matrix->m[3][3] = 1.0f;
}

/**
 * @brief Makes the given matrix a rotation matrix around the Z-axis.
 * 
 * @param matrix Pointer to the matrix structure.
 * @param fAngleRad The rotation angle in radians.
 */
void Matrix_MakeRotationZ(struct mat4x4* matrix, float fAngleRad)
{
    matrix->m[0][0] = cos_approx(fAngleRad);
    matrix->m[0][1] = sin_approx(fAngleRad);
    matrix->m[1][0] = -sin_approx(fAngleRad);
    matrix->m[1][1] = cos_approx(fAngleRad);
    matrix->m[2][2] = 1.0f;
    matrix->m[3][3] = 1.0f;
}
/**
 * Makes a translation matrix.
 *
 * @param matrix The matrix to be modified.
 * @param x The translation along the x-axis.
 * @param y The translation along the y-axis.
 * @param z The translation along the z-axis.
 */
void Matrix_MakeTranslation(struct mat4x4* matrix, float x, float y, float z)
{
    matrix->m[0][0] = 1.0f;
    matrix->m[1][1] = 1.0f;
    matrix->m[2][2] = 1.0f;
    matrix->m[3][3] = 1.0f;
    matrix->m[3][0] = x;
    matrix->m[3][1] = y;
    matrix->m[3][2] = z;
}
/**
 * Makes a projection matrix.
 *
 * @param matrix The matrix to be modified.
 * @param fFovDegrees The field of view in degrees.
 * @param fAspectRatio The aspect ratio of the projection.
 * @param fNear The near clipping plane.
 * @param fFar The far clipping plane.
 */
void Matrix_MakeProjection(struct mat4x4* matrix, float fFovDegrees, float fAspectRatio, float fNear, float fFar)
{
    float fFovRad = 1.0f / tan_approx(fFovDegrees * 0.5f / 180.0f * 3.14159f);
    matrix->m[0][0] = fAspectRatio * fFovRad;
    matrix->m[1][1] = fFovRad;
    matrix->m[2][2] = fFar / (fFar - fNear);
    matrix->m[3][2] = (-fFar * fNear) / (fFar - fNear);
    matrix->m[2][3] = 1.0f;
    matrix->m[3][3] = 0.0f;
}
/**
 * Multiplies two matrices.
 *
 * @param result The resulting matrix.
 * @param m1 The first matrix.
 * @param m2 The second matrix.
 */
void Matrix_MultiplyMatrix(struct mat4x4* result, const struct mat4x4* m1, const struct mat4x4* m2) {
    int c, r;
    for ( c = 0; c < 4; c++)
        for ( r = 0; r < 4; r++)
            result->m[r][c] = m1->m[r][0] * m2->m[0][c] + m1->m[r][1] * m2->m[1][c] + m1->m[r][2] * m2->m[2][c] + m1->m[r][3] * m2->m[3][c];
}

/**
 * Calculates the quick inverse of a matrix.
 *
 * @param result The resulting inverse matrix.
 * @param m The original matrix.
 */
void Matrix_QuickInverse(struct mat4x4* result, const struct mat4x4* m) {
    result->m[0][0] = m->m[0][0]; result->m[0][1] = m->m[1][0]; result->m[0][2] = m->m[2][0]; result->m[0][3] = 0.0f;
    result->m[1][0] = m->m[0][1]; result->m[1][1] = m->m[1][1]; result->m[1][2] = m->m[2][1]; result->m[1][3] = 0.0f;
    result->m[2][0] = m->m[0][2]; result->m[2][1] = m->m[1][2]; result->m[2][2] = m->m[2][2]; result->m[2][3] = 0.0f;
    result->m[3][0] = -(m->m[3][0] * result->m[0][0] + m->m[3][1] * result->m[1][0] + m->m[3][2] * result->m[2][0]);
    result->m[3][1] = -(m->m[3][0] * result->m[0][1] + m->m[3][1] * result->m[1][1] + m->m[3][2] * result->m[2][1]);
    result->m[3][2] = -(m->m[3][0] * result->m[0][2] + m->m[3][1] * result->m[1][2] + m->m[3][2] * result->m[2][2]);
    result->m[3][3] = 1.0f;
}
/**
 * Adds two vectors and stores the result in the given result vector.
 *
 * @param result The vector to store the result in.
 * @param v1 The first vector.
 * @param v2 The second vector.
 */
void Vector_Add(struct vec3d* result, struct vec3d* v1, struct vec3d* v2) {
    result->x = v1->x + v2->x;
    result->y = v1->y + v2->y;
    result->z = v1->z + v2->z;
}

/**
 * Subtracts the second vector from the first vector and stores the result in the given result vector.
 *
 * @param result The vector to store the result in.
 * @param v1 The first vector.
 * @param v2 The second vector.
 */
void Vector_Sub(struct vec3d* result, struct vec3d* v1, struct vec3d* v2) {
    result->x = v1->x - v2->x;
    result->y = v1->y - v2->y;
    result->z = v1->z - v2->z;
}

/**
 * Multiplies the given vector by a scalar value and stores the result in the given result vector.
 *
 * @param result The vector to store the result in.
 * @param v1 The vector to multiply.
 * @param k The scalar value to multiply by.
 */
void Vector_Mul(struct vec3d* result, struct vec3d* v1, float k) {
    result->x = v1->x * k;
    result->y = v1->y * k;
    result->z = v1->z * k;
}

/**
 * Divides the given vector by a scalar value and stores the result in the given result vector.
 *
 * @param result The vector to store the result in.
 * @param v1 The vector to divide.
 * @param k The scalar value to divide by.
 */
void Vector_Div(struct vec3d* result,  struct vec3d* v1, float k) {
    result->x = v1->x / k;
    result->y = v1->y / k;
    result->z = v1->z / k;
}

/**
 * Calculates the dot product of two vectors.
 *
 * @param v1 The first vector.
 * @param v2 The second vector.
 * @return The dot product of the two vectors.
 */
float Vector_DotProduct(struct vec3d* v1, struct vec3d* v2) {
    return v1->x * v2->x + v1->y * v2->y + v1->z * v2->z;
}

/**
 * Calculates the length of the given vector.
 *
 * @param v The vector.
 * @return The length of the vector.
 */
float Vector_Length(struct vec3d* v) {
    return longSqrt2D(Vector_DotProduct(v, v));
}
/**
 * Normalizes a 3D vector.
 *
 * @param result The resulting normalized vector.
 * @param v The vector to be normalized.
 */
void Vector_Normalise(struct vec3d* result, struct vec3d* v) {
    float l = Vector_Length(v);
    result->x = v->x / l;
    result->y = v->y / l;
    result->z = v->z / l;
}

/**
 * Calculates the cross product of two 3D vectors.
 *
 * @param result The resulting cross product vector.
 * @param v1 The first vector.
 * @param v2 The second vector.
 */
void Vector_CrossProduct(struct vec3d* result, const struct vec3d* v1, const struct vec3d* v2) {
    result->x = v1->y * v2->z - v1->z * v2->y;
    result->y = v1->z * v2->x - v1->x * v2->z;
    result->z = v1->x * v2->y - v1->y * v2->x;
}

/**
 * Calculates the transformation matrix that points from a given position to a target,
 * with a specified up direction.
 *
 * @param result The resulting transformation matrix.
 * @param pos The position vector.
 * @param target The target vector.
 * @param up The up direction vector.
 */
void Matrix_PointAt(struct mat4x4* result, struct vec3d* pos,  struct vec3d* target, struct vec3d* up) {
    struct vec3d newForward = {0.0f,0.0f,0.0f,1.0f}, newUp = {0.0f,0.0f,0.0f,1.0f}, newRight = {0.0f,0.0f,0.0f,1.0f};
    //calc forward
    Vector_Sub(&newForward, target, pos);
    Vector_Normalise(&newForward, &newForward);
    struct vec3d a;
    InitializeVector(&a, 0.0f, 0.0f, 0.0f, 1.0f);
    Vector_Mul(&a, &newForward, Vector_DotProduct(up, &newForward) );

    Vector_Sub(&newUp, up, &a);
    Vector_Normalise(&newUp, &newUp);

    Vector_CrossProduct(&newRight, &newUp, &newForward );

    result->m[0][0] = newRight.x;  result->m[0][1] = newRight.y;  result->m[0][2] = newRight.z;  result->m[0][3] = 0.0f;
    result->m[1][0] = newUp.x;     result->m[1][1] = newUp.y;     result->m[1][2] = newUp.z;     result->m[1][3] = 0.0f;
    result->m[2][0] = newForward.x; result->m[2][1] = newForward.y; result->m[2][2] = newForward.z; result->m[2][3] = 0.0f;
    result->m[3][0] = pos->x;       result->m[3][1] = pos->y;       result->m[3][2] = pos->z;       result->m[3][3] = 1.0f;
}

/**
 * Calculates the intersection point between a plane and a line segment.
 *
 * @param result Pointer to the resulting intersection point.
 * @param plane_p Pointer to the point on the plane.
 * @param plane_n Pointer to the normal vector of the plane.
 * @param lineStart Pointer to the starting point of the line segment.
 * @param lineEnd Pointer to the ending point of the line segment.
 */
void Vector_IntersectPlane(struct vec3d* result, struct vec3d* plane_p, struct vec3d* plane_n, struct vec3d* lineStart, struct vec3d* lineEnd) {
    //printf("lineEndReplace: Result (%f, %f, %f)\n", lineEndReplace.x, lineEndReplace.y, lineEndReplace.z);      
    Vector_Normalise(plane_n, plane_n);
    //printf("Vector_Normalise: Result (%f, %f, %f)\n", plane_p_replacement.x, plane_p_replacement.y, plane_p_replacement.z);
    float plane_d = -Vector_DotProduct(plane_n, plane_p);
    //printf("Vector_DotProduct1: Result = %f\n", plane_d);
    float ad = Vector_DotProduct(lineStart, plane_n);
    //printf("Vector_DotProduct2: Result = %f\n", ad);
    float bd = Vector_DotProduct(lineEnd, plane_n);
    //printf("Vector_DotProduct: Result = %f\n", bd);
    float t;
    // if((bd-ad)<1e-6 && (bd-ad)>-1e-6){
    //     t=-plane_d - ad;
    // }
    // else{
    t = (-plane_d - ad) / (bd - ad);//DELAR MED 0 FÖR FAN
    // }
    //printf("t: Result = %f\n", t);
    struct vec3d lineStartToEnd, lineToIntersect; 
    InitializeVector(&lineStartToEnd, 0.0f, 0.0f, 0.0f, 1.0f);
    InitializeVector(&lineToIntersect, 0.0f, 0.0f, 0.0f, 1.0f);

    Vector_Sub(&lineStartToEnd, lineEnd, lineStart);
    Vector_Mul(&lineToIntersect, &lineStartToEnd, t);
    //printf("lineStartToEnd: Result (%f, %f, %f)\n", lineStartToEnd.x, lineStartToEnd.y, lineStartToEnd.z);
    //printf("lineToIntersect: Result (%f, %f, %f)\n", lineToIntersect.x, lineToIntersect.y, lineToIntersect.z);

    struct vec3d resultTemp;
    InitializeVector(&resultTemp, 0.0f, 0.0f, 0.0f, 1.0f);

    Vector_Add(&resultTemp, lineStart, &lineToIntersect);
    result->x = resultTemp.x;
    result->y = resultTemp.y;
    result->z = resultTemp.z;
}

/**
 * Calculates the distance from a point to a plane.
 *
 * @param p         The point for which the distance is calculated.
 * @param plane_p   A point on the plane.
 * @param plane_n   The normal vector of the plane.
 * @return          The distance from the point to the plane.
 */
float DistanceToPlane( struct vec3d* p, struct vec3d* plane_p, struct vec3d* plane_n) {
    struct vec3d normalized_p;
    InitializeVector(&normalized_p, 0.0f, 0.0f, 0.0f, 1.0f);
    Vector_Normalise(&normalized_p, p);
    return (Vector_DotProduct(plane_n, p) - Vector_DotProduct(plane_n, plane_p));
}
/**
 * Clips a triangle against a plane.
 *
 * @param plane_p The point on the plane.
 * @param plane_n The normal vector of the plane.
 * @param in_tri The input triangle to be clipped.
 * @param out_tri1 The first output triangle after clipping.
 * @param out_tri2 The second output triangle after clipping.
 * @return The number of valid output triangles (0, 1, or 2).
 */
int Triangle_ClipAgainstPlane(struct vec3d* plane_p,  struct vec3d* plane_n,  struct triangle* in_tri, struct triangle* out_tri1, struct triangle* out_tri2) {

    Vector_Normalise(plane_n, plane_n);

  
    float d0 = DistanceToPlane(&in_tri->p[0], plane_p, plane_n);
    float d1 = DistanceToPlane(&in_tri->p[1], plane_p, plane_n);
    float d2 = DistanceToPlane(&in_tri->p[2], plane_p, plane_n);

 
    struct vec3d* inside_points[3]= {
        &(struct vec3d){0.0f, 0.0f, 0.0f, 1.0f},
        &(struct vec3d){0.0f, 0.0f, 0.0f, 1.0f},
        &(struct vec3d){0.0f, 0.0f, 0.0f, 1.0f}
    };  int nInsidePointCount = 0;
    //Initilize all the vectors


     struct vec3d* outside_points[3] = {
        &(struct vec3d){0.0f, 0.0f, 0.0f, 1.0f},
        &(struct vec3d){0.0f, 0.0f, 0.0f, 1.0f},
        &(struct vec3d){0.0f, 0.0f, 0.0f, 1.0f}
    };
        int nOutsidePointCount = 0;


    if (d0 >= 0) { inside_points[nInsidePointCount++] = &in_tri->p[0]; }
    else { outside_points[nOutsidePointCount++] = &in_tri->p[0]; }
    if (d1 >= 0) { inside_points[nInsidePointCount++] = &in_tri->p[1]; }
    else { outside_points[nOutsidePointCount++] = &in_tri->p[1]; }
    if (d2 >= 0) { inside_points[nInsidePointCount++] = &in_tri->p[2]; }
    else { outside_points[nOutsidePointCount++] = &in_tri->p[2]; }



    if (nInsidePointCount == 0)
    {

        return 0; // No returned triangles are valid
    }

    if (nInsidePointCount == 3)
    {
   
        //out_tri1 = in_tri; 
        out_tri1->p[0] = in_tri->p[0];
        out_tri1->p[1] = in_tri->p[1];
        out_tri1->p[2] = in_tri->p[2];
        
        return 1; 
    }

    if (nInsidePointCount == 1 && nOutsidePointCount == 2)
    {

        out_tri1->p[0] = *inside_points[0];

        Vector_IntersectPlane(&out_tri1->p[1], plane_p, plane_n, inside_points[0], outside_points[0]);
        Vector_IntersectPlane(&out_tri1->p[2], plane_p, plane_n, inside_points[0], outside_points[1]);

        return 1; // Return the newly formed single triangle
    }

    if (nInsidePointCount == 2 && nOutsidePointCount == 1)
    {

        out_tri1->p[0] = *inside_points[0];
        out_tri1->p[1] = *inside_points[1];
        Vector_IntersectPlane(&out_tri1->p[2], plane_p, plane_n, inside_points[0], outside_points[0]);


        out_tri2->p[0] = *inside_points[1];
        out_tri2->p[1] = out_tri1->p[2];
        Vector_IntersectPlane(&out_tri2->p[2], plane_p, plane_n, inside_points[1], outside_points[0]);

        return 2; 
    }

    return 0; 
}
