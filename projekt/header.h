/* Declare display-related functions from mipslabfunc.c */
void display_image(int x, const uint8_t *data);
void display_init(void);
void display_string(int line, char *s);
void display_update(void);
uint8_t spi_send_recv(uint8_t data);

/* Declare lab-related functions from mipslabfunc.c */
char * itoaconv( int num );
void newlabwork(void);
int nextprime( int inval );
void quicksleep(int cyc);
void tick( unsigned int * timep );


/* Declare display_debug - a function to help debugging.

   After calling display_debug,
   the two middle lines of the display show
   an address and its current contents.

   There's one parameter: the address to read and display.

   Note: When you use this function, you should comment out any
   repeated calls to display_image; display_image overwrites
   about half of the digits shown by display_debug.
*/
void display_debug( volatile int * const addr );

/* Declare bitmap array containing font */
extern const uint8_t const font[128*8];
/* Declare bitmap array containing icon */
extern const uint8_t const icon[128];
/* Declare text buffer for display output */
extern char textbuffer[4][16];

/* Declare functions written by students.
   Note: Since we declare these functions here,
   students must define their functions with the exact types
   specified in the laboratory instructions. */
/* Written as part of asm lab: delay, time2string */
void delay(int);
void time2string( char *, int );
/* Written as part of i/o lab: getbtns, getsw, enable_interrupt */
int getbtns(void);
int getsw(void);
void inputInit(void);
void input(float *fSigma, float *fTheta, float *tx, float *ty, float *tz);
void enable_interrupt(void);




//GEMENSAMMA FUNKTIONER:
#define GRID_WIDTH 128
#define GRID_HEIGHT 32
#define MAX_TRIANGLES 90
#define TABLE_SIZE 720  // One entry per degree
#define PI 3.14159265358979323846
#define GRID_WIDTH 128
#define GRID_HEIGHT 32
#define MAX_TRIANGLES 90
float sine_table[TABLE_SIZE];
float cosine_table[TABLE_SIZE];

int maxTriangles;

struct vec3d
{
    float x, y, z, w;
};


struct triangle
{
    struct vec3d p[3];
};

struct triangle triangles[12];

struct mesh
{
    unsigned int size;
    struct triangle p[12];
};

struct mat4x4
{
    float m[4][4];
};

int abs(int x);
void InitializeMatrix(struct mat4x4 *matrix, float initialValue);
float factorial(int n);
double pow1(double base, int exponent);
double my_fabs(double x);
float longInvSqrt2D(float x);
float taylor_sin(float x);
float taylor_cos(float x);
void initialize_tables();
float sin_approx(float x);
float cos_approx(float x);
float tan_approx(float x);
void MultiplyMatrixVector(struct vec3d *i, struct vec3d *o, struct mat4x4 *m);
void meshInit(struct triangle triangles[], int n, struct triangle *triArr[]);
void Matrix_PointAt(struct mat4x4* result, struct vec3d* pos,  struct vec3d* target, struct vec3d* up);
void Matrix_MakeIdentity(struct mat4x4* matrix);
void Matrix_MakeRotationX(struct mat4x4* matrix, float fAngleRad);
void Matrix_MakeRotationY(struct mat4x4* matrix, float fAngleRad);
void Matrix_MakeRotationZ(struct mat4x4* matrix, float fAngleRad);
void Matrix_MakeTranslation(struct mat4x4* matrix, float x, float y, float z);
void Matrix_MakeProjection(struct mat4x4* matrix, float fFovDegrees, float fAspectRatio, float fNear, float fFar);
void Matrix_MultiplyMatrix(struct mat4x4* result, const struct mat4x4* m1, const struct mat4x4* m2);
void Matrix_QuickInverse(struct mat4x4* result, const struct mat4x4* m);
void Vector_Add(struct vec3d* result, struct vec3d* v1, struct vec3d* v2);
void Vector_Sub(struct vec3d* result, struct vec3d* v1,  struct vec3d* v2);
void Vector_Mul(struct vec3d* result,  struct vec3d* v1, float k);
void Vector_Div(struct vec3d* result, struct vec3d* v1, float k);
float Vector_DotProduct( struct vec3d* v1, struct vec3d* v2);
float Vector_Length( struct vec3d* v);
void Vector_Normalise(struct vec3d* result, struct vec3d* v);
void Vector_CrossProduct(struct vec3d* result, const struct vec3d* v1, const struct vec3d* v2);
void Vector_IntersectPlane(struct vec3d* result, struct vec3d* plane_p, struct vec3d* plane_n, struct vec3d* lineStart, struct vec3d* lineEnd);
float DistanceToPlane( struct vec3d* p, struct vec3d* plane_p, struct vec3d* plane_n);
int Triangle_ClipAgainstPlane(struct vec3d* plane_p,  struct vec3d* plane_n,  struct triangle* in_tri, struct triangle* out_tri1, struct triangle* out_tri2);
void Matrix_MultiplyVector(struct vec3d* i, struct vec3d* o, struct mat4x4* m );
void MoveCamera(float fElapsedTime);
void SetWorldMatrix(struct mat4x4* matWorld, float fElapsedTime);
void SetViewMatrix(struct mat4x4* matView, const struct vec3d* vCamera, const struct vec3d* vLookDir, float fYaw);
void InitializeVector(struct vec3d* vec, float x, float y, float z, float w);
void labinit(int choice);

float fTheta;
float fYaw;
struct vec3d vLookDir;
struct vec3d vCamera;
struct mat4x4 matRotZ, matRotX;
struct mat4x4 matTrans;
float tx;
float ty;
float tz;

float counterTime2;
float timeInSeconds;
float timebefore;
float elapsedTime;
struct mat4x4 matProj;