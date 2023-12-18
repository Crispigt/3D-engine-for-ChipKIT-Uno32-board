#include <stdint.h>  /* Declarations of uint_32 and the like */
#include <stdlib.h>
#include <stdio.h>
#include <pic32mx.h> /* Declarations of system-specific addresses etc */
#include "header.h" /* Declatations for these labs */

/**
 * @brief Retrieves the value of the switches connected to PORTD.
 * 
 * This function reads the value of the switches connected to PORTD,
 * applies a bitmask to isolate the relevant bits, and returns the
 * result after shifting it 8 places to the right.
 * 
 * @return The value of the switches.
 */
int getsw( void )
{
    unsigned int temp1;
    unsigned int temp2;
    unsigned int loadPortD;
    loadPortD = PORTD;
    //bitmask so we get ones everywhere and shift 8 places 
    unsigned int bitmask2 = 0x0F00;
    temp1 = (loadPortD & bitmask2);
    temp2 = (temp1 >> 8);
    return temp2;
}
/**
 * Retrieves the state of the buttons.
 * 
 * @return The state of the buttons as an integer.
 */
int getbtns(void)
{
    unsigned int loadPortD = PORTD;
    unsigned int bitmask = 0x00E0; // Modified bitmask to include button 1
    unsigned int temp = ((loadPortD & bitmask) >> 4) | ((PORTF >> 1) & 0x1);
    return temp;
}
/**
 * @brief Initializes the input settings.
 * 
 * This function sets up the necessary configurations for input handling.
 * It configures the timer, clears flags, and initializes the input ports.
 * 
 * @note This function assumes that the necessary registers and variables are defined.
 */

void inputInit(void){
    T2CON = 0; // Clear time control reg
    T2CON |= (0b111 << 4);
    TMR2 = 0; //clear timer register
    PR2 = 31250; //16 bit timer
    T2CONSET = 0x8000; // start timer
    IFSCLR(0) = (1<<8); //clear flag
    //INPUT---------------------------------------
    unsigned int clearedValueOfPORTE;
    unsigned int valueInPORTE;
    volatile unsigned int *p;
    //Point to PORTE, load value to variabel, bitmask out bit 7-0
    p = (volatile unsigned int *) 0xbf886100;
    valueInPORTE = *p;
    unsigned int bitMask = 0xFFFFFF00;
    clearedValueOfPORTE = valueInPORTE & bitMask;
    //Load new value to PORTE
    *p = clearedValueOfPORTE;

     TRISF |= 0x2; // set bit 1 to input (1) BTN1

    unsigned int tempPortD;
    unsigned int clearedValuePortD;
    tempPortD = TRISD;
    unsigned int bitmask2 = 0x00000FF0;
    clearedValuePortD = tempPortD | bitmask2;
    TRISD = clearedValuePortD;
    //Could have been done in one line TRISD |= 0x0FE0
    // enable_interrupt();
}



/**
 * Sets the world matrix based on the given elapsed time.
 *
 * @param matWorld      Pointer to the world matrix.
 * @param fElapsedTime The elapsed time.
 */
void SetWorldMatrix(struct mat4x4* matWorld, float fElapsedTime) {
    // Set up "World Transform" though not updating theta
    // makes this a bit redundant
    struct mat4x4 matRotZ, matRotX;
    InitializeMatrix(&matRotZ, 0);
    InitializeMatrix(&matRotX, 0);
    //fTheta += 1.0f * fElapsedTime; 
    Matrix_MakeRotationZ(&matRotZ, fTheta * 0.5f);
    // printmatrices(&matRotZ, "matRotZ");
    Matrix_MakeRotationX(&matRotX, fTheta);
    // printmatrices(&matRotX, "matRotX");

    struct mat4x4 matTrans;
    InitializeMatrix(&matTrans, 0);
    Matrix_MakeTranslation(&matTrans, 0.0f, 0.0f, 5.0f);
    // printmatrices(&matTrans, "matTrans");

    Matrix_MakeIdentity(matWorld);	// Form World Matrix
    // printmatrices(matWorld, "matWorld12");
    Matrix_MultiplyMatrix(matWorld, &matRotZ, &matRotX); // Transform by rotation
    // printmatrices(matWorld, "matWorld13");
    Matrix_MultiplyMatrix(matWorld, matWorld, &matTrans); // Transform by translation
    // printmatrices(matWorld, "matWorld14");
}


/**
 * Sets the view matrix based on the camera position, look direction, and yaw angle.
 *
 * @param matView   Pointer to the view matrix to be set.
 * @param vCamera   Pointer to the camera position vector.
 * @param vLookDir  Pointer to the look direction vector.
 * @param fYaw      Yaw angle in radians.
 */

void SetViewMatrix(struct mat4x4* matView, const struct vec3d* vCamera, const struct vec3d* vLookDir, float fYaw) {
        struct vec3d vUp = { 0,1,0,1 };
        struct vec3d vTarget = { 0,0,1,1 };
        struct vec3d vLookDirCopy;
        vLookDirCopy.x = vLookDir->x;
        vLookDirCopy.y = vLookDir->y;
        vLookDirCopy.z = vLookDir->z;
        struct vec3d vCameraCopy;
        vCameraCopy.x = vLookDir->x;
        vCameraCopy.y = vLookDir->y;
        vCameraCopy.z = vLookDir->z;
        struct mat4x4 matCameraRot;
        InitializeMatrix(&matCameraRot,0);
        Matrix_MakeRotationY(&matCameraRot, fYaw);
        Matrix_MultiplyVector(&vTarget, &vLookDirCopy,&matCameraRot);//man kan inte ändra const

        Vector_Add(&vTarget, &vCameraCopy, &vLookDirCopy);
        struct mat4x4 matCamera;
        InitializeMatrix(&matCamera, 0);
        Matrix_PointAt(&matCamera, &vCameraCopy, &vTarget, &vUp);

        // Make view matrix from camera
        Matrix_QuickInverse(matView, &matCamera);
}

