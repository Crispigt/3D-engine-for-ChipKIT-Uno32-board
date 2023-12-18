# 3D-engine-for-ChipKIT-Uno32-board

## Introduction
This project is a 3D graphics engine tailored for the ChipKIT Uno32 microcontroller board, equipped with the basic ChipKIT I/O shield. It was developed for the "IS1500 Datorteknik och komponenter" course project. Note: The program requires the MCB32 toolchain for compilation.

## Functions
The engine can translate objects from 3D space into Cartesian coordinates. It allows for movement and rotation relative to a cube in 3D space, providing a dynamic and interactive experience.



## Rendering Pipeline
We first have hardcoded objects such as a cube. 

- **World Transformation**:
  - Each triangle of the 3D model is transformed from its local space into world space. This is achieved using a world transformation matrix, which encompasses rotations, scaling, and translation of the model.
  - Functions like `Matrix_MakeRotationX`, `Matrix_MakeRotationY`, `Matrix_MakeRotationZ`, and `Matrix_MakeTranslation` are utilized to create the necessary matrices for these transformations.

- **Back-face Culling**:
  - The engine performs back-face culling to determine which triangles are not visible because they face away from the camera.
  - This is done by calculating the dot product of the triangle's normal vector with the vector from the camera to the triangle. Triangles with a dot product less than zero are facing the camera and are rendered; others are skipped.
  - This step enhances performance by reducing the number of triangles processed further.

- **View Transformation**:
  - Triangles are transformed from world space to view space (camera space) using the view matrix, which is based on the camera's position and orientation.
  - Functions like `Matrix_PointAt` and `Matrix_QuickInverse` are used in this process.
  - This involves moving and rotating the world so the camera is at the origin, facing down the negative Z-axis.

- **Projection Transformation**:
  - Triangles are projected onto a 2D plane (the screen) using a projection matrix created by `Matrix_MakeProjection`.
  - This matrix transforms view space coordinates into clip space coordinates and considers the field of view, aspect ratio, and near and far clipping planes.

- **Clipping**:
  - Triangles partially or fully outside the view frustum are clipped.
  - Parts of triangles outside the near and far clipping planes or the screen's boundaries are discarded or split into smaller triangles within the screen.

- **Screen Transformation**:
  - The clipped triangles are transformed to screen coordinates, scaling and translating the coordinates to map the screen's width and height with the center as the origin.

- **Rasterization**:
  - Triangles are converted into screen pixels in this final step.
  - This is done using `drawTriangle`, which employ Bresenham's line-drawing algorithm to draw triangle edges.

- **Display**:
  - The pixel data is sent to the display buffer and rendered on the screen.

Throughout these steps, the pipeline involves various mathematical and geometric calculations, including matrix multiplications, vector operations, and trigonometric approximations. This pipeline facilitates the creation and display of a dynamic 3D scene from a virtual camera's perspective on a 2D screen.

## Background 
Developed as a school assignment for the Computer Hardware Engineering course at the Royal Institute of Technology (KTH) in Stockholm, this project aims to explore the capabilities of microcontroller-based graphics rendering.

## Requirements
- ChipKIT Uno32 development board with ChipKIT Basic I/O shield.
- MCB32 toolchain for building and deploying the code. [Download here](https://github.com/is1200-example-projects/mcb32tools/releases/).

## Inspiration
The project was inspired by Onelonecoder's educational videos on building a 3D engine, which provided valuable insights into the fundamentals of 3D graphics programming.

## Copying
Certain files such as `stubs.c`, `vectors.S`, `func.c`, and parts of `main.c` are reused from a course lab.
