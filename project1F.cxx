// Access my movie at: https://www.youtube.com/watch?v=tHsB8PdgAMY
#include <iostream>
#include <iomanip>
#include <vtkDataSet.h>
#include <vtkImageData.h>
#include <vtkPNGWriter.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPoints.h>
#include <vtkUnsignedCharArray.h>
#include <vtkFloatArray.h>
#include <vtkCellArray.h>
#include <vtkDoubleArray.h>
#include <tgmath.h>
#include <math.h>
#include <string.h>

#define NORMALS

using std::cerr;
using std::endl;




double ceil__441(double f)
{
    return ceil(f - 0.00001);
}

double floor__441(double f)
{
    return floor(f + 0.00001);
}


// Puts the cross product of a and b into the passed-by-reference array, ret
void crossProduct(double a[3], double b[3], double(&ret)[3]) {
    ret[0] = (a[1] * b[2]) - (a[2] * b[1]);
    ret[1] = (b[0] * a[2]) - (a[0] * b[2]);
    ret[2] = (a[0] * b[1]) - (a[1] * b[0]);
}

// normalizes the passed-in vector per the || i || = sqrt(i.x squared) + (i.y squared) + (i.z squared)) formula
void normalizeVector(double in[3], double(&out)[3]) {
    double heehee = (in[0] * in[0]) + (in[1] * in[1]) + (in[2] * in[2]);
    double in_norm = sqrt(heehee);
    out[0] = in[0] / in_norm;
    out[1] = in[1] / in_norm;
    out[2] = in[2] / in_norm;
}

// returns a double that's the dot product of the two input arrays 
double dotProduct(double a[3], double b[3]) {
    double ret;
    ret = (a[0] * b[0]) + (a[1] * b[1]) + (a[2] * b[2]);
    return ret;
}

void scalarProduct(double scalar, double (&a)[3]) {
    a[0] = a[0] * scalar;
    a[1] = a[1] * scalar;
    a[2] = a[2] * scalar;
}

void vectorSubtract(double(&in)[3], double subtractor[3]) {
    in[0] = in[0] - subtractor[0];
    in[1] = in[1] - subtractor[1];
    in[2] = in[2] - subtractor[2];
}

void maxWithZero(double(&in)) {
    if (0 >= in) {
        in = 0;
    }
}

vtkImageData*
NewImage(int height, int width)
{
    vtkImageData* img = vtkImageData::New();
    img->SetDimensions(width, height, 1);
    img->AllocateScalars(VTK_UNSIGNED_CHAR, 3);

    return img;
}

void
WriteImage(vtkImageData* img, const char* filename)
{
    std::string full_filename = filename;
    full_filename += ".png";
    vtkPNGWriter* writer = vtkPNGWriter::New();
    writer->SetInputData(img);
    writer->SetFileName(full_filename.c_str());
    writer->Write();
    writer->Delete();
}

class Matrix
{
public:
    double          A[4][4];  // A[i][j] means row i, column j

    void            TransformPoint(const double* ptIn, double* ptOut);
    static Matrix   ComposeMatrices(const Matrix&, const Matrix&);
    void            Print(ostream& o);
};

void
Matrix::Print(ostream& o)
{
    for (int i = 0; i < 4; i++)
    {
        char str[256];
        sprintf(str, "(%.7f %.7f %.7f %.7f)\n", A[i][0], A[i][1], A[i][2], A[i][3]);
        o << str;
    }
}

Matrix
Matrix::ComposeMatrices(const Matrix& M1, const Matrix& M2)
{
    Matrix rv;
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
        {
            rv.A[i][j] = 0;
            for (int k = 0; k < 4; k++)
                rv.A[i][j] += M1.A[i][k] * M2.A[k][j];
        }

    return rv;
}

void
Matrix::TransformPoint(const double* ptIn, double* ptOut)
{
    ptOut[0] = ptIn[0] * A[0][0]
        + ptIn[1] * A[1][0]
        + ptIn[2] * A[2][0]
        + ptIn[3] * A[3][0];
    ptOut[1] = ptIn[0] * A[0][1]
        + ptIn[1] * A[1][1]
        + ptIn[2] * A[2][1]
        + ptIn[3] * A[3][1];
    ptOut[2] = ptIn[0] * A[0][2]
        + ptIn[1] * A[1][2]
        + ptIn[2] * A[2][2]
        + ptIn[3] * A[3][2];
    ptOut[3] = ptIn[0] * A[0][3]
        + ptIn[1] * A[1][3]
        + ptIn[2] * A[2][3]
        + ptIn[3] * A[3][3];
}

class Camera
{
public:
    double          near, far;
    double          angle;
    double          position[3];
    double          focus[3];
    double          up[3];
    // My additional members
    double          O[3];
    double          w[3];
    double          v[3];
    double          u[3];
    double          m = 1000;
    double          n = 1000;

    Matrix          ViewTransform();
    void            CameraFrameCalc();
    Matrix          CameraTransform();
    Matrix          DeviceTransform();
};


// Potential bug: I'm using w instead of strictly O - focus to calculate v and u... shouldn't be a problem?
void Camera::CameraFrameCalc() {
    for (int i = 0; i < 3; i++) {
        O[i] = position[i];
        w[i] = O[i] - focus[i];
        v[i] = up[i];
    }
    // Passing v and u by reference, so think of it as 
    //u = up X w 
    //for example
    crossProduct(up, w, u);
    crossProduct(w, u, v);

    // Normalize w, u , v
    double norm_w[3];
    double norm_v[3];
    double norm_u[3];

    // Note that normalizeVector does the i / ||i|| step as well
    normalizeVector(w, norm_w);
    normalizeVector(v, norm_v);
    normalizeVector(u, norm_u);

    // Set w, v, u to themselves divided by their norms
    for (int i = 0; i < 3; i++) {
        w[i] = norm_w[i];
        v[i] = norm_v[i];
        u[i] = norm_u[i];
    }

    

}

// Returns a camera transform matrix based on the Camera class variables O, u, v, w, and the dot product function
Matrix Camera::CameraTransform() {
    Matrix ret;

    ret.A[0][0] = u[0];
    ret.A[0][1] = v[0];
    ret.A[0][2] = w[0];
    ret.A[0][3] = 0;

    ret.A[1][0] = u[1];
    ret.A[1][1] = v[1];
    ret.A[1][2] = w[1];
    ret.A[1][3] = 0;

    ret.A[2][0] = u[2];
    ret.A[2][1] = v[2];
    ret.A[2][2] = w[2];
    ret.A[2][3] = 0;

    double t[3];
    t[0] = 0 - O[0];
    t[1] = 0 - O[1];
    t[2] = 0 - O[2];

    ret.A[3][0] = dotProduct(u, t);
    ret.A[3][1] = dotProduct(v, t);
    ret.A[3][2] = dotProduct(w, t);
    ret.A[3][3] = 1;

    return ret;
}

// Returns a view transformation matrix based on the Camera class variables angle, near, and far
Matrix Camera::ViewTransform() {
    Matrix ret;

    ret.A[0][0] = cos(angle/2)/sin(angle/2);
    ret.A[0][1] = 0;
    ret.A[0][2] = 0;
    ret.A[0][3] = 0;

    ret.A[1][0] = 0;
    ret.A[1][1] = cos(angle/2) / sin(angle/2);
    ret.A[1][2] = 0;
    ret.A[1][3] = 0;

    ret.A[2][0] = 0;
    ret.A[2][1] = 0;
    ret.A[2][2] = (far + near) / (far - near);
    ret.A[2][3] = -1;

    ret.A[3][0] = 0;
    ret.A[3][1] = 0;
    ret.A[3][2] = (2 * near * far) / (far - near);
    ret.A[3][3] = 0;

    return ret;
}

// Returns a device transformation matrix based on the Camera class variables n and m, the screen dimensions
Matrix Camera::DeviceTransform() {
    Matrix ret;

    ret.A[0][0] = n / 2;
    ret.A[0][1] = 0;
    ret.A[0][2] = 0;
    ret.A[0][3] = 0;

    ret.A[1][0] = 0;
    ret.A[1][1] = m / 2;
    ret.A[1][2] = 0;
    ret.A[1][3] = 0;

    ret.A[2][0] = 0;
    ret.A[2][1] = 0;
    ret.A[2][2] = 1;
    ret.A[2][3] = 0;
    
    ret.A[3][0] = n / 2;
    ret.A[3][1] = m / 2;
    ret.A[3][2] = 0;
    ret.A[3][3] = 1;

    return ret;

}

struct LightingParameters
{
    LightingParameters(void)
    {
        lightDir[0] = -0.6;
        lightDir[1] = 0;
        lightDir[2] = -0.8;
        Ka = 0.3;
        Kd = 0.7;
        Ks = 2.8;
        alpha = 50.5;
    };


    double lightDir[3]; // The direction of the light source
    double Ka;          // The coefficient for ambient lighting
    double Kd;          // The coefficient for diffuse lighting
    double Ks;          // The coefficient for specular lighting
    double alpha;       // The exponent term for specular lighting
};

LightingParameters lp;

LightingParameters
GetLighting(Camera c)
{
    LightingParameters lp;
    lp.lightDir[0] = c.position[0] - c.focus[0];
    lp.lightDir[1] = c.position[1] - c.focus[1];
    lp.lightDir[2] = c.position[2] - c.focus[2];
    double mag = sqrt(lp.lightDir[0] * lp.lightDir[0]
        + lp.lightDir[1] * lp.lightDir[1]
        + lp.lightDir[2] * lp.lightDir[2]);
    if (mag > 0)
    {
        lp.lightDir[0] /= mag;
        lp.lightDir[1] /= mag;
        lp.lightDir[2] /= mag;
    }

    return lp;
}

double SineParameterize(int curFrame, int nFrames, int ramp)
{
    int nNonRamp = nFrames - 2 * ramp;
    double height = 1. / (nNonRamp + 4 * ramp / M_PI);
    if (curFrame < ramp)
    {
        double factor = 2 * height * ramp / M_PI;
        double eval = cos(M_PI / 2 * ((double)curFrame) / ramp);
        return (1. - eval) * factor;
    }
    else if (curFrame > nFrames - ramp)
    {
        int amount_left = nFrames - curFrame;
        double factor = 2 * height * ramp / M_PI;
        double eval = cos(M_PI / 2 * ((double)amount_left / ramp));
        return 1. - (1 - eval) * factor;
    }
    double amount_in_quad = ((double)curFrame - ramp);
    double quad_part = amount_in_quad * height;
    double curve_part = height * (2 * ramp) / M_PI;
    return quad_part + curve_part;
}

Camera
GetCamera(int frame, int nframes)
{
    double t = SineParameterize(frame, nframes, nframes / 10);
    Camera c;
    c.near = 5;
    c.far = 200;
    c.angle = M_PI / 6;
    c.position[0] = 40 * sin(2 * M_PI * t);
    c.position[1] = 40 * cos(2 * M_PI * t);
    c.position[2] = 40;
    c.focus[0] = 0;
    c.focus[1] = 0;
    c.focus[2] = 0;
    c.up[0] = 0;
    c.up[1] = 1;
    c.up[2] = 0;
    return c;
}


class Triangle
{
public:
    double         X[3];
    double         Y[3];
    double         Z[3];
    double colors[3][3];
    double normals[3][3];
    double shading[3];

    // Sets an arbitrary set of vertices such that, based on X value, they are re-ordered
    // in their array as Left, Middle, and Right.
    void setVertArb(double(&q)[3], double(&Y)[3], double (&Z)[3], double (&triColor)[3][3], double(&shading)[3]) {

        // Case 1A: MLR - works
        if (q[1] < q[0] && q[0] < q[2]) {
            std::swap(q[0], q[1]);
            std::swap(Y[0], Y[1]);
            std::swap(Z[0], Z[1]);
            std::swap(shading[0], shading[1]);
            //color swap
            for (int i = 0; i < 3; i++) {
                std::swap(triColor[0][i], triColor[1][i]);
            }
        }
        // Case 1B: MRL
        else if (q[1] > q[2] && q[2] < q[0] && q[0] < q[1]) {
            std::swap(q[0], q[2]);
            std::swap(Y[0], Y[2]);
            std::swap(Z[0], Z[2]);
            std::swap(shading[0], shading[2]);
            for (int i = 0; i < 3; i++) {
                std::swap(triColor[0][i], triColor[2][i]);
            }


            std::swap(q[1], q[2]);
            std::swap(Y[1], Y[2]);
            std::swap(Z[1], Z[2]);
            std::swap(shading[1], shading[2]);
            for (int i = 0; i < 3; i++) {
                std::swap(triColor[1][i], triColor[2][i]);
            }



        }
        // Case2A: LRM
        else if (q[0] < q[2] && q[2] < q[1]) {
            std::swap(q[1], q[2]);
            std::swap(Y[1], Y[2]);
            std::swap(Z[1], Z[2]);
            std::swap(shading[1], shading[2]);
            for (int i = 0; i < 3; i++) {
                std::swap(triColor[1][i], triColor[2][i]);
            }

        }
        // Case2B: RLM
        else if (q[2] < q[0] && q[1] < q[2]) {
            std::swap(q[1], q[2]);
            std::swap(Y[1], Y[2]);
            std::swap(Z[1], Z[2]);
            std::swap(shading[1], shading[2]);
            for (int i = 0; i < 3; i++) {
                std::swap(triColor[1][i], triColor[2][i]);
            }

            std::swap(q[0], q[2]);
            std::swap(Y[0], Y[2]);
            std::swap(Z[0], Z[2]);
            std::swap(shading[0], shading[2]);
            for (int i = 0; i < 3; i++) {
                std::swap(triColor[0][i], triColor[2][i]);
            }
        }
        // Case 3A: LMR
        else if (q[0] < q[1] && q[1] < q[2]) {
            // Do nothing, already ordered
        }
        // Case 3B: RML
        else if (q[2] < q[1] && q[1] < q[0]) {
            std::swap(q[0], q[2]);
            std::swap(Y[0], Y[2]);
            std::swap(Z[0], Z[2]);
            std::swap(shading[0], shading[2]);
            for (int i = 0; i < 3; i++) {
                std::swap(triColor[0][i], triColor[2][i]);
            }
        }
        else {
            cout << "This should not appear - setVertArb output for a triangle that's not arb" << endl;
        }
    }

};



class Screen
{
public:
    unsigned char* buffer;
    int width, height;
};


std::vector<Triangle>
GetTriangles(void)
{
    vtkPolyDataReader* rdr = vtkPolyDataReader::New();
    rdr->SetFileName("proj1f_geometry.vtk");
    cerr << "Reading" << endl;
    rdr->Update();
    cerr << "Done reading" << endl;
    if (rdr->GetOutput()->GetNumberOfCells() == 0)
    {
        cerr << "Unable to open file!!" << endl;
        exit(EXIT_FAILURE);
    }
    vtkPolyData* pd = rdr->GetOutput();

    int numTris = pd->GetNumberOfCells();
    vtkPoints* pts = pd->GetPoints();
    vtkCellArray* cells = pd->GetPolys();
    vtkDoubleArray* var = (vtkDoubleArray*)pd->GetPointData()->GetArray("hardyglobal");
    double* color_ptr = var->GetPointer(0);
    //vtkFloatArray *var = (vtkFloatArray *) pd->GetPointData()->GetArray("hardyglobal");
    //float *color_ptr = var->GetPointer(0);
    vtkFloatArray* n = (vtkFloatArray*)pd->GetPointData()->GetNormals();
    float* normals = n->GetPointer(0);
    std::vector<Triangle> tris(numTris);
    vtkIdType npts;
    vtkIdType* ptIds;
    int idx;
    for (idx = 0, cells->InitTraversal(); cells->GetNextCell(npts, ptIds); idx++)
    {
        if (npts != 3)
        {
            cerr << "Non-triangles!! ???" << endl;
            exit(EXIT_FAILURE);
        }
        double* pt = NULL;
        pt = pts->GetPoint(ptIds[0]);
        tris[idx].X[0] = pt[0];
        tris[idx].Y[0] = pt[1];
        tris[idx].Z[0] = pt[2];
#ifdef NORMALS
        tris[idx].normals[0][0] = normals[3 * ptIds[0] + 0];
        tris[idx].normals[0][1] = normals[3 * ptIds[0] + 1];
        tris[idx].normals[0][2] = normals[3 * ptIds[0] + 2];
#endif
        pt = pts->GetPoint(ptIds[1]);
        tris[idx].X[1] = pt[0];
        tris[idx].Y[1] = pt[1];
        tris[idx].Z[1] = pt[2];
#ifdef NORMALS
        tris[idx].normals[1][0] = normals[3 * ptIds[1] + 0];
        tris[idx].normals[1][1] = normals[3 * ptIds[1] + 1];
        tris[idx].normals[1][2] = normals[3 * ptIds[1] + 2];
#endif
        pt = pts->GetPoint(ptIds[2]);
        tris[idx].X[2] = pt[0];
        tris[idx].Y[2] = pt[1];
        tris[idx].Z[2] = pt[2];
#ifdef NORMALS
        tris[idx].normals[2][0] = normals[3 * ptIds[2] + 0];
        tris[idx].normals[2][1] = normals[3 * ptIds[2] + 1];
        tris[idx].normals[2][2] = normals[3 * ptIds[2] + 2];
#endif

        // 1->2 interpolate between light blue, dark blue
        // 2->2.5 interpolate between dark blue, cyan
        // 2.5->3 interpolate between cyan, green
        // 3->3.5 interpolate between green, yellow
        // 3.5->4 interpolate between yellow, orange
        // 4->5 interpolate between orange, brick
        // 5->6 interpolate between brick, salmon
        double mins[7] = { 1, 2, 2.5, 3, 3.5, 4, 5 };
        double maxs[7] = { 2, 2.5, 3, 3.5, 4, 5, 6 };
        unsigned char RGB[8][3] = { { 71, 71, 219 },
                                    { 0, 0, 91 },
                                    { 0, 255, 255 },
                                    { 0, 128, 0 },
                                    { 255, 255, 0 },
                                    { 255, 96, 0 },
                                    { 107, 0, 0 },
                                    { 224, 76, 76 }
        };
        for (int j = 0; j < 3; j++)
        {
            float val = color_ptr[ptIds[j]];
            int r;
            for (r = 0; r < 7; r++)
            {
                if (mins[r] <= val && val < maxs[r])
                    break;
            }
            if (r == 7)
            {
                cerr << "Could not interpolate color for " << val << endl;
                exit(EXIT_FAILURE);
            }
            double proportion = (val - mins[r]) / (maxs[r] - mins[r]);
            tris[idx].colors[j][0] = (RGB[r][0] + proportion * (RGB[r + 1][0] - RGB[r][0])) / 255.0;
            tris[idx].colors[j][1] = (RGB[r][1] + proportion * (RGB[r + 1][1] - RGB[r][1])) / 255.0;
            tris[idx].colors[j][2] = (RGB[r][2] + proportion * (RGB[r + 1][2] - RGB[r][2])) / 255.0;
        }
    }

    return tris;
}


// Returns the Y-value for the endpoint of the vertical line that comes off of M, ie
// (X[1], Y[1])
double endpointY(double X[3], double Y[3]) {
    double midy;
    double m;
    double b;

    m = (Y[2] - Y[0]) / (X[2] - X[0]);
    // 2) Find slope intercept
    b = Y[2] - (m * X[2]);
    // 3) Plug in i to get the value we need
    midy = (m * X[1]) + b;
    return midy;
}


void rasterizeGoingR(Triangle tri, unsigned char* buffer, double**depthBuffer) {

    
    double columnMin = ceil__441(tri.X[1]);
    double columnMax = floor__441(tri.X[2]);


    // For every column between columnMin and columnMax, execute the scanline algorithm
    for (double i = columnMin; i <= columnMax; i++) {

        double bottomEnd;
        double topEnd;
        double m;
        double b;
        double m1;
        double b1;
        int index;


        // For topEnd
        // 1) Find slope
        m = (tri.Y[2] - tri.Y[0]) / (tri.X[2] - tri.X[0]);
        // 2) Find slope intercept
        b = tri.Y[2] - (m * tri.X[2]);
        // 3) Plug in i to get the value we need
        topEnd = (m * i) + b;

        //LERP Z, color, and shading for topEnd
        double t1 = (i - tri.X[2]) / (tri.X[0] - tri.X[2]);

        double topEnd_Z = tri.Z[2] + (t1 * (tri.Z[0] - tri.Z[2]));
        double topEnd_colorR = tri.colors[2][0] + (t1 * (tri.colors[0][0] - tri.colors[2][0]));
        double topEnd_colorG = tri.colors[2][1] + (t1 * (tri.colors[0][1] - tri.colors[2][1]));
        double topEnd_colorB = tri.colors[2][2] + (t1 * (tri.colors[0][2] - tri.colors[2][2]));
        double topEnd_shading = tri.shading[2] + (t1 * (tri.shading[0] - tri.shading[2]));


        // Repeat for bottomEnd
        m1 = (tri.Y[1] - tri.Y[2]) / (tri.X[1] - tri.X[2]);
        b1 = tri.Y[2] - (m1 * tri.X[2]);
        bottomEnd = (m1 * i) + b1;

        //Lerp for bottomEnd
        double t2 = (i - tri.X[1]) / (tri.X[2] - tri.X[1]);

        double bottomEnd_Z = tri.Z[1] + (t2 * (tri.Z[2] - tri.Z[1]));
        double bottomEnd_colorR = tri.colors[1][0] + (t2 * (tri.colors[2][0] - tri.colors[1][0]));
        double bottomEnd_colorG = tri.colors[1][1] + (t2 * (tri.colors[2][1] - tri.colors[1][1]));
        double bottomEnd_colorB = tri.colors[1][2] + (t2 * (tri.colors[2][2] - tri.colors[1][2]));
        double bottomEnd_shading = tri.shading[1] + (t2 * (tri.shading[2] - tri.shading[1]));


        // scanline algorithm starting on the leftmost line
        for (double j = ceil__441(bottomEnd); j <= floor__441(topEnd); j++) {
            //Lerp j
            double tj = (j - topEnd) / (bottomEnd - topEnd);

            double jZ = topEnd_Z + (tj * (bottomEnd_Z - topEnd_Z));
            double jColorR = topEnd_colorR + (tj * (bottomEnd_colorR - topEnd_colorR));
            double jColorG = topEnd_colorG + (tj * (bottomEnd_colorG - topEnd_colorG));
            double jColorB = topEnd_colorB + (tj * (bottomEnd_colorB - topEnd_colorB));
            double jshading = topEnd_shading + (tj * (bottomEnd_shading - topEnd_shading));
           
            

            index = ((j * 1000) + i) * 3;

            if (index >= 0 && index + 2 <= (1000 * 1000 * 3) && floor__441(topEnd) <= 1000 && i < 1000 && i >= 0) {

                if (jZ > depthBuffer[int(i)][int(j)]) {
                    depthBuffer[int(i)][int(j)] = jZ;


                    // changed from unsigned char colorG = ceil__441(jColorG * 255);
                    unsigned char colorR;
                    unsigned char colorG;
                    unsigned char colorB;

                    double jcr = jColorR * jshading;
                    if (jcr < 1) {
                        colorR = ceil__441(255 * jcr);
                    }
                    else {
                        colorR = 255;
                    }

                    double jcg = jColorG * jshading;
                    if (jcg < 1) {
                        colorG = ceil__441(255 * jcg);
                    }
                    else {
                        colorG = 255;
                    }

                    double jcb = jColorB * jshading;
                    if (jcb < 1) {
                        colorB = ceil__441(255 * jcb);
                    }
                    else {
                        colorB = 255;
                    }
                    buffer[index + 0] = colorR;
                    buffer[index + 1] = colorG;
                    buffer[index + 2] = colorB;
                }
            }

        }

    }
}



void rasterizeGoingL(Triangle tri, unsigned char* buffer, double **depthBuffer) {

    int columnMax = floor__441(tri.X[1]);
    int columnMin = ceil__441(tri.X[2]);
    
    // For every column between columnMin and columnMax, execute the scanline algorithm
    for (double i = columnMax; i >= columnMin; i--) {
        

        double bottomEnd;
        double topEnd;
        double m;
        double b;
        double m1;
        double b1;
        int index;

        
        // For topEnd
        // 1) Find slope
        m = (tri.Y[2] - tri.Y[0]) / (tri.X[2] - tri.X[0]);
        // 2) Find slope intercept
        b = tri.Y[2] - (m * tri.X[2]);
        // 3) Plug in i to get the value we need
        topEnd = (m * i) + b;

        double t1 = (i - tri.X[0]) / (tri.X[2] - tri.X[0]);

        double topEnd_Z = tri.Z[0] + (t1 * (tri.Z[2] - tri.Z[0]));
        double topEnd_colorR = tri.colors[0][0] + (t1 * (tri.colors[2][0] - tri.colors[0][0]));
        double topEnd_colorG = tri.colors[0][1] + (t1 * (tri.colors[2][1] - tri.colors[0][1]));
        double topEnd_colorB = tri.colors[0][2] + (t1 * (tri.colors[2][2] - tri.colors[0][2]));
        double topEnd_shading = tri.shading[0] + (t1 * (tri.shading[2] - tri.shading[0]));


        // Repeat for bottomEnd
        m1 = (tri.Y[1] - tri.Y[2]) / (tri.X[1] - tri.X[2]);
        b1 = tri.Y[2] - (m1 * tri.X[2]);
        bottomEnd = (m1 * i) + b1;

        double t2 = (i - tri.X[1]) / (tri.X[2] - tri.X[1]);

        double bottomEnd_Z = tri.Z[1] + (t2 * (tri.Z[2] - tri.Z[1]));
        double bottomEnd_colorR = tri.colors[1][0] + (t2 * (tri.colors[2][0] - tri.colors[1][0]));
        double bottomEnd_colorG = tri.colors[1][1] + (t2 * (tri.colors[2][1] - tri.colors[1][1]));
        double bottomEnd_colorB = tri.colors[1][2] + (t2 * (tri.colors[2][2] - tri.colors[1][2]));
        double bottomEnd_shading = tri.shading[1] + (t1 * (tri.shading[2] - tri.shading[1]));

        // scanline algorithm starting on the leftmost line             
        for (double j = ceil__441(bottomEnd); j <= floor__441(topEnd); j++) {

            double tj = (j - topEnd) / (bottomEnd - topEnd);

            double jZ = topEnd_Z + (tj * (bottomEnd_Z - topEnd_Z));
            double jColorR = topEnd_colorR + (tj * (bottomEnd_colorR - topEnd_colorR));
            double jColorG = topEnd_colorG + (tj * (bottomEnd_colorG - topEnd_colorG));
            double jColorB = topEnd_colorB + (tj * (bottomEnd_colorB - topEnd_colorB));
            double jshading = topEnd_shading + (tj * (bottomEnd_shading - topEnd_shading));

            index = (j * 1000 + i) * 3;

            if (index >= 0 && index + 2 <= (1000 * 1000 * 3) && floor__441(topEnd) <= 1000 && i <1000 && i >=0) {

                if (jZ > depthBuffer[int(i)][int(j)]) {
                    depthBuffer[int(i)][int(j)] = jZ;

                    // changed from unsigned char colorG = ceil__441(jColorG * 255);
                    unsigned char colorR;
                    unsigned char colorG; 
                    unsigned char colorB;

                    double jcr = jColorR * jshading;
                    if (jcr < 1) {
                        colorR = ceil__441(255 * jcr);
                    }
                    else {
                        colorR = 255;
                    }

                    double jcg = jColorG * jshading;
                    if (jcg < 1) {
                        colorG = ceil__441(255 * jcg);
                    }
                    else {
                        colorG = 255;
                    }

                    double jcb = jColorB * jshading;
                    if (jcb < 1) {
                        colorB = ceil__441(255 * jcb);
                    }
                    else {
                        colorB = 255;
                    }
                    

                    buffer[index + 0] = colorR;
                    buffer[index + 1] = colorG;
                    buffer[index + 2] = colorB;
                }
                
            }

        }

    }

    
}



bool isArb(double q[3], double Y[3]) {
    // Case 1A: MLR - works
    if (q[1] < q[0] && q[0] < q[2]) {
        return true;
    }
    // Case 1B: MRL
    else if (q[1] > q[2] && q[2] < q[0] && q[0] < q[1]) {
        return true;
    }
    // Case2A: LRM
    else if (q[0] < q[2] && q[2] < q[1]) {
        return true;
    }
    // Case2B: RLM
    else if (q[2] < q[0] && q[1] < q[2]) {
        return true;
    }
    // Case 3A: LMR
    else if (q[0] < q[1] && q[1] < q[2]) {
        return true;
    }
    // Case 3B: RML
    else if (q[2] < q[1] && q[1] < q[0]) {
        return true;
    }
    else {
        return false;
    }
}


// Pass in the triangle to split (original), and then assign the passed-in L and R
// pointers to the data that makes them represent original split in half
//
// Also LERP M2's 3 color values and its z value. Assign each L and R triangle
// the right color and Z value appropriately 
void split(Triangle original, Triangle* triL, Triangle* triR) {
    // Find M2
    double m2_y = endpointY(original.X, original.Y);
    double m2_x = original.X[1];
    double m2_colorR;
    double m2_colorG;
    double m2_colorB;
    double m2_z;
    double m2_shading;

    // Lerp Z for M2
    double t = (m2_x - original.X[0]) / (original.X[2] - original.X[0]);
    m2_z = original.Z[0] + (t * (original.Z[2] - original.Z[0]));

    // Lerp colorR, colorG, colorB for M2
    m2_colorR = original.colors[0][0] + (t * (original.colors[2][0] - original.colors[0][0]));
    m2_colorG = original.colors[0][1] + (t * (original.colors[2][1] - original.colors[0][1]));
    m2_colorB = original.colors[0][2] + (t * (original.colors[2][2] - original.colors[0][2]));
    // 1F: same for shading values
    m2_shading = original.shading[0] + (t * (original.shading[2] - original.shading[0]));
    

    // Right side
    if (m2_y > original.Y[1]) {
        triR->X[0] = m2_x;
        triR->Y[0] = m2_y;
        triR->X[1] = original.X[1];
        triR->Y[1] = original.Y[1];
        triR->X[2] = original.X[2];
        triR->Y[2] = original.Y[2];

        //Z stuff
        triR->Z[0] = m2_z;
        triR->Z[1] = original.Z[1];
        triR->Z[2] = original.Z[2];

        //R
        triR->colors[0][0] = m2_colorR;
        triR->colors[0][1] = m2_colorG;
        triR->colors[0][2] = m2_colorB;
        //G
        triR->colors[1][0] = original.colors[1][0];
        triR->colors[1][1] = original.colors[1][1];
        triR->colors[1][2] = original.colors[1][2];
        //B
        triR->colors[2][0] = original.colors[2][0];
        triR->colors[2][1] = original.colors[2][1];
        triR->colors[2][2] = original.colors[2][2];

        // Shading
        triR->shading[0] = m2_shading;
        triR->shading[1] = original.shading[1];
        triR->shading[2] = original.shading[2];
        
    }
    else {
        triR->X[0] = original.X[1];
        triR->Y[0] = original.Y[1];
        triR->X[1] = m2_x;
        triR->Y[1] = m2_y;
        triR->X[2] = original.X[2];
        triR->Y[2] = original.Y[2];

        //Z stuff
        triR->Z[0] = original.Z[1];
        triR->Z[1] = m2_z;
        triR->Z[2] = original.Z[2];

        //R
        triR->colors[0][0] = original.colors[1][0];
        triR->colors[0][1] = original.colors[1][1];
        triR->colors[0][2] = original.colors[1][2];
        //G
        triR->colors[1][0] = m2_colorR;
        triR->colors[1][1] = m2_colorG;
        triR->colors[1][2] = m2_colorB;
        //B
        triR->colors[2][0] = original.colors[2][0];
        triR->colors[2][1] = original.colors[2][1];
        triR->colors[2][2] = original.colors[2][2];

        // Shading
        triR->shading[0] = original.shading[1];
        triR->shading[1] = m2_shading;
        triR->shading[2] = original.shading[2];

    }
    

    // Left side
    if (m2_y > original.Y[1]) {
        triL->X[0] = m2_x;
        triL->Y[0] = m2_y;
        triL->X[1] = original.X[1];
        triL->Y[1] = original.Y[1];
        triL->X[2] = original.X[0];
        triL->Y[2] = original.Y[0];

        //Z stuff
        triL->Z[0] = m2_z;
        triL->Z[1] = original.Z[1];
        triL->Z[2] = original.Z[0];

        //R
        triL->colors[0][0] = m2_colorR;
        triL->colors[0][1] = m2_colorG;
        triL->colors[0][2] = m2_colorB;
        //G
        triL->colors[1][0] = original.colors[1][0];
        triL->colors[1][1] = original.colors[1][1];
        triL->colors[1][2] = original.colors[1][2];
        //B
        triL->colors[2][0] = original.colors[0][0];
        triL->colors[2][1] = original.colors[0][1];
        triL->colors[2][2] = original.colors[0][2];

        // Shading
        triL->shading[0] = m2_shading;
        triL->shading[1] = original.shading[1];
        triL->shading[2] = original.shading[0];

    }
    else {
        triL->X[0] = original.X[1];
        triL->Y[0] = original.Y[1];
        triL->X[1] = m2_x;
        triL->Y[1] = m2_y;
        triL->X[2] = original.X[0];
        triL->Y[2] = original.Y[0];

        //Z stuff
        triL->Z[0] = original.Z[1];
        triL->Z[1] = m2_z;
        triL->Z[2] = original.Z[0];

        //R
        triL->colors[0][0] = original.colors[1][0];
        triL->colors[0][1] = original.colors[1][1];
        triL->colors[0][2] = original.colors[1][2];
        //G
        triL->colors[1][0] = m2_colorR;
        triL->colors[1][1] = m2_colorG;
        triL->colors[1][2] = m2_colorB;
        //B
        triL->colors[2][0] = original.colors[0][0];
        triL->colors[2][1] = original.colors[0][1];
        triL->colors[2][2] = original.colors[0][2];

        // Shading
        triL->shading[0] = original.shading[1];
        triL->shading[1] = m2_shading;
        triL->shading[2] = original.shading[0];

    }
    
}



/*
// Check to see if two points make a vertical line (have the same X value)
// Then, see if the last point's X value is less than this line, making 
// a going left triangle
//
// If it is, swap the coordinate storage such that each vertex is stored as:
// Upper right = X0, Y0        Lower Right = X1, Y1,       Left = X2, Y2
//
*/
bool isGoingL_setVert(Triangle& tri) {
    if (tri.X[0] == tri.X[1]) {
        if (tri.X[2] <= tri.X[0] && tri.X[2] <= tri.X[1]) {
            if (tri.Y[0] > tri.Y[1]) {
                return true;
            }
            else if (tri.Y[1] > tri.Y[0]) {
                std::swap(tri.X[0], tri.X[1]);
                std::swap(tri.Y[0], tri.Y[1]);
                std::swap(tri.Z[0], tri.Z[1]);
                std::swap(tri.shading[0], tri.shading[1]);
                for (int i = 0; i < 3; i++) {
                    std::swap(tri.colors[0][i], tri.colors[1][i]);
                }
                return true;
            }
            else {
                return false;
            }
        }
        else {
            return false;
        }
    }
    else if (tri.X[1] == tri.X[2]) {
        if (tri.X[0] <= tri.X[1] && tri.X[0] <= tri.X[2]) {
            if (tri.Y[1] > tri.Y[2]) {

                std::swap(tri.X[1], tri.X[2]);
                std::swap(tri.Y[1], tri.Y[2]);
                std::swap(tri.Z[1], tri.Z[2]);
                std::swap(tri.shading[1], tri.shading[2]);
                for (int i = 0; i < 3; i++) {
                    std::swap(tri.colors[1][i], tri.colors[2][i]);
                }

                std::swap(tri.X[0], tri.X[2]);
                std::swap(tri.Y[0], tri.Y[2]);
                std::swap(tri.Z[0], tri.Z[2]);
                std::swap(tri.shading[0], tri.shading[2]);
                for (int i = 0; i < 3; i++) {
                    std::swap(tri.colors[0][i], tri.colors[2][i]);
                }

                return true;
            }
            else if (tri.Y[2] > tri.Y[1]) {
                std::swap(tri.X[0], tri.X[2]);
                std::swap(tri.Y[0], tri.Y[2]);
                std::swap(tri.Z[0], tri.Z[2]);
                std::swap(tri.shading[0], tri.shading[2]);
                for (int i = 0; i < 3; i++) {
                    std::swap(tri.colors[0][i], tri.colors[2][i]);
                }
                return true;
            }
            else {
                return false;
            }
            
        }
        else {
            return false;
        }
    }
    else if (tri.X[2] == tri.X[0]) {
        if (tri.X[1] <= tri.X[2] && tri.X[1] <= tri.X[0]) {
            if (tri.Y[2] > tri.Y[0]) {
                std::swap(tri.X[0], tri.X[2]);
                std::swap(tri.Y[0], tri.Y[2]);
                std::swap(tri.Z[0], tri.Z[2]);
                std::swap(tri.shading[0], tri.shading[2]);
                for (int i = 0; i < 3; i++) {
                    std::swap(tri.colors[0][i], tri.colors[2][i]);
                }

                std::swap(tri.X[1], tri.X[2]);
                std::swap(tri.Y[1], tri.Y[2]);
                std::swap(tri.Z[1], tri.Z[2]);
                std::swap(tri.shading[1], tri.shading[2]);
                for (int i = 0; i < 3; i++) {
                    std::swap(tri.colors[1][i], tri.colors[2][i]);
                }
                return true;
            }
            else if (tri.Y[0] > tri.Y[2]) {
                std::swap(tri.X[1], tri.X[2]);
                std::swap(tri.Y[1], tri.Y[2]);
                std::swap(tri.Z[1], tri.Z[2]);
                std::swap(tri.shading[1], tri.shading[2]);
                for (int i = 0; i < 3; i++) {
                    std::swap(tri.colors[1][i], tri.colors[2][i]);
                }
                return true;
            }
            else {
                return false;
            }
            
        }
        else {
            return false;
        }
    }
    else {
        return false;
    }
    
}

/*
// Check to see if two points make a vertical line (have the same X value)
// Then, see if the last point's X value is greater than this line, making
// a going right triangle
// 
// If it is, swap the coordinates such that they are sored in the following format:
// Upper left = X0, Y0    Lower left = X1, Y1,    Right = X2, Y2
//
*/

bool isGoingR_setVert(Triangle& tri) {
    if (tri.X[0] == tri.X[1]) {
        if (tri.X[2] >= tri.X[0] && tri.X[2] >= tri.X[0]) {
            if (tri.Y[0] > tri.Y[1]) {

                return true;
            }
            else if (tri.Y[1] > tri.Y[0]) {
                std::swap(tri.X[0], tri.X[1]);
                std::swap(tri.Y[0], tri.Y[1]);
                std::swap(tri.Z[0], tri.Z[1]);
                std::swap(tri.shading[0], tri.shading[2]);
                for (int i = 0; i < 3; i++) {
                    std::swap(tri.colors[0][i], tri.colors[1][i]);
                }
                return true;
            }
            else {
                return false;
            }
        }
        else {
            return false;
        }
    }
    else if (tri.X[1] == tri.X[2]) {
        if (tri.X[0] >= tri.X[1] && tri.X[0] >= tri.X[2]) {
            if (tri.Y[1] > tri.Y[2]) {
                std::swap(tri.X[1], tri.X[2]);
                std::swap(tri.Y[1], tri.Y[2]);
                std::swap(tri.Z[1], tri.Z[2]);
                std::swap(tri.shading[1], tri.shading[2]);
                for (int i = 0; i < 3; i++) {
                    std::swap(tri.colors[1][i], tri.colors[2][i]);
                }

                std::swap(tri.X[0], tri.X[2]);
                std::swap(tri.Y[0], tri.Y[2]);
                std::swap(tri.Z[0], tri.Z[2]);
                std::swap(tri.shading[0], tri.shading[2]);
                for (int i = 0; i < 3; i++) {
                    std::swap(tri.colors[0][i], tri.colors[2][i]);
                }
                return true;
            }
            else if (tri.Y[2] > tri.Y[1]) {
                std::swap(tri.X[0], tri.X[2]);
                std::swap(tri.Y[0], tri.Y[2]);
                std::swap(tri.Z[0], tri.Z[2]);
                std::swap(tri.shading[0], tri.shading[2]);
                for (int i = 0; i < 3; i++) {
                    std::swap(tri.colors[0][i], tri.colors[2][i]);
                }

                return true;
            }
            else {
                return false;
            }

        }
        else {
            return false;
        }
    }
    else if (tri.X[2] == tri.X[0]) {
        if (tri.X[1] >= tri.X[2] && tri.X[1] >= tri.X[0]) {
            if (tri.Y[2] > tri.Y[0]) {
                std::swap(tri.X[0], tri.X[2]);
                std::swap(tri.Y[0], tri.Y[2]);
                std::swap(tri.Z[0], tri.Z[2]);
                std::swap(tri.shading[0], tri.shading[2]);
                for (int i = 0; i < 3; i++) {
                    std::swap(tri.colors[0][i], tri.colors[2][i]);
                }
                std::swap(tri.X[1], tri.X[2]);
                std::swap(tri.Y[1], tri.Y[2]);
                std::swap(tri.Z[1], tri.Z[2]);
                std::swap(tri.shading[1], tri.shading[2]);
                for (int i = 0; i < 3; i++) {
                    std::swap(tri.colors[1][i], tri.colors[2][i]);
                }
                return true;
            }
            else if (tri.Y[0] > tri.Y[2]) {
                std::swap(tri.X[1], tri.X[2]);
                std::swap(tri.Y[1], tri.Y[2]);
                std::swap(tri.Z[1], tri.Z[2]);
                std::swap(tri.shading[1], tri.shading[2]);
                for (int i = 0; i < 3; i++) {
                    std::swap(tri.colors[1][i], tri.colors[2][i]);
                }
                return true;
            }
            else {
                return false;
            }

        }
        else {
            return false;
        }
    }
    else {
        return false;
    }
}

void transformTrianglesToDeviceSpace(Triangle &tri, Matrix m) {
    // The triangle member arrays for coordinates are by coordinate, then point instead
    // of point, then coordinate. Create a 2D array that represents as ret[point number][x, y, z, or w coordinate]
    double** ret = new double*[3];
    for (int i = 0; i < 3; i++) {
        ret[i] = new double[4];
    }
    // ret[0] is the first point in (x, y, z, w) format
    ret[0][0] = tri.X[0];
    ret[0][1] = tri.Y[0];
    ret[0][2] = tri.Z[0];
    ret[0][3] = 1;
    // ret[1] is the second point of the triangle in the same format
    ret[1][0] = tri.X[1];
    ret[1][1] = tri.Y[1];
    ret[1][2] = tri.Z[1];
    ret[1][3] = 1;
    // ret[2] is the third point of the triangle in the same format
    ret[2][0] = tri.X[2];
    ret[2][1] = tri.Y[2];
    ret[2][2] = tri.Z[2];
    ret[2][3] = 1;

    // Multiply each point in the triangle by the m matrix to get its place in
    // device space
    double *p1prime = new double[4];
    double *p2prime = new double[4];
    double *p3prime = new double[4];
    m.TransformPoint(ret[0], p1prime);
    m.TransformPoint(ret[1], p2prime);
    m.TransformPoint(ret[2], p3prime);

    // Now, simultaneously,
    // set the triangle X Y and Z arrays to their new device space values
    // while also dividing by m, the last coordinate
    tri.X[0] = p1prime[0] / p1prime[3];
    tri.Y[0] = p1prime[1] / p1prime[3];
    tri.Z[0] = p1prime[2] / p1prime[3];

    tri.X[1] = p2prime[0] / p2prime[3];
    tri.Y[1] = p2prime[1] / p2prime[3];
    tri.Z[1] = p2prime[2] / p2prime[3];

    tri.X[2] = p3prime[0] / p3prime[3];
    tri.Y[2] = p3prime[1] / p3prime[3];
    tri.Z[2] = p3prime[2] / p3prime[3];
}


double calculatePhongShading(LightingParameters lightP, double diffuse, double specular) {
    double ret = lightP.Ka + (lightP.Kd * diffuse) + (lightP.Ks * specular);
    
    return ret;
}

// Return a pointer to an array that holds the Specular value for phong shading, in the 0, 1, 2 vertex order
double* calculateSpecularShading(LightingParameters lightP, Camera cam, Triangle tri) {
    // Let's do the "right branch" of my drawing for this step first
    double N0[3];
    double N1[3];
    double N2[3];
    double L[3];
    double viewDir0[3];
    double viewDir1[3];
    double viewDir2[3];
    double* ret = new double[3];
    double alpha = lightP.alpha;

    // Make copies of normals and lightDir as N0, N1, N2 and L, matching the naming conventions
    N0[0] = tri.normals[0][0];
    N0[1] = tri.normals[0][1];
    N0[2] = tri.normals[0][2];

    N1[0] = tri.normals[1][0];
    N1[1] = tri.normals[1][1];
    N1[2] = tri.normals[1][2];

    N2[0] = tri.normals[2][0];
    N2[1] = tri.normals[2][1];
    N2[2] = tri.normals[2][2];


    L[0] = lightP.lightDir[0];
    L[1] = lightP.lightDir[1];
    L[2] = lightP.lightDir[2];

    // Normalize the N's
    double N1norm[3];
    double N2norm[3];
    double N0norm[3];

    normalizeVector(N0, N0norm);
    normalizeVector(N1, N1norm);
    normalizeVector(N2, N2norm);

    // Dot produc the N's and L
    double LNdot0 = dotProduct(L, N0norm);
    double LNdot1 = dotProduct(L, N1norm);
    double LNdot2 = dotProduct(L, N2norm);

    // 2 * (dotProduct(L, N's)
    LNdot0 *= 2;
    LNdot1 *= 2;
    LNdot2 *= 2;

    // 2 * (dotProduct(L, N's) * N
    scalarProduct(LNdot0, N0norm);
    scalarProduct(LNdot1, N1norm);
    scalarProduct(LNdot2, N2norm);

    // 2 * (dotProduct(L, N's) * N - L
    vectorSubtract(N0norm, L);
    vectorSubtract(N1norm, L);
    vectorSubtract(N2norm, L);

    // Now, the N's are the R's by my convention in my notes. Let's do the left branch now.

    // Calculate the viewDir's
    viewDir0[0] = cam.position[0] - tri.X[0];
    viewDir0[1] = cam.position[1] - tri.Y[0];
    viewDir0[2] = cam.position[2] - tri.Z[0];

    viewDir1[0] = cam.position[0] - tri.X[1];
    viewDir1[1] = cam.position[1] - tri.Y[1];
    viewDir1[2] = cam.position[2] - tri.Z[1];

    viewDir2[0] = cam.position[0] - tri.X[2];
    viewDir2[1] = cam.position[1] - tri.Y[2];
    viewDir2[2] = cam.position[2] - tri.Z[2];

    // Normalize the viewDir's
    double viewDir0norm[3];
    double viewDir1norm[3];
    double viewDir2norm[3];
    normalizeVector(viewDir0, viewDir0norm);
    normalizeVector(viewDir1, viewDir1norm);
    normalizeVector(viewDir2, viewDir2norm);

    // Now the two branches from notes meet
    double shading0 = dotProduct(viewDir0norm, N0norm);
    double shading1 = dotProduct(viewDir1norm, N1norm);
    double shading2 = dotProduct(viewDir2norm, N2norm);

    maxWithZero(shading0);
    maxWithZero(shading1);
    maxWithZero(shading2);

    ret[0] = pow(shading0, alpha);
    ret[1] = pow(shading1, alpha);
    ret[2] = pow(shading2, alpha);

    return ret;
}

// Returns a pointer to an array that holds the Diffuse values for each vertex in order of 0, 1, 2
// These diffuse values are what we use in calculatePhongShading() to calculate the overall shading value for each vertex
double* calculateDiffuseShading(LightingParameters lightP, Triangle tri) {
    double N0[3];
    double N1[3];
    double N2[3];
    double L[3];
    double* ret = new double[3];

    // Make copies of normals and lightDir as N0, N1, N2 and L, matching the naming conventions
    N0[0] = tri.normals[0][0];
    N0[1] = tri.normals[0][1];
    N0[2] = tri.normals[0][2];

    N1[0] = tri.normals[1][0];
    N1[1] = tri.normals[1][1];
    N1[2] = tri.normals[1][2];

    N2[0] = tri.normals[2][0];
    N2[1] = tri.normals[2][1];
    N2[2] = tri.normals[2][2];
    

    L[0] = lightP.lightDir[0];
    L[1] = lightP.lightDir[1];
    L[2] = lightP.lightDir[2];

    // Normalize the N's
    double N1norm[3];
    double N2norm[3];
    double N0norm[3];

    normalizeVector(N0, N0norm);
    normalizeVector(N1, N1norm);
    normalizeVector(N2, N2norm);

    // Dot Product N's with L
    double shade0 = dotProduct(L, N0norm);
    double shade1 = dotProduct(L, N1norm);
    double shade2 = dotProduct(L, N2norm);

    // Set to zero if zero is larger than shade values
    maxWithZero(shade0);
    maxWithZero(shade1);
    maxWithZero(shade2);

    ret[0] = shade0;
    ret[1] = shade1;
    ret[2] = shade2;

    return ret;
}



int main()
{
    // Image is 1000x1000
    vtkImageData* image = NewImage(1000, 1000);
    unsigned char* buffer =
        (unsigned char*)image->GetScalarPointer(0, 0, 0);
    int npixels = 1000 * 1000;
    for (int i = 0; i < npixels * 3; i++)
        buffer[i] = 0;

    std::vector<Triangle> triangles = GetTriangles();

    Screen screen;
    screen.buffer = buffer;
    screen.width = 1000;
    screen.height = 1000;

    // Set the depthBuffer to be a 2D array of size [screenwidth][screenheight]
    // And set every element in it to -1
    double **depthBuffer = new double*[screen.width+1];
    for (int i = 0; i < screen.width; i++) {
        depthBuffer[i] = new double[screen.height+1];
        for (int j = 0; j < screen.height; j++) {
            depthBuffer[i][j] = -1;
        }
    }

    int numTriangles = triangles.size();

   
    //std::string s;
    //const char* ccc;


    for (int k = 0; k < 1000; k++) {



        // Rest depth buffer and buffer for a fresh new image
        for (int i = 0; i < screen.width; i++) {
            depthBuffer[i] = new double[screen.height + 1];
            for (int j = 0; j < screen.height; j++) {
                depthBuffer[i][j] = -1;
            }
        }

        for (int i = 0; i < npixels * 3; i++)
            buffer[i] = 0;


        int frame = k;
        Camera c = GetCamera(frame, 1000);

        //The four transformation steps, one Camera Frame setup and 3 matrix creations
        c.CameraFrameCalc();
        Matrix camTrans = c.CameraTransform();
        Matrix viewTrans = c.ViewTransform();
        Matrix devTrans = c.DeviceTransform();


        // Multiply all three matrices together to create m
        Matrix temp = temp.ComposeMatrices(camTrans, viewTrans);
        Matrix m = m.ComposeMatrices(temp, devTrans);

        // 1F
        LightingParameters lightP = GetLighting(c);


        for (int i = 0; i < numTriangles; i++) {

            Triangle tri = triangles[i];


            // add CalculatePhongShading() here?
            double* diffuse = calculateDiffuseShading(lightP, tri);
            double* specular = calculateSpecularShading(lightP, c, tri);


            // Calculate the final value of Phong shading and set it to the shading array in the triangle object
            tri.shading[0] = calculatePhongShading(lightP, diffuse[0], specular[0]);
            tri.shading[1] = calculatePhongShading(lightP, diffuse[1], specular[1]);
            tri.shading[2] = calculatePhongShading(lightP, diffuse[2], specular[2]);


            // The culmination of all the transforming work, apply it to the triangle
            transformTrianglesToDeviceSpace(tri, m);



            // Is the triangle arbitrary?
            if (isArb(tri.X, tri.Y)) {

                tri.setVertArb(tri.X, tri.Y, tri.Z, tri.colors, tri.shading); // in the 0 1 2 or L M R format. not split yet

                Triangle* triL = new Triangle();
                Triangle* triR = new Triangle();
                split(tri, triL, triR);

                rasterizeGoingL(*triL, buffer, depthBuffer);
                rasterizeGoingR(*triR, buffer, depthBuffer);
            }
            // Ok, it's not, it's already going R or L
            else {
                // Is it going L? Checking if it is will also set the vertices to the right format (check function comments for more)
                if (isGoingL_setVert(tri)) {

                    rasterizeGoingL(tri, buffer, depthBuffer);
                }
                // Is it going R? Checking if it is will also set the vertices to the right format (check function comments for more)
                else if (isGoingR_setVert(tri)) {

                    rasterizeGoingR(tri, buffer, depthBuffer);
                }
                else {
                }
            }

        }

        // Terrible code for getting multiple file name outputs
        /*
        std::string s = "";
        int j = k;
        std::string s1 = std::to_string(j);

        if (k == 0) {
            s1 = "0";
        }
        s += s1;
        */
        //s = std::to_string(k);
        //ccc = s.c_str();
        char str[128];
        sprintf(str, "%03d", k);
        WriteImage(image, str);
    }

}
