#include <math.h>
#include <iostream>
#include <algorithm>
#include <vtkDataSet.h>
#include <vtkImageData.h>
#include <vtkPNGWriter.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPoints.h>
#include <vtkUnsignedCharArray.h>
#include <vtkFloatArray.h>
#include <vtkDoubleArray.h>
#include <vtkCellArray.h>
#include <vtkDataSetWriter.h>

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

using std::cerr;
using std::endl;
using std::min;
using std::max;

class Triangle {
  public:
      double          X[3];
      double          Y[3];
      double          Z[3];
      double  colors[3][3];
      double normals[3][3];
};

struct LightingParameters {
    LightingParameters(void) {
         lightDir[0] = -0.6;
         lightDir[1] = 0;
         lightDir[2] = -0.8;
         Ka = 0.3;
         Kd = 0.7;
         Ks = 5.3;
         alpha = 7.5;
    };
  
    double lightDir[3]; // The direction of the light source
    double Ka;           // The coefficient for ambient lighting.
    double Kd;           // The coefficient for diffuse lighting.
    double Ks;           // The coefficient for specular lighting.
    double alpha;        // The exponent term for specular lighting.
};

LightingParameters lp;

std::vector<Triangle> GetTriangles(void) {
    vtkPolyDataReader *rdr = vtkPolyDataReader::New();
    rdr->SetFileName("proj1e_geometry.vtk");
    cerr << "Reading" << endl;
    rdr->Update();
    cerr << "Done reading" << endl;
    if (rdr->GetOutput()->GetNumberOfCells() == 0) {
        cerr << "Unable to open file!!" << endl;
        exit(EXIT_FAILURE);
    }
    vtkPolyData *pd = rdr->GetOutput();

    int numTris = pd->GetNumberOfCells();
    vtkPoints *pts = pd->GetPoints();
    vtkCellArray *cells = pd->GetPolys();
    vtkDoubleArray *var = (vtkDoubleArray *) pd->GetPointData()->GetArray("hardyglobal");
    double *color_ptr = var->GetPointer(0);
    vtkFloatArray *n = (vtkFloatArray *) pd->GetPointData()->GetNormals();
    float *normals = n->GetPointer(0);
    std::vector<Triangle> tris(numTris);
    vtkIdType npts;
    vtkIdType *ptIds;
    int idx;
    for (idx = 0, cells->InitTraversal() ; cells->GetNextCell(npts, ptIds) ; idx++)
    {
        if (npts != 3)
        {
            cerr << "Non-triangles!! ???" << endl;
            exit(EXIT_FAILURE);
        }
        double *pt = NULL;
        pt = pts->GetPoint(ptIds[0]);
        tris[idx].X[0] = (pt[0]+10)*50.0;
        tris[idx].Y[0] = (pt[1]+10)*50.0;
        tris[idx].Z[0] = (pt[2]-10)*0.05;
        tris[idx].normals[0][0] = normals[3*ptIds[0]+0];
        tris[idx].normals[0][1] = normals[3*ptIds[0]+1];
        tris[idx].normals[0][2] = normals[3*ptIds[0]+2];
        pt = pts->GetPoint(ptIds[1]);
        tris[idx].X[1] = (pt[0]+10)*50.0;
        tris[idx].Y[1] = (pt[1]+10)*50.0;
        tris[idx].Z[1] = (pt[2]-10)*0.05;
        tris[idx].Y[1] = (pt[1]+10)*50.0;
        tris[idx].Z[1] = (pt[2]-10)*0.05;
        tris[idx].normals[1][0] = normals[3*ptIds[1]+0];
        tris[idx].normals[1][1] = normals[3*ptIds[1]+1];
        tris[idx].normals[1][2] = normals[3*ptIds[1]+2];
        pt = pts->GetPoint(ptIds[2]);
        tris[idx].X[2] = (pt[0]+10)*50.0;
        tris[idx].Y[2] = (pt[1]+10)*50.0;
        tris[idx].Z[2] = (pt[2]-10)*0.05;
        tris[idx].normals[2][0] = normals[3*ptIds[2]+0];
        tris[idx].normals[2][1] = normals[3*ptIds[2]+1];
        tris[idx].normals[2][2] = normals[3*ptIds[2]+2];

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
        for (int j = 0 ; j < 3 ; ++j)
        {   
            float val = color_ptr[ptIds[j]];
            int r;
            for (r = 0 ; r < 7 ; r++)
            {
                if (mins[r] <= val && val < maxs[r])
                    break;
            }
            if (r == 7)
            {
                cerr << "Could not interpolate color for " << val << endl;
                exit(EXIT_FAILURE);
            }
            double proportion = (val-mins[r]) / (maxs[r]-mins[r]);
            tris[idx].colors[j][0] = (RGB[r][0]+proportion*(RGB[r+1][0]-RGB[r][0]))/255.0;
            tris[idx].colors[j][1] = (RGB[r][1]+proportion*(RGB[r+1][1]-RGB[r][1]))/255.0;
            tris[idx].colors[j][2] = (RGB[r][2]+proportion*(RGB[r+1][2]-RGB[r][2]))/255.0;
        }
    }

    return tris;
}

double dotProduct(double A[], double B[]) {
    double dot = 0.0;
    for (int i = 0; i < 3; ++i) {
        dot += (A[i] * B[i]);
    }
    return dot;
}

void normalize(double A[]) {
    double normal = 0.0;
    for (int i = 0; i < 3; ++i) {
        normal += (A[i] * A[i]);
    }
    normal = sqrt(normal);
    for (int i = 0; i < 3; ++i) {
        A[i] /= normal;
    }
}

double ceil441(double f) {
    return ceil(f-0.00001);
}

double floor441(double f) {
    return floor(f+0.00001);
}

double MIN(double *a) {
    double smallest = a[0];
    for (int i = 0; i < 3; ++i) {
        if (a[i] < smallest) {
            smallest = a[i];
        }
    }
    return smallest;
}

double MAX(double *a) {
    double max = a[0];
    for (int i = 0; i < 3; ++i) {
        if (a[i] > max) {
            max = a[i];
        }
    }
    return max;
}

vtkImageData *
NewImage(int width, int height)
{
    vtkImageData *img = vtkImageData::New();
    img->SetDimensions(width, height, 1);
    img->AllocateScalars(VTK_UNSIGNED_CHAR, 3);

    return img;
}

void
WriteImage(vtkImageData *img, const char *filename)
{
   std::string full_filename = filename;
   full_filename += ".png";
   vtkPNGWriter *writer = vtkPNGWriter::New();
   writer->SetInputData(img);
   writer->SetFileName(full_filename.c_str());
   writer->Write();
   writer->Delete();
}

class Screen {
  public:
      unsigned char   *buffer;
      double          *zbuffer;
      int width, height;
};

double interpolate(double start, double end, double step) {
    return start + (end - start) * (step);
}

double specular(double L[3], double N[3]) {
    double V[3] = { 0, 0, -1.0 };
    double s = 2 * dotProduct(L, N); 
    double R[3] = {(s * N[0]) - L[0], (s * N[1]) - L[1], (s * N[2]) - L[2]};
    double RdotV = dotProduct(R, V);
    double spec = lp.Ks * pow(RdotV, lp.alpha);
    return std::max(0., spec);
}

double calculateShading(double *viewDirection, double *normal) {
    normalize(normal);
    double diffuse = fabs(dotProduct(lp.lightDir, normal));
    double temp = 2 * dotProduct(lp.lightDir, normal);
    double R[3];
    for (int i = 0; i < 3; ++i) {
        R[i] = temp * normal[i] - lp.lightDir[i];
    }
    double specular = fmax(0, pow(dotProduct(R, viewDirection), lp.alpha));
    return lp.Ka + lp.Kd * diffuse + lp.Ks * specular;
}

void draw(double leftEnd, double yMin, Screen s, double colors[3], double norm[3]) {
     int idx = yMin * s.width * 3 + leftEnd * 3;
     double V[3] = { 0, 0, -1.0 };
     colors[0] = std::max(0.0, (calculateShading(V, norm)) * colors[0]);
     colors[1] = std::max(0.0, (calculateShading(V, norm)) * colors[1]);
     colors[2] = std::max(0.0, (calculateShading(V, norm)) * colors[2]);
     s.buffer[idx]   = std::min(255.0, ceil441(colors[0] * 255));
     s.buffer[idx+1] = std::min(255.0, ceil441(colors[1] * 255));
     s.buffer[idx+2] = std::min(255.0, ceil441(colors[2] * 255));
}

void sortVertices(double peak[], double right[], double left[], double pNorm[], double rNorm[], double lNorm[], Triangle t) {
    for (int i = 0; i < 3; ++i) {
        if (t.Y[i] == peak[1]) {
            peak[0] = t.X[i];
            peak[1] = t.Y[i];
            peak[2] = t.Z[i];
            peak[3] = t.colors[i][0];
            peak[4] = t.colors[i][1];
            peak[5] = t.colors[i][2];
            pNorm[0] = t.normals[i][0];
            pNorm[1] = t.normals[i][1];
            pNorm[2] = t.normals[i][2];
        }
    } 
    for (int i = 0; i < 3; ++i) {
        if (t.Y[i] != peak[1]) {
            if (t.X[i] > right[0]) {
                right[0] = t.X[i];
                right[1] = t.Y[i];
                right[2] = t.Z[i];
                right[3] = t.colors[i][0];
                right[4] = t.colors[i][1];
                right[5] = t.colors[i][2];
                rNorm[0] = t.normals[i][0];
                rNorm[1] = t.normals[i][1];
                rNorm[2] = t.normals[i][2];
            }
        }
    }
    for (int i = 0; i < 3; ++i) {
        if (t.Y[i] != peak[1]) {
            if (t.X[i] < right[0]) {
                left[0] = t.X[i];
                left[1] = t.Y[i];
                left[2] = t.Z[i];
                left[3] = t.colors[i][0];
                left[4] = t.colors[i][1];
                left[5] = t.colors[i][2];
                lNorm[0] = t.normals[i][0];
                lNorm[1] = t.normals[i][1];
                lNorm[2] = t.normals[i][2];
            }
        }
    }
} 

void flatBottom(Triangle t, Screen s) {
    double top[6]   = {MAX(t.X), MAX(t.Y), 0, 0, 0, 0};
    double right[6] = {MIN(t.X), MIN(t.Y), 0, 0, 0, 0};
    double left[6]  = {MIN(t.X), MIN(t.Y), 0, 0, 0, 0};
    double tNorm[3];
    double rNorm[3];
    double lNorm[3];
    sortVertices(top, right, left, tNorm, rNorm, lNorm, t);
    normalize(tNorm); 
    normalize(rNorm); 
    normalize(lNorm); 
    double leftEnd;
    double rightEnd;
    double yMin = ceil441(left[1]);
    double yMax = floor441(top[1]);

    if (top[1] >= s.height) { yMax = s.height-1; }
    if (yMin < 0) { yMin = 0; }
    for (; yMin <= yMax && yMin < s.height; ++yMin) {
        double step = (yMin - left[1]) / (top[1] - left[1]);
        double leftZ = interpolate(left[2], top[2], step);
        double rightZ = interpolate(right[2], top[2], step);
        double leftBound = interpolate(left[0], top[0], step);
        double rightBound = interpolate(right[0], top[0], step);
        double leftNorm[3] = { interpolate(lNorm[0], tNorm[0], step),
                               interpolate(lNorm[1], tNorm[1], step),
                               interpolate(lNorm[2], tNorm[2], step)};
        double rightNorm[3] = { interpolate(rNorm[0], tNorm[0], step),
                                interpolate(rNorm[1], tNorm[1], step),
                                interpolate(rNorm[2], tNorm[2], step)};
        normalize(leftNorm);
        normalize(rightNorm);
        leftEnd = ceil441(leftBound);
        rightEnd = floor441(rightBound);
        if (rightEnd >= s.width) { rightEnd = s.width - 1; }
        if (leftEnd < 0) { leftEnd = 0; }
        double leftColor[3]  =  { interpolate(left[3], top[3], step),
                                  interpolate(left[4], top[4], step),
                                  interpolate(left[5], top[5], step)};
        double rightColor[3] = { interpolate(right[3], top[3], step),
                                 interpolate(right[4], top[4], step),
                                 interpolate(right[5], top[5], step)};
        for (; leftEnd <= rightEnd && leftEnd < s.width; ++leftEnd) {
            step = (leftEnd - leftBound) / (rightBound - leftBound);
            double zValue = interpolate(leftZ, rightZ, step);
            double pNorm[3] = { interpolate(leftNorm[0], rightNorm[0], step),
                                interpolate(leftNorm[1], rightNorm[1], step),
                                interpolate(leftNorm[2], rightNorm[2], step)};
            normalize(pNorm);
            int idx = (int)yMin * s.width + (int)leftEnd;
            if (s.zbuffer[idx] > zValue) {
                continue;
            }
            double pColor[3] = { interpolate(leftColor[0], rightColor[0], step),
                                 interpolate(leftColor[1], rightColor[1], step),
                                 interpolate(leftColor[2], rightColor[2], step)};
            s.zbuffer[idx] = zValue;
            draw(leftEnd, yMin, s, pColor, pNorm);
        }
    }
}

void flatTop(Triangle t, Screen s) {
    double bot[6]   = {MIN(t.X), MIN(t.Y), 0, 0, 0, 0};
    double right[6] = {MIN(t.X), MIN(t.Y), 0, 0, 0, 0};
    double left[6]  = {MIN(t.X), MIN(t.Y), 0, 0, 0, 0};
    double bNorm[3];
    double rNorm[3];
    double lNorm[3];
    sortVertices(bot, right, left, bNorm, rNorm, lNorm, t);
    normalize(bNorm); 
    normalize(rNorm); 
    normalize(lNorm); 
    double leftEnd;
    double rightEnd;
    double yMax = floor441(left[1]); 
    double yMin = ceil441(bot[1]);

    if (yMax >= s.height) { yMax = s.height - 1; }
    if (yMin < 0) { yMin = 0; }
    for (; yMin <= yMax && yMin < s.height ; ++yMin) {
        double step = (yMin - bot[1]) / (left[1] - bot[1]);
        double leftZ = interpolate(bot[2], left[2], step);
        double rightZ = interpolate(bot[2], right[2], step);
        double leftBound = interpolate(bot[0], left[0], step);
        double rightBound = interpolate(bot[0], right[0], step);
        double leftNorm[3] = { interpolate(bNorm[0], lNorm[0], step),
                               interpolate(bNorm[1], lNorm[1], step),
                               interpolate(bNorm[2], lNorm[2], step)};
        double rightNorm[3] = { interpolate(bNorm[0], rNorm[0], step),
                                interpolate(bNorm[1], rNorm[1], step),
                                interpolate(bNorm[2], rNorm[2], step)};
        normalize(leftNorm);
        normalize(rightNorm);
        leftEnd = ceil441(leftBound);
        rightEnd = floor441(rightBound);
        if (rightEnd >= s.width) { rightEnd = s.width - 1; }
        if (leftEnd < 0) { leftEnd = 0; }
        double leftColor[3]  =  { interpolate(bot[3], left[3], step),
                                  interpolate(bot[4], left[4], step),
                                  interpolate(bot[5], left[5], step)};
        double rightColor[3] = { interpolate(bot[3], right[3], step),
                                 interpolate(bot[4], right[4], step),
                                 interpolate(bot[5], right[5], step)};
        for (; leftEnd <= rightEnd && leftEnd < s.width; ++leftEnd) {
            step = (leftEnd - leftBound) / (rightBound - leftBound);
            double zValue = interpolate(leftZ, rightZ, step);
            double pNorm[3] = { interpolate(leftNorm[0], rightNorm[0], step),
                                interpolate(leftNorm[1], rightNorm[1], step),
                                interpolate(leftNorm[2], rightNorm[2], step)};
            normalize(pNorm);
            int idx = (int)yMin * s.width + (int)leftEnd;
            if (s.zbuffer[idx] > zValue) {
                continue;
            }
            double pColor[3] = { interpolate(leftColor[0], rightColor[0], step),
                                 interpolate(leftColor[1], rightColor[1], step),
                                 interpolate(leftColor[2], rightColor[2], step)};
            s.zbuffer[idx] = zValue;
            draw(leftEnd, yMin, s, pColor, pNorm);
        }
    }
}

void split(Triangle t, Screen s) {
    double top[6] = {MIN(t.X), MIN(t.Y), 0, 0, 0, 0};
    double split[6] = {MIN(t.X), MIN(t.Y), 0, 0, 0, 0};
    double nVertex[6] = {MIN(t.X), MIN(t.Y), 0, 0, 0, 0};
    double bot[6] = {MAX(t.X), MAX(t.Y), 0, 0, 0, 0};
    double tNorm[3];
    double sNorm[3];
    double nNorm[3];
    double bNorm[3];

    for (int i = 0; i < 3; ++i) {
        if (t.Y[i] > top[1]) {
            top[0] = t.X[i];
            top[1] = t.Y[i];
            top[2] = t.Z[i];
            top[3] = t.colors[i][0];
            top[4] = t.colors[i][1];
            top[5] = t.colors[i][2];
            tNorm[0] = t.normals[i][0];
            tNorm[1] = t.normals[i][1];
            tNorm[2] = t.normals[i][2];
        }
        if (t.Y[i] < bot[1]) {
            bot[0] = t.X[i];
            bot[1] = t.Y[i];
            bot[2] = t.Z[i];
            bot[3] = t.colors[i][0];
            bot[4] = t.colors[i][1];
            bot[5] = t.colors[i][2];
            bNorm[0] = t.normals[i][0];
            bNorm[1] = t.normals[i][1];
            bNorm[2] = t.normals[i][2];
        }
    }
    for (int i = 0; i < 3; ++i) {
        if (t.Y[i] < top[1] && t.Y[i] > bot[1]) {
            split[0] = t.X[i];
            split[1] = t.Y[i];
            split[2] = t.Z[i];
            split[3] = t.colors[i][0];
            split[4] = t.colors[i][1];
            split[5] = t.colors[i][2];
            sNorm[0] = t.normals[i][0];
            sNorm[1] = t.normals[i][1];
            sNorm[2] = t.normals[i][2];
        }
    }
    Triangle splitTriangle;
    double step = (split[1] - bot[1]) / (top[1] - bot[1]);
    nVertex[0] = interpolate(bot[0], top[0], step);
    nVertex[1] = split[1];
    nVertex[2] = interpolate(bot[2], top[2], step);
    nVertex[3] = interpolate(bot[3], top[3], step);
    nVertex[4] = interpolate(bot[4], top[4], step);
    nVertex[5] = interpolate(bot[5], top[5], step);
    nNorm[0] = interpolate(bNorm[0], tNorm[0], step);
    nNorm[1] = interpolate(bNorm[1], tNorm[1], step);
    nNorm[2] = interpolate(bNorm[2], tNorm[2], step);
    splitTriangle.X[0] = top[0];
    splitTriangle.X[1] = split[0];
    splitTriangle.X[2] = nVertex[0]; 
    splitTriangle.Y[0] = top[1];
    splitTriangle.Y[1] = split[1];
    splitTriangle.Y[2] = nVertex[1];
    splitTriangle.Z[0] = top[2];
    splitTriangle.Z[1] = split[2];
    splitTriangle.Z[2] = nVertex[2]; 
    splitTriangle.colors[0][0] = top[3];
    splitTriangle.colors[0][1] = top[4];
    splitTriangle.colors[0][2] = top[5];
    splitTriangle.colors[1][0] = split[3];
    splitTriangle.colors[1][1] = split[4];
    splitTriangle.colors[1][2] = split[5];
    splitTriangle.colors[2][0] = nVertex[3];
    splitTriangle.colors[2][1] = nVertex[4];
    splitTriangle.colors[2][2] = nVertex[5];
    splitTriangle.normals[0][0] = tNorm[0];
    splitTriangle.normals[0][1] = tNorm[1];
    splitTriangle.normals[0][2] = tNorm[2];
    splitTriangle.normals[1][0] = sNorm[0];
    splitTriangle.normals[1][1] = sNorm[1];
    splitTriangle.normals[1][2] = sNorm[2];
    splitTriangle.normals[2][0] = nNorm[0];
    splitTriangle.normals[2][1] = nNorm[1];
    splitTriangle.normals[2][2] = nNorm[2];
    flatBottom(splitTriangle, s);
    splitTriangle.X[0] = bot[0];
    splitTriangle.Y[0] = bot[1];
    splitTriangle.Z[0] = bot[2];
    splitTriangle.colors[0][0] = bot[3];
    splitTriangle.colors[0][1] = bot[4];
    splitTriangle.colors[0][2] = bot[5];
    splitTriangle.normals[0][0] = bNorm[0];
    splitTriangle.normals[0][1] = bNorm[1];
    splitTriangle.normals[0][2] = bNorm[2];
    flatTop(splitTriangle, s);
}

bool triangleType(Triangle t) {
	if (t.Y[0] != t.Y[1] && t.Y[0] != t.Y[2] && t.Y[1] != t.Y[2]) {
		return true;
	}
	return false;
}

bool lineCheck(Triangle t) {
	if (t.X[0] == t.X[1] && t.X[1] == t.X[2] && t.Y[0] == t.Y[1] && t.Y[1] == t.Y[2]) {
		return false;
	}
	if (t.X[0] == t.X[1] && t.Y[0] == t.Y[1]) {
		return false;
	}
	if (t.X[0] == t.X[2] && t.Y[0] == t.Y[2]) {
		return false;
	}
	if (t.X[1] == t.X[2] && t.Y[1] == t.Y[2]) {
		return false;
	}
	return triangleType(t);
}

bool valid(Triangle t, Screen s) {
    if (t.X[0] < 0 && t.X[1] < 0 && t.X[2] < 0) {
        return false;
    }
    if (t.X[0] > s.width && t.X[1] > s.width && t.X[2] > s.width) {
        return false;    
    }
    if (t.Y[0] < 0 && t.Y[1] < 0 && t.Y[2] < 0) {
        return false;
    }
    if (t.Y[0] > s.height && t.Y[1] > s.height && t.Y[2] > s.height) {
        return false;
    }
    return lineCheck(t);
}

bool isflatbottom(Triangle t) {
    if (t.Y[0] == t.Y[1] && t.Y[2] > t.Y[0]) {
        return true;
    }
    if (t.Y[0] == t.Y[2] && t.Y[1] > t.Y[0]) {
        return true;
    }
    if (t.Y[1] == t.Y[2] && t.Y[0] > t.Y[1]) {
        return true;
    }
    return false;
}

void rasterize(Triangle t, Screen s) {
    if (valid(t, s)) {
    	split(t, s);
    } else if (isflatbottom(t)) {
    	flatBottom(t, s);
    } else {
    	flatTop(t, s);
    }
}

int main() {
   vtkImageData *image = NewImage(1000, 1000);
   int npixels = 1000*1000;
   unsigned char *buffer = 
     (unsigned char *) image->GetScalarPointer(0,0,0);
   double *zbuffer = new double[npixels];
   for (int i = 0 ; i < npixels*3 ; ++i) {
       buffer[i] = 0;
   }
   for (int i = 0 ; i < npixels ; ++i) {
       zbuffer[i] = -1;
   }
   std::vector<Triangle> triangles = GetTriangles();   
   Screen screen;
   screen.buffer = buffer;
   screen.zbuffer = zbuffer;
   screen.width = 1000;
   screen.height = 1000;
   for (int i = 0; i < triangles.size(); ++i) {
   	   rasterize(triangles[i], screen);
   }
   WriteImage(image, "project1E");
}