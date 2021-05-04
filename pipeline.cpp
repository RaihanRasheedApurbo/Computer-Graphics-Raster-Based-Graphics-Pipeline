#include <bits/stdc++.h>
#include "bmp_image_codes/bitmap_image.hpp"
using namespace std;

class Vector
{
public:
    double x,y,z;
    friend ostream & operator << (ostream &out, const Vector &c);

    friend istream & operator >> (istream &in,  Vector &c);

    static double dotProduct(const Vector &x, const Vector &y)
    {
        return x.x*y.x + x.y*y.y + x.z*y.z;
    }

    static Vector crossProduct(const Vector &x, const Vector &y)
    {
        // (a2b3−a3b2)i−(a1b3−a3b1)j+(a1b2−a2b1)k.
        Vector result;
        result.x = x.y*y.z-x.z*y.y;
        result.y = -(x.x*y.z-x.z*y.x);
        result.z = x.x*y.y-x.y*y.x;
        return result;

    }

    static Vector normalize(const Vector& a)
    {
        Vector result;
        double divisor = sqrt(a.x*a.x+a.y*a.y+a.z*a.z);
        result.x = a.x/divisor;
        result.y = a.y/divisor;
        result.z = a.z/divisor;
        
        if(abs(result.x) < 0.0001)
        {
            result.x = 0;
        }
        if(abs(result.y) < 0.0001)
        {
            result.y = 0;
        }
        if(abs(result.z) < 0.0001)
        {
            result.z = 0;
        }


        return result;

    }

    static Vector rodrigezFormula(const Vector &x, const Vector &a, const double angle)
    {
        Vector normalizedA = Vector::normalize(a);
        Vector result;
        double constant1 = 0; // calculateion of (1-cos(theta)(a dot x))
        const double PIE = 2*acos(0);
        double angleInRedian = (PIE/180)*angle;
        constant1 = (1-cos(angleInRedian))*Vector::dotProduct(normalizedA,x);
        Vector aCrossX = Vector::crossProduct(normalizedA,x);
        double cosTheta = cos(angleInRedian);
        double sinTheta = sin(angleInRedian);

        result.x = cosTheta*x.x + constant1*normalizedA.x + sinTheta* aCrossX.x;
        result.y = cosTheta*x.y + constant1*normalizedA.y + sinTheta* aCrossX.y;
        result.z = cosTheta*x.z + constant1*normalizedA.z + sinTheta* aCrossX.z;

        if(abs(result.x) < 0.0001)
        {
            result.x = 0;
        }
        if(abs(result.y) < 0.0001)
        {
            result.y = 0;
        }
        if(abs(result.z) < 0.0001)
        {
            result.z = 0;
        }


        return result;
    }


};

ostream & operator << (ostream &out, const Vector &c)
{
    out<<c.x<<","<<c.y<<","<<c.z<<endl;
    // out << c.real;
    // out << "+i" << c.imag << endl;
    // return out;
    return out;
}
  
istream & operator >> (istream &in,  Vector &c)
{
    in>>c.x;
    in>>c.y;
    in>>c.z;
    // cout << "Enter Real Part ";
    // in >> c.real;
    // cout << "Enter Imaginary Part ";
    // in >> c.imag;
    // return in;
    return in;

}

class Triangle
{
public:
    Vector a,b,c;
    double topScanLine(const double topY, const double dy)
    {
        double maxY = max(a.y,max(b.y,c.y));
        if(topY>maxY)
        {
            double dif = topY-maxY; // .9 - .35 = .55
            int distance = dif/dy;  // .55/.2 = 2.75(double) = 2 int
            return topY-distance*dy; // .9 - .2*2 = .5; .35 is between .3 and .5 hence .5 is upper boundary
        }
        else
        {
            return topY;
        }
    }

    double bottomScanLine(const double bottomY, const double dy)
    {
        double minY = min(a.y,min(b.y,c.y));
        if(bottomY<minY)
        {
            double dif = minY-bottomY; //  - .35 - (-.9) = .55
            int distance = dif/dy;  // .55/.2 = 2.75(double) = 2 int
            return bottomY+distance*dy; // -.9 + .2*2 = .5; -.35 is between -.3 and -.5 hence .5 is upper boundary
        }
        else
        {
            return bottomY;
        }
    }
};

class Matrix
{
public:
    double table[4][4];

    Vector multiplyVector(const Vector &v)
    {
        double homogenousPoint[4];
        homogenousPoint[0] = v.x;
        homogenousPoint[1] = v.y;
        homogenousPoint[2] = v.z;
        homogenousPoint[3] = 1;

        return multiplyPoint(homogenousPoint);
           
    }
    
    Vector multiplyPoint(const double* point)
    {
        double result[4];
        for(int i=0;i<4;i++)
        {
            result[i] = 0;
            for(int j=0;j<4;j++)
            {
                result[i] += table[i][j] * point[j];
            }
        }
        Vector res;
        res.x = result[0]/result[3];
        res.y = result[1]/result[3];
        res.z = result[2]/result[3];
        return res;

    }

    

    Triangle multiplyTriangle(const Triangle &t)
    {
        Triangle result;
        result.a = multiplyVector(t.a);
        result.b = multiplyVector(t.b);
        result.c = multiplyVector(t.c);
        return result;
    }

    vector<Triangle> multiplyTriangles(const vector<Triangle> &triangles)
    {
        vector<Triangle> result;
        for(int i=0;i<triangles.size();i++)
        {
            Triangle t = multiplyTriangle(triangles[i]);
            result.push_back(t);
        }
        return result;
    }

    static Matrix multiply(const Matrix& matrix1, const Matrix& matrix2)
    {
        Matrix result;
        for(int i=0;i<4;i++)
        {
            for(int j=0;j<4;j++)
            {
                result.table[i][j] = 0;
                for(int k=0;k<4;k++)
                {
                    result.table[i][j] += matrix1.table[i][k]*matrix2.table[k][j];
                }
            }
        }
        return result;
    }
    static void makeIdentity(Matrix& matrix)
    {
        for(int i=0;i<4;i++)
        {
            for(int j=0;j<4;j++)
            {
                matrix.table[i][j] = 0;
                if(i==j)
                {
                    matrix.table[i][j] = 1;
                }
            }
        }
    }
};

class TransformationMachine
{
public:
    stack<Matrix> matrixStack;
    TransformationMachine()
    {
        Matrix t;
        Matrix::makeIdentity(t);
        matrixStack.push(t);
    }
    Vector transform(const Vector& point)
    {
        Matrix &t = matrixStack.top();
        double homogeneousPoint[4];
        homogeneousPoint[0] = point.x;
        homogeneousPoint[1] = point.y;
        homogeneousPoint[2] = point.z;
        homogeneousPoint[3] = 1;
        return t.multiplyPoint(homogeneousPoint);
    }

    Triangle transform(const Triangle& triangle)
    {
        Triangle transformedTriangle;
        transformedTriangle.a = transform(triangle.a);
        transformedTriangle.b = transform(triangle.b);
        transformedTriangle.c = transform(triangle.c);
        return transformedTriangle;
    }
    void add(Matrix &matrix)
    {
        Matrix &t = matrixStack.top();
        Matrix afterMultiplication = Matrix::multiply(t,matrix);
        pop();
        matrixStack.push(afterMultiplication);

    }
    void push()
    {
        if(matrixStack.size()>1)
        {
            matrixStack.push(matrixStack.top());
        }
    }
    void pop()
    {
        if(matrixStack.size()>1)
        {
            matrixStack.pop();
        }
    }


    
};

void printStringAscii(const string &s)
{
    for(int i=0;i<s.size();i++)
    {
        cout<<int(s[i])<<" ";
    }
    cout<<endl;
}

double toRedian(double angle)
{
    double PIE = 2 * acos(0);
    return (PIE/180)*angle;
}

class Color
{
public:
    int r,g,b;

};

int main()
{
    
    
    
    ifstream inputFile("Test Cases (Updated)/4/scene.txt");
    ofstream outputFile("output1.txt");

    Vector eye;
    inputFile>>eye;
    // cout<<eye;

    Vector look;
    inputFile>>look;
    // cout<<look;

    Vector up;
    inputFile>>up;
    // cout<<up;
    double fovY,aspectRatio,near,far;
    inputFile>>fovY>>aspectRatio>>near>>far;
    // cout<<fovy<<endl;
    // cout<<aspectRatio<<endl;
    // cout<<near<<endl;
    // cout<<far<<endl;

    // initializing matrix stack
    // stage 1: Modeling Transformation
    TransformationMachine machine;
    vector<Triangle> triangles; 

    string s;
    getline(inputFile,s);// after getting all input before first command

    int counter = 0;
    while(true)
    {
        getline(inputFile,s);
        // cout<<s<<endl;
        // printStringAscii(s);
        if(s.find("end") != std::string::npos)
        {
            cout<<"total command: "<<counter<<endl;
            cout<<"exiting from loop"<<endl;
            break;
        }
        else if(s.find("triangle") != std::string::npos)
        {
            Triangle t1;
            inputFile>>t1.a>>t1.b>>t1.c;
            Triangle afterTransform = machine.transform(t1);
            outputFile<<afterTransform.a<<afterTransform.b<<afterTransform.c<<endl;
            triangles.push_back(afterTransform);
            counter++;
        }
        else if(s.find("translate") != std::string::npos)
        {
            Vector p1;
            inputFile>>p1;
            Matrix matrix;
            Matrix::makeIdentity(matrix);
            matrix.table[0][3] = p1.x;
            matrix.table[1][3] = p1.y;
            matrix.table[2][3] = p1.z;
            machine.add(matrix);
            counter++;
        }
        else if(s.find("scale") != std::string::npos)
        {
            Vector p1;
            inputFile>>p1;
            Matrix matrix;
            Matrix::makeIdentity(matrix);
            matrix.table[0][0] = p1.x;
            matrix.table[1][1] = p1.y;
            matrix.table[2][2] = p1.z;
            machine.add(matrix);
            counter++;
        }
        else if(s.find("rotate") != std::string::npos)
        {
            double angle;
            Vector p1;
            inputFile>>angle>>p1;
            // c1=R(i,a,angle)
            // c2=R(j,a,angle)
            // c3=R(k,a,angle)
            Vector t;
            
            t.x = 1, t.y = 0, t.z = 0; // x = i vector
            Vector c1 = Vector::rodrigezFormula(t,p1,angle);

            t.x = 0, t.y = 1, t.z = 0; // x = j vector
            Vector c2 = Vector::rodrigezFormula(t,p1,angle);

            t.x = 0, t.y = 0, t.z = 1; // x = k vector
            Vector c3 = Vector::rodrigezFormula(t,p1,angle);

            Matrix matrix;
            Matrix::makeIdentity(matrix);
            matrix.table[0][0] = c1.x;
            matrix.table[1][0] = c1.y;
            matrix.table[2][0] = c1.z;
            matrix.table[0][1] = c2.x;
            matrix.table[1][1] = c2.y;
            matrix.table[2][1] = c2.z;
            matrix.table[0][2] = c3.x;
            matrix.table[1][2] = c3.y;
            matrix.table[2][2] = c3.z;

            machine.add(matrix);

            counter++;
        }
        else if(s.find("push") != std::string::npos)
        {
            machine.push();
            counter++;
        }
        else if(s.find("pop") != std::string::npos)
        {
            machine.pop();
            counter++;
        }
         
    }
    

    
    inputFile.close();
    outputFile.close();

    // stage 1: Modeling Transformation completed

    // stage 2: View Transformation
    Vector l,r,u;
    l.x = look.x-eye.x;
    l.y = look.y-eye.y;
    l.z = look.z-eye.z;
    l = Vector::normalize(l);
    r = Vector::crossProduct(l,up);
    r = Vector::normalize(r);
    u = Vector::crossProduct(r,l);
    u = Vector::normalize(u);

    Matrix trans;
    Matrix::makeIdentity(trans);
    trans.table[0][3] = -eye.x;
    trans.table[1][3] = -eye.y;
    trans.table[2][3] = -eye.z;

    Matrix rot;
    Matrix::makeIdentity(rot);
    rot.table[0][0] = r.x;
    rot.table[0][1] = r.y;
    rot.table[0][2] = r.z;
    
    rot.table[1][0] = u.x;
    rot.table[1][1] = u.y;
    rot.table[1][2] = u.z;

    rot.table[2][0] = -l.x;
    rot.table[2][1] = -l.y;
    rot.table[2][2] = -l.z;

    Matrix v = Matrix::multiply(rot,trans);

    vector<Triangle> trianglesStage2 = v.multiplyTriangles(triangles);

    ofstream outputFile2("output2.txt");

    for(int i=0;i<trianglesStage2.size();i++)
    {
        outputFile2<<trianglesStage2[i].a<<trianglesStage2[i].b<<trianglesStage2[i].c<<endl;
    }

    outputFile2.close();

    // stage 2: View Transformation complted

    //stage 3: Projection Transformation

    // double fovy,aspectRatio,near,far;
    // inputFile>>fovy>>aspectRatio>>near>>far;
    double fovX = fovY * aspectRatio;
    double tConstant = near * tan(toRedian(fovY/2));
    double rConstant = near * tan(toRedian(fovX/2));
    Matrix p;
    Matrix::makeIdentity(p);
    p.table[0][0] = near/rConstant;
    p.table[1][1] = near/tConstant;
    p.table[2][2] = -(far+near)/(far-near);
    p.table[2][3] = -(2*far*near)/(far-near);
    p.table[3][2] = -1;
    p.table[3][3] = 0;


    vector<Triangle> trianglesStage3 = p.multiplyTriangles(trianglesStage2);

    ofstream outputFile3("output3.txt");

    for(int i=0;i<trianglesStage3.size();i++)
    {
        outputFile3<<trianglesStage3[i].a<<trianglesStage3[i].b<<trianglesStage3[i].c<<endl;
    }

    outputFile3.close();

    ifstream inputFile2("Test Cases (Updated)/4/config.txt");

    double screenWidth, screenHeight, leftLimit, rightLimit, bottomLimit,
           topLimit, zFront, zNear;
        
    inputFile2>>screenWidth>>screenHeight;
    inputFile2>>leftLimit>>bottomLimit;
    rightLimit = -leftLimit;
    topLimit = - bottomLimit;
    inputFile2>>zFront>>zNear;

    vector<Color> triangleColors;
    for(int i=0;i<trianglesStage3.size();i++)
    {
        int RGB_LIMIT = 255;
        Color t;
        t.r = rand()%RGB_LIMIT;
        t.g = rand()%RGB_LIMIT;
        t.b = rand()%RGB_LIMIT;
        triangleColors.push_back(t);
    }
    cout<<"kill meh"<<endl;
    
    double dx = (rightLimit-leftLimit)/screenWidth;
    double dy = (topLimit-bottomLimit)/screenHeight;
    double topY = topLimit - dy/2;
    double leftX = leftLimit + dx/2;
    

    vector<vector<double>> zBuffer;
    for(int i=0;i<screenHeight;i++)
    {
        vector<double> t;
        for(int j=0;j<screenWidth;j++)
        {
            t.push_back(zNear);
        }
        zBuffer.push_back(t);
    }

    bitmap_image image(screenWidth,screenHeight);

    for(int i=0;i<trianglesStage3.size();i++)
    {
        Triangle &t = trianglesStage3[i];
        double topScanLine = t.topScanLine(topY,dy);
        double bottomScanLine = t.bottomScanLine(-topY,dy);
        int rowStart = (topY-topScanLine)/dy;
        int rowEnd = (topY-bottomScanLine)/dy;
        for(int j=rowStart;j<=rowEnd;j++)
        {
            
            double ys = topY-j*dy;
            Vector p1 = t.a;
            Vector p2 = t.b;
            Vector p3 = t.c;
            if(p1.y==p2.y)
            {
                Vector t = p1;
                p1 = p3;
                p3 = t;
            }
            if(p1.y==p3.y)
            {
                Vector t = p1;
                p1 = p2;
                p2 = t;
            }
            double x1 = p1.x, y1 =p1.y, z1 = p1.z;
            double x2 = p2.x, y2 =p2.y, z2 = p2.z;
            double x3 = p3.x, y3 =p3.y, z3 = p3.z;
            // double y1 = t.a.y, y2 = t.b.y, y3 = t.c.y;
            // double z1 = t.a.z, z2 = t.b.z, z3 = t.c.z;
            // double x1 = t.a.x, x2 = t.b.x, x3 = t.c.x;
            double za = z1 - (z1-z2) * ((y1-ys)/(y1-y2));
            double zb = z1 - (z1-z3) * ((y1-ys)/(y1-y3));
            double xa = x1 - (x1-x2) * ((y1-ys)/(y1-y2));
            double xb = x1 - (x1-x3) * ((y1-ys)/(y1-y3));

            int colStart = (xa-leftX)/dx;
            int colEnd = (xb-leftX)/dx;
            for(int k=colStart;k<=colEnd;k++)
            {
                double xp = leftX + k*dx;
                double zp = zb - (zb-za) * ((xb-xp)/(xb-xa));
                // cout<<zp<<" "<<xb<<" "<<xa<<endl;
                if(zFront<=zp && zp<=zNear)
                {
                    if(zp<zBuffer[j][k])
                    {
                        zBuffer[j][k] = zp;
                        Color &t = triangleColors[i];
                        image.set_pixel(k,j,t.r,t.b,t.g);
                    }
                }
            }
        }
    }

    ofstream outputFile4("z_buffer.txt");
    outputFile4.precision(6);
    for(int i=0;i<screenHeight;i++)
    {
        for(int j=0;j<screenWidth;j++)
        {
            if(zBuffer[i][j]==2)
            {
                outputFile4<<setw(8)<<setw(4);
            }
            else
            {
                outputFile4<<zBuffer[i][j]<<setw(4);
            }
        }
        outputFile4<<endl;
    
    }
    outputFile4.close();

    image.save_image("test.bmp");


}