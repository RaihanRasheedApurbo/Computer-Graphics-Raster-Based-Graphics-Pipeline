#include <bits/stdc++.h>
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
};

class Matrix
{
public:
    double table[4][4];

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
        res.x = result[0];
        res.y = result[1];
        res.z = result[2];
        return res;

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

int main()
{
    
    
    
    ifstream inputFile("scene.txt");
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
    double fovy,aspectRatio,near,far;
    inputFile>>fovy>>aspectRatio>>near>>far;
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
    ofstream outputFile("output2.txt");
    Vector l,r,u;
    l.x = look.x-eye.x;
    l.y = look.y-eye.y;
    l.z = look.z-eye.z;
    l = Vector::normalize(l);
    r = Vector::crossProduct(l,up);
    u = Vector::crossProduct(r,l);
    
    



}