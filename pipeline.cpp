#include <bits/stdc++.h>
using namespace std;

class Vector
{
public:
    double x,y,z;
    friend ostream & operator << (ostream &out, const Vector &c);
    friend istream & operator >> (istream &in,  Vector &c);


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

int main()
{
    

    ifstream inputFile("scene.txt");

    Vector eye;
    inputFile>>eye;
    cout<<eye;

    Vector look;
    inputFile>>look;
    cout<<look;

    Vector up;
    inputFile>>up;
    cout<<up;

    double fovy,aspectRatio,near,far;
    inputFile>>fovy>>aspectRatio>>near>>far;
    cout<<fovy<<endl;
    cout<<aspectRatio<<endl;
    cout<<near<<endl;
    cout<<far<<endl;


    
    inputFile.close();

}