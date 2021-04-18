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
    
    // stage 1: Modeling Transformation
    
    ifstream inputFile("scene.txt");

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
    string s;
    getline(inputFile,s);// after getting all input before first command

    int counter = 0;
    while(true)
    {
        getline(inputFile,s);
        cout<<s<<endl;
        // printStringAscii(s);
        if(s.find("end") != std::string::npos)
        {

            cout<<"exiting from loop"<<endl;
            break;
        }
        else if(s.find("triangle") != std::string::npos)
        {
            counter++;
        }
        else if(s.find("translate") != std::string::npos)
        {
            counter++;
        }
        else if(s.find("scale") != std::string::npos)
        {
            counter++;
        }
        else if(s.find("rotate") != std::string::npos)
        {
            counter++;
        }
        else if(s.find("push") != std::string::npos)
        {
            counter++;
        }
        else if(s.find("pop") != std::string::npos)
        {
            counter++;
        }
        // else
        // {
        //     cout<<"Something is wrong in input!"<<endl;
        //     break;
        // }   
    }
    cout<<counter<<" hey"<<endl;


    
    inputFile.close();

}