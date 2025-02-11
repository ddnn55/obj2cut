#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>
#include <cstdlib>
#include <iomanip>

struct Vec3 {
    double x, y, z;
    Vec3() : x(0), y(0), z(0) {}
    Vec3(double a, double b, double c): x(a), y(b), z(c) {}
};

Vec3 operator-(const Vec3& a, const Vec3& b) {
    return Vec3(a.x-b.x, a.y-b.y, a.z-b.z);
}

Vec3 operator+(const Vec3& a, const Vec3& b) {
    return Vec3(a.x+b.x, a.y+b.y, a.z+b.z);
}

Vec3 operator*(const Vec3& a, double s) {
    return Vec3(a.x*s, a.y*s, a.z*s);
}

double dot(const Vec3 &a, const Vec3 &b){
    return a.x*b.x + a.y*b.y + a.z*b.z;
}

Vec3 cross(const Vec3 &a, const Vec3 &b){
    return Vec3(a.y*b.z - a.z*b.y,
                a.z*b.x - a.x*b.z,
                a.x*b.y - a.y*b.x);
}

Vec3 normalize(const Vec3 &a){
    double len = std::sqrt(dot(a,a));
    if(len == 0) return a;
    return a * (1.0 / len);
}

struct Vec2 {
    double x, y;
    Vec2(): x(0), y(0) {}
    Vec2(double a, double b): x(a), y(b){}
};

struct Triangle2D {
    Vec2 p0, p1, p2;
};

struct Face {
    int v[3]; // indices to vertices (0-indexed)
};

int main(int argc, char* argv[]){
    if(argc < 2){
        std::cerr << "Usage: " << argv[0] << " path/to/file.obj\n";
        return 1;
    }
    
    std::ifstream ifs(argv[1]);
    if(!ifs){
        std::cerr << "Error opening file " << argv[1] << "\n";
        return 1;
    }
    
    std::vector<Vec3> vertices;
    std::vector<Face> faces;
    std::string line;
    while(std::getline(ifs, line)){
        if(line.size()<2) continue;
        std::istringstream iss(line);
        std::string type;
        iss >> type;
        if(type == "v"){
            double x, y, z;
            iss >> x >> y >> z;
            vertices.push_back(Vec3(x, y, z));
        } else if(type == "f"){
            Face face;
            int idx[3];
            // OBJ indices are 1-indexed; we only take first part if slashes are present.
            std::string token;
            for (int i = 0; i < 3 && iss >> token; ++i) {
                std::istringstream tokenStream(token);
                std::string indexStr;
                std::getline(tokenStream, indexStr, '/');
                idx[i] = std::stoi(indexStr);
            }
            // store as 0-indexed
            face.v[0] = idx[0] - 1;
            face.v[1] = idx[1] - 1;
            face.v[2] = idx[2] - 1;
            faces.push_back(face);
        }
    }
    
    // For each face, project the triangle in 2D using its plane.
    std::vector<Triangle2D> triangles2D;
    for(auto& f : faces){
        Vec3 v0 = vertices[f.v[0]];
        Vec3 v1 = vertices[f.v[1]];
        Vec3 v2 = vertices[f.v[2]];
        
        Vec3 edge1 = v1 - v0;
        Vec3 edge2 = v2 - v0;
        Vec3 normal = normalize(cross(edge1, edge2));
        // Create in-plane basis:
        Vec3 u = normalize(edge1);
        Vec3 v = cross(normal, u); // Ensures u,v are in the triangle's plane.
        
        Triangle2D tri;
        Vec3 d0 = v0 - v0; // zero
        Vec3 d1 = v1 - v0;
        Vec3 d2 = v2 - v0;
        tri.p0 = Vec2(dot(d0, u), dot(d0, v));
        tri.p1 = Vec2(dot(d1, u), dot(d1, v));
        tri.p2 = Vec2(dot(d2, u), dot(d2, v));
        triangles2D.push_back(tri);
    }
    
    // Determine grid dimensions
    int count = triangles2D.size();
    int gridCols = std::ceil(std::sqrt(count));
    int gridRows = std::ceil((double)count / gridCols);
    double cellWidth = 200.0;
    double cellHeight = 200.0;
    
    // Start SVG output
    double svgWidth = gridCols * cellWidth;
    double svgHeight = gridRows * cellHeight;
    
    std::cout << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << "\n";
    std::cout << "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"" 
              << svgWidth << "\" height=\"" << svgHeight << "\">\n";
    
    // For each triangle, center it in its cell
    for (size_t i = 0; i < triangles2D.size(); ++i) {
        int col = i % gridCols;
        int row = i / gridCols;
        double offsetX = col * cellWidth;
        double offsetY = row * cellHeight;
        
        // Get triangle bounding box
        Triangle2D tri = triangles2D[i];
        double minX = std::min(tri.p0.x, std::min(tri.p1.x, tri.p2.x));
        double maxX = std::max(tri.p0.x, std::max(tri.p1.x, tri.p2.x));
        double minY = std::min(tri.p0.y, std::min(tri.p1.y, tri.p2.y));
        double maxY = std::max(tri.p0.y, std::max(tri.p1.y, tri.p2.y));
        double triWidth = maxX - minX;
        double triHeight = maxY - minY;
        
        // Compute translation to center the triangle in the cell.
        double translateX = offsetX + (cellWidth - triWidth)/2 - minX;
        double translateY = offsetY + (cellHeight - triHeight)/2 - minY;
        
        // Define points after translation.
        auto ptToStr = [&](const Vec2 &p) {
            std::ostringstream oss;
            oss << std::fixed << std::setprecision(2) 
                << (p.x + translateX) << "," << (p.y + translateY);
            return oss.str();
        };
        
        std::string points = ptToStr(tri.p0) + " " + ptToStr(tri.p1) + " " + ptToStr(tri.p2);
        std::cout << "  <polygon points=\"" << points << "\" stroke=\"black\" fill=\"none\" />\n";
    }
    
    std::cout << "</svg>\n";
    return 0;
}