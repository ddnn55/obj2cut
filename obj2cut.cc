#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <limits>
#include <unordered_map> // added for edge label mapping

/**************
 * command line tool to convert an OBJ file to an SVG file with
 * the OBJ's triangles packed into a 2D space.
 */

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

// Add a custom hash for std::pair<int,int>
struct PairHash {
    std::size_t operator()(const std::pair<int,int>& p) const {
        return std::hash<int>()(p.first) ^ (std::hash<int>()(p.second) << 1);
    }
};

int main(int argc, char* argv[]){
    if(argc < 2){
        std::cerr << "Usage: " << argv[0] << " path/to/file.obj [scale]\n";
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
    
    // Build edge-to-label mapping
    std::unordered_map<std::pair<int,int>, int, PairHash> edgeLabels;
    int labelCounter = 0;
    for(const Face &face : faces){
        for (int i = 0; i < 3; i++) {
            int j = (i+1)%3;
            int a = face.v[i], b = face.v[j];
            if(a > b) std::swap(a, b);
            std::pair<int,int> key = {a, b};
            if(edgeLabels.find(key) == edgeLabels.end()){
                edgeLabels[key] = labelCounter++;
            }
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
    
    // Compute the global bounding box minimum.
    double minX = std::numeric_limits<double>::max();
    double minY = std::numeric_limits<double>::max();
    for (const Triangle2D& tri : triangles2D) {
        minX = std::min(minX, std::min(tri.p0.x, std::min(tri.p1.x, tri.p2.x)));
        minY = std::min(minY, std::min(tri.p0.y, std::min(tri.p1.y, tri.p2.y)));
    }

    // Translate all triangles so that the lower-left point is at (0,0)
    for (Triangle2D& tri : triangles2D) {
        tri.p0.x -= minX; tri.p0.y -= minY;
        tri.p1.x -= minX; tri.p1.y -= minY;
        tri.p2.x -= minX; tri.p2.y -= minY;
    }
    
    // Pack triangles using a greedy row packing algorithm.
    double margin = 20.0;
    double maxRowWidth = 1000.0;
    double currentX = margin, currentY = margin;
    double rowHeight = 0.0;
    for (auto& tri : triangles2D) {
        // Compute triangle bounding box.
        double tMinX = std::min(tri.p0.x, std::min(tri.p1.x, tri.p2.x));
        double tMaxX = std::max(tri.p0.x, std::max(tri.p1.x, tri.p2.x));
        double tMinY = std::min(tri.p0.y, std::min(tri.p1.y, tri.p2.y));
        double tMaxY = std::max(tri.p0.y, std::max(tri.p1.y, tri.p2.y));
        double width = tMaxX - tMinX;
        double height = tMaxY - tMinY;
        
        // If current row is full, move to next row.
        if (currentX + width > maxRowWidth - margin) {
            currentX = margin;
            currentY += rowHeight + margin;
            rowHeight = 0.0;
        }
        
        // Compute translation for this triangle.
        double offsetX = currentX - tMinX;
        double offsetY = currentY - tMinY;
        tri.p0.x += offsetX; tri.p0.y += offsetY;
        tri.p1.x += offsetX; tri.p1.y += offsetY;
        tri.p2.x += offsetX; tri.p2.y += offsetY;
        
        currentX += width + margin;
        rowHeight = std::max(rowHeight, height);
    }
    double svgWidth = maxRowWidth;
    double svgHeight = currentY + rowHeight + margin;

// Apply scaling factor if provided as second command-line argument.
    double scale = 10.0;
    if(argc >= 3) {
        scale = std::atof(argv[2]);
    }
    for (auto& tri : triangles2D) {
        tri.p0.x *= scale; tri.p0.y *= scale;
        tri.p1.x *= scale; tri.p1.y *= scale;
        tri.p2.x *= scale; tri.p2.y *= scale;
    }
    svgWidth *= scale;
    svgHeight *= scale;

    // Start SVG output
    std::cout << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << "\n";
    std::cout << "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"" 
              << svgWidth << "\" height=\"" << svgHeight << "\">\n";
    
    // Output each packed triangle.
    for (const Triangle2D &tri : triangles2D) {
        auto ptToStr = [&](const Vec2 &p) {
            std::ostringstream oss;
            oss << std::fixed << std::setprecision(2) << p.x << "," << p.y;
            return oss.str();
        };
        std::string points = ptToStr(tri.p0) + " " + ptToStr(tri.p1) + " " + ptToStr(tri.p2);
        std::cout << "  <polygon points=\"" << points << "\" stroke=\"black\" fill=\"none\" />\n";
    }
    
    // For each triangle, add text labels for each edge.
    for (size_t i = 0; i < triangles2D.size(); i++) {
        const Face &face = faces[i];
        const Triangle2D &tri = triangles2D[i];
        // Define lambda to output text element.
        auto printLabel = [&](const Vec2 &p1, const Vec2 &p2, int edgeLabel) {
            double midX = (p1.x + p2.x) / 2.0;
            double midY = (p1.y + p2.y) / 2.0;
            std::cout << "  <text x=\"" << midX << "\" y=\"" << midY
                      << "\" font-size=\"20\" fill=\"red\" text-anchor=\"middle\" dominant-baseline=\"middle\">"
                      << edgeLabel << "</text>\n";
        };
        // For each edge in the triangle.
        for (int j = 0; j < 3; j++) {
            int k = (j+1)%3;
            int a = face.v[j], b = face.v[k];
            if(a > b) std::swap(a,b);
            std::pair<int,int> key = {a, b};
            int edgeLabel = edgeLabels[key];
            if(j==0)
                printLabel(tri.p0, tri.p1, edgeLabel);
            else if(j==1)
                printLabel(tri.p1, tri.p2, edgeLabel);
            else
                printLabel(tri.p2, tri.p0, edgeLabel);
        }
    }
    
    std::cout << "</svg>\n";
    return 0;
}