#include<vector>
using namespace std;
typedef double CoordinateType;



class Point
{
private:
    CoordinateType xCoord;
    CoordinateType yCoord;
    CoordinateType zCoord;

public:
    Point( CoordinateType x, CoordinateType y, CoordinateType z);
    const void GetAxes(CoordinateType& x, CoordinateType& y, CoordinateType& z) const;
};
class TSPGenome
{
private:
    vector<int> sequence;
    double  value;
public:
    TSPGenome(int numPoints);
    TSPGenome(const vector<int>& order);
    vector<int> getOrder() const;
    void   computeCircuitLength(const vector<Point>& points); 
    double getCircuitLength() const;
    void mutate();
};

