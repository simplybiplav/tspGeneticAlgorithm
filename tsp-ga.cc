#include "tsp-ga.hh"
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <set>
//bool debug = true;
bool debug = false;
void printPointVector(const vector<Point>& points)
{
    int loop = 0;
    for ( vector<Point>::const_iterator iter = points.begin(); iter != points.end(); iter ++, loop++)
    {
        CoordinateType x , y , z;
        iter->GetAxes(x,y,z);
        cout << "Point " << (loop+1) << " =[" << x << " "<<y<<" " <<z<<" ]"<<endl;
    }

}
void printIntVector(const vector<int>& vertices)
{

    for ( vector<int>::const_iterator iter = vertices.begin(); iter != vertices.end(); iter ++)
    {
        cout << *iter << " ";
    }

}

Point::Point(CoordinateType x, CoordinateType y, CoordinateType z)
{
    xCoord = x;
    yCoord = y;
    zCoord = z;
}

const void Point::GetAxes(CoordinateType& x, CoordinateType& y, CoordinateType& z) const
{
   x  = xCoord;
   y  = yCoord;
   z  = zCoord;
}
TSPGenome::TSPGenome(int numPoints):value(-1)
{
    for(int loop=0; loop <numPoints; loop++)
    {
        sequence.push_back(loop);
    }
    random_shuffle(sequence.begin(),sequence.end());
}

TSPGenome::TSPGenome(const vector<int>& order):sequence(order),value(-1)
{
}

vector<int> TSPGenome:: getOrder() const
{
    return sequence;
}

void  TSPGenome::computeCircuitLength(const vector<Point>& points)
{

    double totalLength = 0;
    if (debug)
    {
        cout << __func__ << "order"<<endl;
        printIntVector(sequence);
        cout << __func__ << "order"<<endl;
    }
    for ( vector<int>::const_iterator iter = sequence.begin(); iter != sequence.end() ; )
    {
        double length =0 ;
        CoordinateType x1,y1,z1;
        CoordinateType x2,y2,z2;
        points[*iter].GetAxes(x1,y1,z1);
        iter++;
        if (iter == sequence.end() )
        {
            iter = sequence.begin();
            points[*iter].GetAxes(x2,y2,z2);
            iter = sequence.end();
        }  
        else
        {
            points[*iter].GetAxes(x2,y2,z2);
        }
        length = sqrt ( pow((x2 - x1),2) + pow ((y2 - y1), 2) + pow ( (z2 - z1), 2)); 
        totalLength += length;

        if (debug)
        {
            cout << "Point 1 [" << x1 << " "<< y1 <<" "<< z1 << "] "<< "Point 2 [" << x2 << " "<< y2 <<" "<< z2 <<"] "<<"distance" << length<<endl;
        } 
    }

    value = totalLength;

} 


double TSPGenome::getCircuitLength() const
{
    return value;
}

void TSPGenome::mutate()
{
    int index1 = rand() % sequence.size(); 
    int index2 = rand() % sequence.size(); 
    if (index1 == index2)
    {
        index2++;
        index2 = index2 >= sequence.size()?0 : index2 ;
    }
    int c = sequence[ index1 ];
    sequence [index1 ] = sequence [index2];
    sequence [index2 ] = c;
}

TSPGenome crosslink(const TSPGenome &g1, const TSPGenome &g2)
{
    vector<int> g1Order = g1.getOrder();
    vector<int> g2Order = g2.getOrder();
    vector<int> siblingOrder;
    set<int> OrderSet;
    int index = rand() % g1Order.size();
    if ( g1Order.size() != g2Order.size())
    {
        cout<< "Logic error, g1 and g2 size changed "<< g1Order.size() << " g2:"<< g2Order.size() <<endl;
    }
    for (int loop = 0; loop < index ; loop++)
    {
        if ( OrderSet.find(g1Order[loop]) != OrderSet.end() )
        {
            cout<< "Logic error, duplicate"<< g1Order[loop] <<endl;
        }
        siblingOrder.push_back( g1Order[loop]);
        OrderSet.insert( g1Order[loop]);
    }
    for (int loop = 0; loop < g2Order.size() ;loop++)
    {
        if ( OrderSet.find(g2Order[loop]) == OrderSet.end() )
        {
            siblingOrder.push_back( g2Order[loop]);
            OrderSet.insert( g2Order[loop]);
        }
    }
    if ( g1Order.size() > siblingOrder.size())
    {
        cout<< "Logic error, size changed"<< g1Order.size() << " sibling:"<< siblingOrder.size() <<endl;
    }
    TSPGenome sibling(siblingOrder);
    return sibling;
}

bool isShorterPath( const TSPGenome &g1, const TSPGenome & g2)
{
    return g1.getCircuitLength() <  g2.getCircuitLength();
}

TSPGenome findAShortPath(const vector<Point> &points,

  int populationSize, int numGenerations,
  int keepPopulation, int numMutations)
{
    vector<TSPGenome> tspPopulation;
    for (int loop =0; loop < populationSize; loop++)
    {
       TSPGenome genome(points.size());
       tspPopulation.push_back(genome); 
    }

    for (vector<TSPGenome>::iterator iter = tspPopulation.begin();
            iter != tspPopulation.end();
            iter++)
    {
        iter->computeCircuitLength(points);
    }

    for (int genLoop = 1; genLoop <= numGenerations; genLoop++)
    {
        sort(tspPopulation.begin(),tspPopulation.end(),isShorterPath);
        for (int sibLoop = (populationSize - keepPopulation)-1; sibLoop < populationSize; sibLoop++)
        {
           TSPGenome g1 = tspPopulation[rand() % keepPopulation];
           TSPGenome g2 = tspPopulation[ rand() % keepPopulation ];
           tspPopulation[sibLoop] = crosslink( g1,g2);
           tspPopulation[sibLoop].computeCircuitLength(points);
        }
        for(int mutateLoop = 0; mutateLoop < numGenerations; mutateLoop++)
        {
            // don't mutate the best solution so far i.e 0 index
            int index = 1 + rand() % (populationSize-1);
            tspPopulation[index].mutate();
            tspPopulation[index].computeCircuitLength(points);
        }
        if( genLoop % 10 == 0)
        {
            cout << "Generation "<<genLoop<<" : shortest path is "<<tspPopulation[0].getCircuitLength() << " ";
            printIntVector(tspPopulation[0].getOrder());
            cout <<endl;
        }
    }

    return tspPopulation[0];

}

int main( int argc, char**argv)
{

    int noOfPoints;
    CoordinateType x, y, z;
    vector<Point> VerticesVector;
    vector<int> shortestVectorPath;
    double shortestDistance = 0;
    int populationSize = 1000;
    int numGenerations = 100;
    double populationPercent = .2;
    double mutateRatio = 1.5;
    int keepPopulation = 6; int numMutations = 3;
    if (argc < 5 )
    {
        cout << "Invalid no args : Enter population generation population% mutation" <<endl; 
    }
    
    populationSize = atoi (argv[1]);
    numGenerations = atoi ( argv[2]);
    populationPercent = atof (argv[3]);
    mutateRatio = atof(argv[4]);
   
    cout << populationSize << "  population size "<<endl;
    cout << numGenerations << " no of gens"<<endl;
    cout << populationPercent << " population %" <<endl;
    cout << mutateRatio << " mutate ratio" <<endl;
    cout << "Enter the number of points" << endl;
    cin >> noOfPoints;
    for (int loop =0; loop < noOfPoints;loop++)
    {
        cout << "Enter Point[" << loop<<"]:"<<endl;  
        cin >> x;
        cin >> y;
        cin >> z;
        VerticesVector.push_back(Point(x,y,z));
    }
    keepPopulation = populationSize * populationPercent;
    numMutations = mutateRatio * populationSize;
    printPointVector(VerticesVector);
    TSPGenome shortGenome(findAShortPath(VerticesVector,populationSize,numGenerations,keepPopulation,numMutations));
    shortestVectorPath = shortGenome.getOrder();
    shortestDistance = shortGenome.getCircuitLength();
    cout << "Shortest Path is [ ";
    printIntVector(shortestVectorPath);
    cout << "] with distance = " <<  shortestDistance <<endl;

}

