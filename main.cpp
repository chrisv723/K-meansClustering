/*
 *  main.cpp
 *
 *  Created on: Mar 31, 2020
 *  Author: Christopher Valerio
 * 
 *  Given a list of 2-D points, and a cluster number, K, 
 *  the task is to partition the input point set to K clusters such that all points
 *  within one cluster are closer to the centroid of their own cluster (in distance) 
 *  than to the centroids of all other clusters.
 */

#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <cmath>
#include <ctime>
using namespace std;

class Kmean
{
public:
	class point
	{
	public:
		double xCoord;				// x component of given point
		double yCoord;				// y component of given point
		int label = 0;				// given points label. initialized to zero
		double distance = 99999.00; // distance to its own cluster centroid. initialized to 99999.00
	};

	int k;			 // total number of clusters, given in argv[2]
	int numPts;		 // total number of points
	point *pointset; // array of point class objects of size numPoints

	int numRows;
	int numCols;

	int **displayAry;	 // 2D array size of numRows by numCols. used for displaying purposes
	point *kCentroidAry; // Array of point class of size k, used to store information on centroids
	int change;			 // used to track label changes, initialized to zero

	/**
	 * parameterized constructor for objects of type Kmean
	 * 
	 * @param r - total number of rows 
	 * @param c - total number of columns
	 * @param numPts - total number of points
	 * @param k - number of clusters by which to group points
	 */
	Kmean(int r, int c, int numPts, int k)
	{
		this->k = k;
		this->numPts = numPts;
		pointset = new point[numPts];
		numRows = r;
		numCols = c;

		displayAry = new int *[numRows];
		for (int i = 0; i < numRows; i++)
		{
			displayAry[i] = new int[numCols];
		}

		kCentroidAry = new point[k + 1];
		change = 0;
	}

	/**
	 * 
	 * function reads in data from given input file and pointset array is loaded
	 * with data points
	 * 
	 *  @param in - reference to input file stream used to read in all given points from specified input file
	 *  @param pointset - Array of all points specified in file
	*/
	void loadPointSet(ifstream &in, point *pointset);

	/**
	 * 
	 *  performs k means clustering algorithm and outputs resulting information to specified output file
	 * 
	 *  @param ptSet - array of point class objects of size numPoints
	 *  @param k - number representing how many cluster centroids to construct
	 *  @param kCentroidAry - array used to store information about constructed centroids
	 *  @param out - output file stream used to display data
	*/
	void kMeansClustering(point *ptSet, int k, point *kCentroidAry, ofstream &out);

	/**
	 * 
	 * randomly select k number of points to be used as data centroids
	 * does not select identical centroids
	 * 
	 * @param ptSet - array of all points specified from given input file
	 * @param kCentroidAry - array of randomly chosen centroid points
	*/
	void selectKcentroids(point *ptSet, point *kCentroidAry);

	/**
	 * 
	 * This functions goes through the entire pointSet array only once to compute the centroids for all K clusters
	 * Stores the computed centroids in kCentroidAry
	 * 
	 * @param ptSet - array of all points specified from given input file
	 * @param kCentroidAry - array of randomly chosen centroid points
	*/
	void computeCentroids(point *ptSet, point *kCentroidAry);

	/**
	 * 
	 * this function computes the distance from a point pt to each of the K centroids.
	 * This methods returns minLabel 
	 * 
	 * @param pt - point from which to compute distance
	 * @param kCentroidAry - kCentroidAry holding each of the K centroids to which to compute distance
	 * @param minDist - calculated min distance
	 * 
	 * @return minimum label found
	 * 
	*/
	int distanceMinLabel(point pt, point *kCentroidAry, int minDist);

	/**
	 * 
	 * This function computes the distance from a point to a given centroid
	 * 
	 * @param pt - given point object containing an x & y component, label, and distance to its own cluster
	 * @param whichCentroid - given centroid of type point containing an x & y component, label, and distance to its own cluster
	 * 
	*/
	double computeDist(point pt, point whichCentroid);

	/**
	 * 
	 * for each point i, in pointSet, plot the pointSet[i].label onto the display Ary at location of pointSet[i]'s xCoord & yCoord
	 * 
	 * @param ptSet - array of all points given
	 * @param displayAry - 2D array used for displaying purposes
	 * 
	*/
	void plotDisplayAry(point *ptSet, int **displayAry);

	/**
	 * 
	 * prints results of k means clustering showing how all points were grouped together 
	 * 
	 * @param displayAry - array used to loop through and display k means clustering information
	 * @param out - reference to output file stream used to output results of k means clustering up to k times
	 * @param iteration - counter variable used to keep track of loop iteration
	 * 
	*/
	void prettyPrint(int **displayAry, ofstream &out, int iteration);

	/**
	 * 
	 * Output numRows numCols, numPts, and the x-coord, y-coord of points in
	 * pointSet, using the outFile2 file format given in the above output specs.
	 * 
	 * @param ptSet - array of all points given
	 * @param out - reference to output file stream used to output results of all of the given points
	 * 
	*/
	void printResult(point *ptSet, ofstream &out);
};

// controls program operation
int main(int argc, char **argv)
{
	ifstream inFile;
	inFile.open(argv[1]);
	int k = stoi(argv[2]);
	ofstream outFile1, outFile2;
	outFile1.open(argv[3]);
	outFile2.open(argv[4]);

	string rows, cols, numPts;
	inFile >> rows;
	inFile >> cols;
	inFile >> numPts;

	Kmean *image = new Kmean(stoi(rows), stoi(cols), stoi(numPts), k);
	//cout << rows << " " << cols << " " << numPts << " " << k << endl;

	image->loadPointSet(inFile, image->pointset);

	image->kMeansClustering(image->pointset, image->k, image->kCentroidAry, outFile1);
	image->printResult(image->pointset, outFile2);

	return 0;
}

void Kmean::loadPointSet(ifstream &in, point *pointset)
{
	int index = 0;
	string x, y;
	while (!in.eof())
	{
		in >> x;
		in >> y;
		pointset[index].xCoord = stod(x);
		pointset[index].yCoord = stod(y);
		pointset[index].label = 0; // meaning  no label yet
		pointset[index].distance = 99999.00;
		index++;
	}
}

void Kmean::kMeansClustering(point *ptSet, int k, point *kCentroidAry, ofstream &out)
{

	//ofstream outt; outt.open("test.txt");
	int changes = 0;
	this->selectKcentroids(ptSet, kCentroidAry);
	bool first = true;
	int iteration = 0;
	while (changes > 2 || first == true)
	{
		first = false;
		changes = 0;

		for (int i = 0; i < this->numPts; i++)
		{
			int minLabel = this->distanceMinLabel(ptSet[i], kCentroidAry, (int)ptSet->distance);

			if (ptSet[i].label != minLabel)
			{
				ptSet[i].label = minLabel;
				changes++;
			}
		}

		if (changes > 2)
		{
			this->computeCentroids(ptSet, kCentroidAry);
			//changes = 0;
		}

		plotDisplayAry(ptSet, displayAry);
		prettyPrint(this->displayAry, out, iteration++);
	}
}

void Kmean::plotDisplayAry(point *ptSet, int **displayAry)
{

	for (int i = 0; i < this->numPts; i++)
	{
		int row = (int)ptSet[i].xCoord;
		int col = (int)ptSet[i].yCoord;
		displayAry[row - 1][col - 1] = ptSet[i].label;
	}
}

int Kmean::distanceMinLabel(point pt, point *kCentroidAry, int minDist)
{
	int min = minDist;
	int minLabel = 0;
	int label = 1;
	double dist;
	while (label <= this->k)
	{
		point whichCentroid = kCentroidAry[label];
		dist = computeDist(pt, whichCentroid);

		if (dist <= min)
		{
			minLabel = label;
			min = dist;
		}
		label++;
	}
	return minLabel;
}
double Kmean::computeDist(point pt, point whichCentroid)
{
	double x1 = pt.xCoord;
	double y1 = pt.yCoord;
	double x2 = whichCentroid.xCoord;
	double y2 = whichCentroid.yCoord;
	//cout << x1 << " " << y1 << endl;
	//cout << x2 << " " << y2 << endl;

	return sqrt(pow((x2 - x1), 2.0) + pow((y2 - y1), 2.0));
}

bool checkRepeat(int index, int *valAry, int size)
{
	for (int i = 0; i < size; i++)
	{
		if (valAry[i] == index)
			return true;
	}
	return false;
}
void Kmean::selectKcentroids(point *ptSet, point *kCentroidAry)
{

	int cnt = 0;
	int *usedVals = new int[this->k];
	int kCnt = 0;
	srand(time(NULL));

	for (; kCnt < this->k; kCnt++)
	{
		int index = rand() % this->numPts;
		bool used = checkRepeat(index, usedVals, this->k);
		while (used == true)
		{
			index = rand() % this->numPts;
			used = checkRepeat(index, usedVals, this->k);
		}
		usedVals[kCnt] = index;
		kCentroidAry[kCnt].xCoord = ptSet[index].xCoord;
		kCentroidAry[kCnt].yCoord = ptSet[index].yCoord;
		kCentroidAry[kCnt].label = kCnt;
		kCentroidAry[kCnt].distance = 0.0;
		//cout << index << endl;
	}
}

void Kmean::computeCentroids(point *ptSet, point *kCentroidAry)
{
	double *sumX = new double[(this->k + 1)];
	double *sumY = new double[(this->k + 1)];
	int *totalPt = new int[(this->k + 1)];
	for (int i = 0; i < k + 1; i++)
	{
		sumX[i] = 0.0;
		sumY[i] = 0.0;
		totalPt[i] = 0;
	}

	int label;
	int index = 0;
	while (index < this->numPts)
	{
		label = ptSet[index].label;
		sumX[label] += ptSet[index].xCoord;
		sumY[label] += ptSet[index].yCoord;
		totalPt[label] += 1;

		index++;
	}
	label = 1;
	while (label <= this->k)
	{
		if (totalPt[label] > 0.0)
		{
			kCentroidAry[label].xCoord = (sumX[label] / totalPt[label]);
			kCentroidAry[label].yCoord = (sumY[label] / totalPt[label]);
		}
		label++;
	}
}

void Kmean::prettyPrint(int **displayAry, ofstream &out, int iteration)
{

	out << "******************  Result of iteration " << iteration << "  K=" << this->k << "  ******************" << endl;
	for (int i = 0; i < this->numRows; i++)
	{
		for (int j = 0; j < this->numCols; j++)
		{

			if (displayAry[i][j] > 0)
				out << displayAry[i][j] << " ";
			else
				out << "  ";
		}
		out << endl;
	}
}

void Kmean::printResult(point *ptSet, ofstream &out)
{
	out << this->numRows << " " << this->numCols << endl
		<< this->numPts << endl;
	for (int i = 0; i < this->numPts; i++)
	{
		out << ptSet[i].xCoord << " " << ptSet[i].yCoord << " " << ptSet[i].label << endl;
	}
}
