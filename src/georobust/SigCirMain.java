package georobust;

import java.io.IOException;
import java.lang.reflect.Array;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.Random;
import java.util.function.IntPredicate;

import javax.management.openmbean.ArrayType;
import javax.print.attribute.Size2DSyntax;

public class SigCirMain {
	
	static long[] time_prune = new long[5];
	static long[] time_refine = new long[5]; // used for record Monte Carlo simulation time cost
	
	String dataset_path_200K = "/Users/tangxun/Desktop/200K.csv"; // dataset used for experiment
	String dataset_path_300K = "/Users/tangxun/Desktop/300K.csv";
	String dataset_path_400K = "/Users/tangxun/Desktop/400K.csv";
	String dataset_path_500K = "/Users/tangxun/Desktop/500K.csv";
	String dataset_path_600K = "/Users/tangxun/Desktop/600K.csv";
	String dataset_toy = "/Users/tangxun/Desktop/ICDM_2015/datasets/toy.csv";
	
	double theta = 100;
	int cellSize = 5;
	double xmax = 50;
	double ymax = 50;
	double xmin = 0;
	double ymin = 0;
	/*
	void meshgrid(int[][] columnsCircleGrid, int[][] rowsCircleGrid, int top_x, int top_y){
		// same as "meshgrid()" in Matlab
		for (int i = 0; i < top_y; i++){
			for (int j = 0; j < top_x; j++){
				columnsCircleGrid[i][j] = j;
				rowsCircleGrid[i][j] = i;
			}
		}
	}
	*/
	
	/*
	void meshgrid(int[][] columnsCircleGrid, int[][] rowsCircleGrid, int yCells, int xCells){
		// same as "meshgrid()" in Matlab
		for (int i = 0; i < xCells; i++){
			for (int j = 0; j < yCells; j++){
				columnsCircleGrid[i][j] = j;
				rowsCircleGrid[i][j] = i;
			}
		}
	}
	*/
	
	double lrupper(int ctot, double uc, double ub, double lc, double lb){
		// equivalent to "lrupper" function in Matlab	
		double logcin = Math.log(uc);
		double logbin = Math.log(lb);
		double logcout = Math.log(ctot - lc);
		double logbout = Math.log(ctot - ub);
		double lru;
		
		if( (uc > lb) && (lc >= ub) ) {
			lru = (uc * (logcin - logbin)) + ((ctot - uc) * (logcout - logbout));
		}
		else if ((uc > lb) && (lc < ub)){
			lru = (uc * (logcin - logbin)) + 1;
		}
		else{
			lru = 0;
		}
		return lru;
	}
	int unique(double[][] array, int nRows, int nColumns){
		// return the number of unique rows, first and second and third columns being same means not unique
		int[] index = new int[nRows];
		int current_index = 0;
		
		boolean flag = false; // "true" indicates duplicate happens
		
		for (int i = nRows - 1; i >= 0; i--){
			for (int j = i-1; j >= 0; j--){
				if ((array[i][0] == array[j][0]) && (array[i][1] == array[j][1]) && (array[i][2] == array[j][2])){
					flag = true;
					break;
				}
			}
			if (flag == false){ // i is unique
				index[i] = current_index;
				current_index++;
			}
			else {
				index[i] = -1; // to be deleted since its duplicate
				flag = false;
			}
		}
		//System.out.println("current_index = "+current_index);

		double[][] temp_array = new double[current_index][nColumns];
		
		for (int i = 0; i < nRows; i++){ // index is reversed before, not re-order them
			if (index[i] != -1){
				for (int j=0; j < nColumns; j++){
					temp_array[current_index-1-index[i]][j] = array[i][j];
				}
			}
			for (int k = 0; k < nColumns; k++){
				array[i][k] = 0; // clean the old 2d array
			}
		}
		// copy data from temp array to array 
		for (int i = 0; i < current_index; i++){
			for (int j = 0; j < nColumns; j++){
				array[i][j] = temp_array[i][j];
			}
		}
		return current_index;
	}
	
	void writePointList (double[] circleGridULR_first, CountGridCells[][] countGridCells, int[][] filledCircle, ArrayList<ArrayList<SSTD_Point>> prunedSets, int size_x, int size_y, int counter){
		// same as writePointList in "Matlab"
		
		for (int i = 0; i < size_x; i++){
			for (int j = 0; j < size_y; j++){
				if (((i-circleGridULR_first[0]) * (i-circleGridULR_first[0]) + (j-circleGridULR_first[1]) * (j-circleGridULR_first[1])) < (circleGridULR_first[2] * circleGridULR_first[2])){
					filledCircle[i][j] = 1; // cell [x, y] is in the circle
					for (int k = 0; k < countGridCells[i][j].points.size(); k++){
						prunedSets.get(counter).add(countGridCells[i][j].points.get(k));
					}
					countGridCells[i][j].points.clear();
					countGridCells[i][j].countGrid = 0;
				}
				else{
					filledCircle[i][j] = 0; // cell [x, y] is not in the circle
				}
			}
		}
	}
	
	void setdiff (ArrayList<ArrayList<SSTD_Point>> prunedSets, int counter){
		// same as "setdiff" in Matlab
		int size = prunedSets.get(counter).size(); // number of points in this cell
		for (int i = size-1; i >= 0; i--){
			Dataset.pointSet.remove(prunedSets.get(counter).get(i));
		}
	}
	
	double minCircle_SEC (ArrayList<SSTD_Point> points, ArrayList<SSTD_Point> tmpSet, SSTD_Point center, ArrayList<SSTD_Point> circlePoints) {
		double radius = 0;
		int n = points.size();
		ArrayList<SSTD_Point> firstNPoints = new ArrayList<SSTD_Point>();
		if (tmpSet.size() == 1){
			radius = 0;
			center.x = tmpSet.get(0).x;
			center.y = tmpSet.get(0).y;
		}
		else if (tmpSet.size() == 2){
			radius = Math.sqrt( (tmpSet.get(0).x - tmpSet.get(1).x)*(tmpSet.get(0).x - tmpSet.get(1).x) + (tmpSet.get(0).y - tmpSet.get(1).y)*(tmpSet.get(0).y - tmpSet.get(1).y) ) / 2;
			center.x = (tmpSet.get(0).x + tmpSet.get(1).x) / 2;
			center.y = (tmpSet.get(0).y + tmpSet.get(1).y) / 2;
		}
		else if (tmpSet.size() == 3){
			enc3(tmpSet.get(0), tmpSet.get(1), tmpSet.get(2), center);
			radius = Math.sqrt( (tmpSet.get(0).x - center.x)*(tmpSet.get(0).x - center.x) );
			for (int i = 0; i < tmpSet.size(); i++){
				circlePoints.add(tmpSet.get(i));
			}
			return radius;
		}
		for (int i = 0; i < n; i++){
			if ((points.get(i).x - center.x)*(points.get(i).x - center.x) + (points.get(i).y - center.y)*(points.get(i).y - center.y) > radius*radius){
				tmpSet.add(points.get(i));
				firstNPoints.clear();
				for (int j = 0; j < i; j++){
					firstNPoints.add(points.get(i));
				}	
				radius = minCircle_SEC(firstNPoints, tmpSet, center, circlePoints);
				System.out.println(i);
			}
		}
		return radius;
	}

	
	double minCircle (ArrayList<SSTD_Point> points, boolean hullflag, SSTD_Point center, ArrayList<SSTD_Point> circlePoints){
		// same as "minCircle" in Matlab code. returns the radius
		QuickHull myQuickHull = new QuickHull();
		double radius = -1;
		double tol;

		int n = points.size(); // find minimal bounding circle from n points
		//System.out.println("in minCircle n = " + n);
		
		if ((hullflag == true) && (n > 3)){
			points = myQuickHull.quickHull(points);
		}
		n = points.size(); // number of points on the convex hull
		//System.out.println("in minCircle n = " + n);

		if (n == 0){
			center.x = -1;
			center.y = -1;
			return -1; // no point
		}
		else if (n == 1){
			center.x = points.get(0).x;
			center.y = points.get(0).y;
			radius = 0;
			//System.out.println("n = 1, center: ("+center.x + ", " + center.y+ ")");
			return radius;
		}
		else if (n == 2) {
			center.x = ( points.get(0).x + points.get(1).x ) / 2;
			center.y = ( points.get(0).y + points.get(1).y ) / 2;
			radius = Math.sqrt( (center.x-points.get(0).x)*(center.x-points.get(0).x) + (center.y-points.get(0).y)*(center.y-points.get(0).y));
			//System.out.println("n = 2, center: ("+center.x + ", " + center.y+ ")");
			return radius;
		}
		else if (n == 3){
			enc3(points.get(0), points.get(1), points.get(2), center);
			radius = Math.sqrt( (center.x-points.get(0).x)*(center.x-points.get(0).x) + (center.y-points.get(0).y)*(center.y-points.get(0).y));
			circlePoints.clear();
			circlePoints.add(points.get(0));
			circlePoints.add(points.get(1));
			circlePoints.add(points.get(2));
			//System.out.println("n = 3, center: ("+center.x + ", " + center.y+ ")");

			return radius;
		}
		else{ // if n >= 4
			int[] aset = new int[3];
			aset[0] = 0;
			aset[1] = 1;
			aset[2] = 2;
			
			int[] iset = new int[n-3];
			for (int i = 3; i < n; i++){
				iset[i-3] = i;
			}		
			tol = computeTol(points);
				
			int[] s1 = new int[3];
			SSTD_Point c1 = new SSTD_Point();
			double r1;
		
			int size_old = 10;
			int[][] old_sets = new int[size_old][3];
			double[] old_rads = new double[size_old];
			double[][] old_centers = new double[size_old][2];
			for (int i = 0; i < size_old; i++){
				old_sets[i][0] = -1; // -1 here represents NaN in Matlab code
				old_sets[i][1] = -1;
				old_sets[i][2] = -1;
				old_rads[i] = Double.MAX_VALUE;
				old_centers[i][0] = -1;
				old_centers[i][1] = -1;
			}
			
			boolean flag = true;
			while (flag){
				int i;
				Arrays.sort(aset);
				for (i = 0; i < size_old; i++){
					if ( (aset[0] == old_sets[i][0]) && (aset[1] == old_sets[i][1]) && (aset[2] == old_sets[i][2]) ){
						break; 
					}				
				}

				if (i < size_old){	// old.sets contains aset
					center.x = old_centers[0][0];
					center.y = old_centers[0][1];
					radius = old_rads[0];
					flag = false;
					continue;
				}
				enc3(points.get(aset[0]), points.get(aset[1]), points.get(aset[2]), center);
				//System.out.println(center.x + ", "+center.y);
				radius = Math.sqrt( (center.x-points.get(aset[0]).x)*(center.x-points.get(aset[0]).x) + (center.y-points.get(aset[0]).y)*(center.y-points.get(aset[0]).y));
				circlePoints.clear();
				circlePoints.add(points.get(aset[0]));
				circlePoints.add(points.get(aset[1]));
				circlePoints.add(points.get(aset[2]));
				
				if (radius < old_rads[size_old-1]){
					Arrays.sort(aset);
					old_sets[size_old-1][0] = aset[0];
					old_sets[size_old-1][1] = aset[1];
					old_sets[size_old-1][2] = aset[2];
					old_rads[size_old-1] = radius;
					old_centers[size_old-1][0] = center.x;
					old_centers[size_old-1][1] = center.y;
					
					// sorting;
					for (int j = 0; j < size_old-1; j++){
						for (int k = j+1; k < size_old-1; k++){
							if (old_rads[j] > old_rads[k]){
								
								double swap = old_rads[k]; //swap
								old_rads[k] = old_rads[j];
								old_rads[j] = swap;
								
								swap = old_centers[k][0];
								old_centers[k][0] = old_centers[j][0];
								old_centers[j][0] = swap;
								swap = old_centers[k][1];
								old_centers[k][1] = old_centers[j][1];
								old_centers[j][1] = swap;
								
								int swap2 = old_sets[k][0];
								old_sets[k][0] = old_sets[j][0];
								old_sets[j][0] = swap2;
								swap2 = old_sets[k][1];
								old_sets[k][1] = old_sets[j][1];
								old_sets[j][1] = swap2;
								swap2 = old_sets[k][2];
								old_sets[k][2] = old_sets[j][2];
								old_sets[j][2] = swap2;
							}
						}
					}
				}
				double[] r = new double[n-3];
				for (int j = 0; j < n-3; j++){
					r[j] = Math.sqrt( (points.get(iset[j]).x - center.x)*(points.get(iset[j]).x - center.x) + (points.get(iset[j]).y - center.y)*(points.get(iset[j]).y - center.y) );
				}
				
				double rmax = r[0]; // find rmax and its index
				int k_index = 0;
				for (int k = 0; k < n-3; k++){
					if (r[k] > rmax){
						rmax = r[k];
						k_index = k;
					}
				}
				if ((rmax - radius) <= tol){
					flag = false;
				}
				else {
					s1[0] = aset[1];
					s1[1] = aset[2];
					s1[2] = iset[k_index];		
					enc3(points.get(s1[0]), points.get(s1[1]), points.get(s1[2]), c1);
					r1 = Math.sqrt( (c1.x-points.get(s1[0]).x)*(c1.x-points.get(s1[0]).x) + (c1.y-points.get(s1[0]).y)*(c1.y-points.get(s1[0]).y) );
					circlePoints.clear();
					circlePoints.add(points.get(s1[0]));
					circlePoints.add(points.get(s1[1]));
					circlePoints.add(points.get(s1[2]));
					if ( (c1.x-points.get(aset[0]).x)*(c1.x-points.get(aset[0]).x) + (c1.y-points.get(aset[0]).y)*(c1.y-points.get(aset[0]).y)
							<= r1*r1 ){  //compare these two norm^2
						center.x = c1.x;
						center.y = c1.y;
						//System.out.println(center.x + ", "+center.y);
						radius = r1;
						
						int swap = aset[0];
						aset[0] = iset[k_index];
						iset[k_index] = swap; 
						continue;
					}
					
					s1[0] = aset[0];
					s1[1] = aset[2];
					s1[2] = iset[k_index];
					enc3(points.get(s1[0]), points.get(s1[1]), points.get(s1[2]), c1);
					r1 = Math.sqrt( (c1.x-points.get(s1[0]).x)*(c1.x-points.get(s1[0]).x) + (c1.y-points.get(s1[0]).y)*(c1.y-points.get(s1[0]).y) );
					circlePoints.clear();
					circlePoints.add(points.get(s1[0]));
					circlePoints.add(points.get(s1[1]));
					circlePoints.add(points.get(s1[2]));
					if ( (c1.x-points.get(aset[1]).x)*(c1.x-points.get(aset[1]).x) + (c1.y-points.get(aset[1]).y)*(c1.y-points.get(aset[1]).y)
							<= r1*r1 ){  //compare these two norm^2
						center.x = c1.x;
						center.y = c1.y;
						radius = r1;
						
						int swap = aset[1];
						aset[1] = aset[0];
						aset[0] = iset[k_index];
						iset[k_index] = swap; 
						continue;
					}
					
					s1[0] = aset[0];
					s1[1] = aset[1];
					s1[2] = iset[k_index];
					enc3(points.get(s1[0]), points.get(s1[1]), points.get(s1[2]), c1);
					r1 = Math.sqrt( (c1.x-points.get(s1[0]).x)*(c1.x-points.get(s1[0]).x) + (c1.y-points.get(s1[0]).y)*(c1.y-points.get(s1[0]).y) );
					circlePoints.clear();
					circlePoints.add(points.get(s1[0]));
					circlePoints.add(points.get(s1[1]));
					circlePoints.add(points.get(s1[2]));
					if ( (c1.x-points.get(aset[2]).x)*(c1.x-points.get(aset[2]).x) + (c1.y-points.get(aset[2]).y)*(c1.y-points.get(aset[2]).y)
							<= r1*r1 ){  //compare these two norm^2
						center.x = c1.x;
						center.y = c1.y;
						radius = r1;
						
						int swap = aset[2];
						aset[2] = aset[1];
						aset[1] = aset[0];
						aset[0] = iset[k_index];
						iset[k_index] = swap;
						continue;
					}
					tol *= 2;
				}
			}		
		}
		//System.out.println("lalala"+center.x+", "+center.y);
		return radius; // "center is already stored"
	}
	
	double computeTol (ArrayList<SSTD_Point> points){
		// compute the tolarance
		int size = points.size();
		double mean_x = 0;
		double mean_y = 0;
		double temp_x, temp_y;
		double max_x = 0;
		double max_y = 0;
		
		for (int i = 0; i < size; i++){
			mean_x += points.get(i).x;
			mean_y += points.get(i).y;
		}
		mean_x /= size;
		mean_y /= size;
		
		for (int i = 0; i < size; i++){
			temp_x = Math.abs(mean_x - points.get(i).x);
			if (temp_x > max_x){
				max_x = temp_x;
			}
			temp_y = Math.abs(mean_y - points.get(i).y);
			if (temp_y > max_y){
				max_y = temp_y;
			}
		}
		return 0.00000000000000001 * (max_x+max_y);
	}
	
	void enc3(SSTD_Point A, SSTD_Point B, SSTD_Point C, SSTD_Point center){
		// given 3 points, return the center of their minimal bounding circle
		double yDelta_a = B.y - A.y;
		double xDelta_a = B.x - A.x;
		double yDelta_b = C.y - B.y;
		double xDelta_b = C.x - B.x;

		double aSlope = yDelta_a / xDelta_a;
		double bSlope = yDelta_b / xDelta_b;

		SSTD_Point AB_Mid = new SSTD_Point((A.x+B.x)/2, (A.y+B.y)/2);
		SSTD_Point BC_Mid = new SSTD_Point((B.x+C.x)/2, (B.y+C.y)/2);

		if(yDelta_a == 0)         //aSlope == 0
		{
		    center.x = AB_Mid.x;
		    if (xDelta_b == 0)         //bSlope == INFINITY
		    {
		        center.y = BC_Mid.y;
		    }
		    else
		    {
		        center.y = BC_Mid.y + (BC_Mid.x-center.x)/bSlope;
		    }
		}
		else if (yDelta_b == 0)               //bSlope == 0
		{
		    center.x = BC_Mid.x;
		    if (xDelta_a == 0)             //aSlope == INFINITY
		    {
		        center.y = AB_Mid.y;
		    }
		    else
		    {
		        center.y = AB_Mid.y + (AB_Mid.x-center.x)/aSlope;
		    }
		}
		else if (xDelta_a == 0)        //aSlope == INFINITY
		{
		    center.y = AB_Mid.y;
		    center.x = bSlope*(BC_Mid.y-center.y) + BC_Mid.x;
		}
		else if (xDelta_b == 0)        //bSlope == INFINITY
		{
		    center.y = BC_Mid.y;
		    center.x = aSlope*(AB_Mid.y-center.y) + AB_Mid.x;
		}
		else
		{
		    center.x = (aSlope*bSlope*(AB_Mid.y-BC_Mid.y) - aSlope*BC_Mid.x + bSlope*AB_Mid.x)/(bSlope-aSlope);
		    center.y = AB_Mid.y - (center.x - AB_Mid.x)/aSlope;
		}

		//if (center.x == 0 && center.y == 0){System.out.println("00000");}
	}
	
	double LikelihoodRatio(int ctot, double c, double B){
		double logc = Math.log(c);
		double logb = Math.log(B);
		double logctotc = Math.log((double)ctot-c);
		double logctotb = Math.log((double)ctot -B);
		double logLR;
		
		if (c > B){
			logLR = c * (logc - logb) + (double)(ctot - c) * (logctotc-logctotb);
		}
		else{
			logLR = 0;
		}
		return logLR;
	}
	
	int pdist2(ArrayList<SSTD_Point> tmpSet, double appX, double appY){
		// "pdist2" in Matlab code, returns index of point in tmpSet farthest from the center
		int size = tmpSet.size();
		int poi = 0;
		double max_square = 0;
		double temp_square;
		
		for (int i = 0; i < size; i++){
			temp_square = (tmpSet.get(i).x - appX)*(tmpSet.get(i).x - appX) + (tmpSet.get(i).y - appY)*(tmpSet.get(i).y - appY);
			if (temp_square > max_square ){
				max_square = temp_square;
				poi = i;
			}
		}
		return poi;
	}
	
	public void circle (String dataset_path) throws NumberFormatException, IOException{
	// main method of the algorithm
		
		//for (int loop = 0; loop < 1000; loop++){
		
		Dataset myDataset = new Dataset();

		int ctot; // amount of points in the dataset
		double studyArea = (xmax - xmin) * (ymax - ymin);
	
		int maxRadius;	
		
		ctot = myDataset.readDataset(dataset_path); // read data from .csv file
		System.out.println("Activity Set |A| = " + ctot);
		System.out.println("Study Area S = " + studyArea);
		
		int xCells = (int) Math.ceil((xmax - xmin) / cellSize);
		//double xmaxNew = xCells * cellSize + xmin;
		int yCells = (int) Math.ceil((ymax - ymin) / cellSize);
		//double ymaxNew = yCells * cellSize + ymin;
		
		//System.out.println("x: "+xCells+" y: "+yCells+" ctot: "+ctot);
		
		int cellArea = cellSize * cellSize;
		double cellBaseline = ((double)(ctot * cellArea)) / studyArea; // expected number of points in a cell
		boolean flag = true;	
		// pointSet is "pointSet_x and pointSet_y here"
		int counter = 0;
		
		ArrayList<ArrayList<SSTD_Point>> prunedSets = new ArrayList<ArrayList<SSTD_Point>>(); // prunedSet in Matlab code
		
		// create and initialize N*N count grid, where N is number of cells on one axis
		CountGridCells[][] countGridCells = new CountGridCells[xCells][yCells];
		for (int i = 0; i < xCells; i++)
			for (int j = 0; j < yCells; j++)
				countGridCells[i][j] = new CountGridCells();
						
		for (int i = 0; i < ctot; i++){
			int xCoord = Math.abs((int)Math.floor((Dataset.pointSet.get(i).x - xmin) / cellSize)); // which cell a point goes to
			int yCoord = Math.abs((int)Math.floor((Dataset.pointSet.get(i).y - ymin) / cellSize));
			//System.out.println("i:" + i+" xCoord = "+xCoord+" yCoord = "+yCoord);

			if (xCoord >= xCells)
				xCoord--;
			if (yCoord >= yCells)
				yCoord--;

			countGridCells[xCoord][yCoord].points.add(Dataset.pointSet.get(i));
			countGridCells[xCoord][yCoord].countGrid++;
		}
		
		// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	    // Circle_Grid Generation. Circles in param. space will be stored in here.
		if (xCells < yCells)
			maxRadius = xCells / 2;
		else
			maxRadius = yCells / 2;
		
		ArrayList<double[]> ulrList = new ArrayList<double[]>();
		
		long startTime1 = System.currentTimeMillis();

		while (flag == true){
			double[][][] circleGridB = new double[xCells][yCells][maxRadius]; // number of filled cells
			double[][][] circleGridc = new double[xCells][yCells][maxRadius]; // number of points in filled cells
			double[][] circleGridULR = new double[xCells * yCells * maxRadius][8];
			int[][] filledCircle = new int[xCells][yCells];
			//System.out.println("xcell: "+xCells+" ycell: "+yCells+ " maxRadius: "+maxRadius);
			int ii = 0;
			for (int i = 0; i < xCells; i++){
				for (int j = 0; j < yCells; j++){
					for (int k = 0; k < maxRadius; k++){
						circleGridB[i][j][k] = 0;	
						circleGridc[i][j][k] = 0;
						for (int l = 0; l < 8; l++){
							circleGridULR[ii][l] = 0;
						}
						ii++;
					}
				}
			}
			
			//int[][] columnsCircleGrid = new int[xCells][yCells]; // indicates column index of a grid
			//int[][] rowsCircleGrid = new int[xCells][yCells];    // indicates row index of a grid
			//meshgrid(columnsCircleGrid, rowsCircleGrid, yCells, xCells); // order of parameters swapped from origin


			ii = 0;
			for (int x = 0; x < xCells; x++){ // enumerate circle center
				for (int y = 0; y < yCells; y++){
					for (int radius = 1; radius <= maxRadius; radius++){
						// the following 2-layer loop is to compute "filledCircle"
						for (int k = 0; k < xCells; k++){
							for (int l = 0; l < yCells; l++){
								/*
								if ((rowsCircleGrid[k][l] - x) * (rowsCircleGrid[k][l] - x) + (columnsCircleGrid[k][l] - y) * (columnsCircleGrid[k][l] - y) < radius * radius){
									filledCircle[k][l] = 1;
									circleGridB[x][y][radius-1]++; // be cautious about this shift, 0 indicates a radius of 1
									circleGridc[x][y][radius-1] += countGridCells[k][l].countGrid;
								}
								else {
									filledCircle[k][l] = 0;
								}
								*/
								if ( ((k - x)*(k - x) + (l - y)*(l - y)) < radius * radius ){
									filledCircle[k][l] = 1;
									circleGridB[x][y][radius-1]++;
									circleGridc[x][y][radius-1] += countGridCells[k][l].countGrid;
								}
							}
						}
						if (radius >= 2){
							double lB = circleGridB[x][y][radius-2] * cellBaseline; // maximum bounded
							double lc = circleGridc[x][y][radius-2];
							double uB = circleGridB[x][y][radius-1] * cellBaseline; // minimal bounding
							double uc = circleGridc[x][y][radius-1];
							
							if (lc < uB){
								lc = uB;
							}
							if ((uB <= ctot) && (lB > 0) && (uc >= lB)){
								double lru = lrupper(ctot, uc, uB, lc, lB);				
								if (lru >= theta){
									circleGridULR[ii][0] = x; // cell index of the center
									circleGridULR[ii][1] = y;
									circleGridULR[ii][2] = radius;
									circleGridULR[ii][3] = lc;
									circleGridULR[ii][4] = lB;
									circleGridULR[ii][5] = uc;
									circleGridULR[ii][6] = uB;
									circleGridULR[ii][7] = lru;
									ii++;
								}
							}
						}
					}
				}
			}	

			//System.out.println("ii = "+ii);
			
			for (int i = 0; i < ii; i++){
				if (circleGridULR[i][7] <= theta){
					//System.out.println("fjdkfjkd\ndfd\nfdfd\ndfd\n");
					for (int j = 0; j < 8; j++){
						circleGridULR[i][j] = 0;
					}
				}
			}
			// equivalent to "sortrows with the 8th column"
			Arrays.sort(circleGridULR, new Comparator<double[]>() {
				public int compare(double[] a, double[] b) {
					return -Double.compare(a[7], b[7]);
				}
			});
			
			/*
			for (int i = 0; i < ii; i++){
				System.out.println("No."+i+" circleGridULR[0] = "+circleGridULR[i][0] + " circleGridULR[2] = "+circleGridULR[i][2]);
			}
			*/
			
			// nRows is equal to unique rows in circleGridULR
			int nRows = unique(circleGridULR, xCells * yCells * maxRadius, 8);
					
			if (nRows > 0){
				prunedSets.add(new ArrayList<SSTD_Point>());
				writePointList (circleGridULR[0], countGridCells, filledCircle, prunedSets, xCells, yCells, counter);
				ulrList.add(circleGridULR[0]);
				// new line??
				setdiff(prunedSets, counter); // delete points that are in prunedSet[counter] from Pointset
			}
			if ( (Dataset.pointSet.size() <= 0) || ((circleGridULR[0][0] == 0)&&(circleGridULR[0][1] == 0)&&(circleGridULR[0][2] == 0)&&(circleGridULR[0][3] == 0)&&(circleGridULR[0][4] == 0)&&(circleGridULR[0][5] == 0)&&(circleGridULR[0][6] == 0)&&(circleGridULR[0][7] == 0)) ) {
				flag = false;
			}
			counter++; // each loop find one pruned set, which is one potential circle
		}

		// this loop is used as "~cellfun(@isempty, prunedSets)" in Matlab code
		for (int i = prunedSets.size()-1; i >= 0; i--){
			if (prunedSets.get(i).size() == 0){
				prunedSets.remove(i);
			}
		}
		
		/*
		// for debugging
		for (int i = 0; i < prunedSets.size(); i++){
			System.out.println("pruneSet + " +i+ " size: " + prunedSets.get(i).size());
			for (int j = 0; j < prunedSets.get(i).size(); j++){
				System.out.print("("+prunedSets.get(i).get(j).x+", "+prunedSets.get(i).get(j).y+") ");
			}
			System.out.println();
		}
		*/
		
		// ------------------------------- refine phase starts -----------------------
		System.out.println("Refine phase starts!");
		System.out.println("pruned set size: "+ prunedSets.size());
		
		long startTime2 = System.currentTimeMillis();
		
		if (prunedSets.size() > 0){ // one element in prunedSet indicates one potential circle
			double[][] finalLR = new double[prunedSets.size()][7];
			double radius = Double.MAX_VALUE; // this "radius" different from the previous one?
			
			for (int i = 0; i < prunedSets.size(); i++){ // check each candidate circle
				
				ArrayList<SSTD_Point> tmpSet = new ArrayList<SSTD_Point>(); // clone prunedSets{i}
				ArrayList<SSTD_Point> circlePoints = new ArrayList<SSTD_Point>();
				
				tmpSet.clear();
				for (int j = 0; j < prunedSets.get(i).size(); j++){
					tmpSet.add(prunedSets.get(i).get(j)); // clone prunedSets element
				}
				
				double previousLR = 0; // used to get the max LR for a pruned set
				boolean newflag = true;
				while (newflag == true){
					SSTD_Point center = new SSTD_Point();
					radius = minCircle(tmpSet, true, center, circlePoints);
					//System.out.println("jijiji " + center.x+", "+center.y);
					//System.out.println("radius = " + radius);

					double c = (double)tmpSet.size();
					double area = radius * radius * Math.PI;
					double B = area * ctot / studyArea;
					double likelihoodRatio = LikelihoodRatio(ctot, c, B);
					if ( (likelihoodRatio > previousLR) && (c >= 2) ){ // new max LR circle

						finalLR[i][0] = center.x;
						finalLR[i][1] = center.y;
						finalLR[i][2] = radius;
						finalLR[i][3] = area;
						finalLR[i][4] = B;
						finalLR[i][5] = c;
						finalLR[i][6] = likelihoodRatio;
						previousLR = likelihoodRatio;
					}
					if ( (tmpSet.size() > 0) && (circlePoints.size() > 0) ){
						double appX = xmin + ulrList.get(i)[0] * cellSize - ((double)cellSize) / 2; // center of circle candidate
						double appY = ymin + ulrList.get(i)[1] * cellSize - ((double)cellSize) / 2;
						int poi = pdist2(tmpSet, appX, appY);
						tmpSet.remove(poi);
					}
					else {
						newflag = false;
					}
				}			
			}	
			for (int i = 0; i < prunedSets.size(); i++){
				if (finalLR[i][6] > theta){
					System.out.println("No." + i + ":center: ("+finalLR[i][0]+", "+finalLR[i][1]+") radius: "+finalLR[i][2]
							+" area: "+finalLR[i][3]+ " nPoints: " + finalLR[i][5] +" LR: "+finalLR[i][6]);
				}
			}
		}
		else{
			System.out.println("nothing found!");
		}
		long startTime3 = System.currentTimeMillis();
		
		/*
		for (int i = loop / 200; i < 5; i++){
			SigCirMain.time_prune[i] += (startTime2-startTime1);
			SigCirMain.time_refine[i] += (startTime3 - startTime2);
		}
		*/
		
		System.out.println("prune: "+(startTime2-startTime1)+" milliseconds");
		System.out.println("refine: "+(startTime3-startTime2)+" milliseconds");
		System.out.println("GPR total: "+(startTime3-startTime1)+" milliseconds");
		
	}	
	
	public void test(SSTD_Point p){
		// no meaning in this program, just for test
		p.x = 10;
		p.y = 20;
	}
	
	public static void main(String[] args) throws NumberFormatException, IOException {
		
		/*
		for (int i = 0; i < 5; i++){
			SigCirMain.time_prune[i] = 0;
			SigCirMain.time_refine[i] = 0;
		}
		*/
		
		SigCirMain mySigCirMain = new SigCirMain();
		//System.out.println("200K: ");
		//mySigCirMain.circle(mySigCirMain.dataset_path_200K);
		//System.out.println("300K: ");
		//mySigCirMain.circle(mySigCirMain.dataset_path_300K);
		//System.out.println("400K: ");
		//mySigCirMain.circle(mySigCirMain.dataset_path_400K);
		//System.out.println("500K: ");
		//mySigCirMain.circle(mySigCirMain.dataset_path_500K);
		//System.out.println("600K: ");
		
		//mySigCirMain.circle(mySigCirMain.dataset_toy);
		
		/*
		SSTD_Point A = new SSTD_Point(0,2);
		SSTD_Point B = new SSTD_Point(0,0);
		SSTD_Point C = new SSTD_Point(0,1);
		SSTD_Point center = new SSTD_Point();
		mySigCirMain.enc3(A, B, C, center);
		System.out.println(center.x+" "+center.y);
		*/
		
		
		ArrayList<SSTD_Point> points = new ArrayList<SSTD_Point>();
		SSTD_Point center = new SSTD_Point();
		ArrayList<SSTD_Point> tmpSet = new ArrayList<SSTD_Point>();
		ArrayList<SSTD_Point> circlePoints = new ArrayList<SSTD_Point>();	

		points.add(new SSTD_Point(0,2));
		points.add(new SSTD_Point(2,2));
		points.add(new SSTD_Point(0,0));
		points.add(new SSTD_Point(2,0));

		points.add(new SSTD_Point(1,1));
		
		double radius = mySigCirMain.minCircle_SEC(points, tmpSet, center, circlePoints);
		System.out.println("SEC result: " + center.x + ", " + center.y + ", " + radius);
		
		for (int i = 0; i < circlePoints.size(); i++){
			System.out.println(circlePoints.get(i).x + ", "+circlePoints.get(i).y);
		}
		
		
		/*
		for (int i = 0; i < 5; i++){
			System.out.println((i+1)*200+" MC simulations: prune: "+time_prune[i]+" refine: "+time_refine[i]);
		}
		*/
		/*
		SSTD_Point cen = new SSTD_Point();
		cen = mySigCirMain.enc3(new SSTD_Point(0,0), new SSTD_Point(2,0), new SSTD_Point(1,3));
		System.out.println("x: "+cen.x+"y: "+cen.y);
		double[] aset = new double[3];
		aset[0] = 2;aset[1] = 1;aset[2] = 4;
		Arrays.sort(aset);
		System.out.println(aset[0]+" "+aset[1]+""+aset[2]);
		*/
	}
}
