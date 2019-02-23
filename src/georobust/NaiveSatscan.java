package georobust;
import java.awt.Point;
import java.awt.font.NumericShaper.Range;
import java.io.IOException;
import java.io.ObjectInputStream.GetField;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.Set;


public class NaiveSatscan {
	
	String dataset_path = "/Users/tangxun/Desktop/ellipse_hotspot_detection/datasets/circletoy150.csv"; // dataset used for experiment
	double theta = 10000;

	void knnsearch(int pointIndex, double[] distance){
		// same as "knnsearch" in Matlab code
		int size = Dataset.pointSet.size();
		for (int i = 0; i < size; i++){
			distance[i] = Math.sqrt( (Dataset.pointSet.get(i).x-Dataset.pointSet.get(pointIndex).x)*(Dataset.pointSet.get(i).x-Dataset.pointSet.get(pointIndex).x) + (Dataset.pointSet.get(i).y-Dataset.pointSet.get(pointIndex).y)*(Dataset.pointSet.get(i).y-Dataset.pointSet.get(pointIndex).y) ); // distance stores distances from each point to the objective point		
		}
		Arrays.sort(distance);
	}
	
	void rangesearch(SSTD_Point point, double radius, ArrayList<Integer> idx){
		// same as "rangesearch" in Matlab code, but order of idx is not associated with distances
		int size = Dataset.pointSet.size();
		for (int i = 0; i < size; i++){
			if ( Math.sqrt( (point.x - Dataset.pointSet.get(i).x)*(point.x - Dataset.pointSet.get(i).x) + (point.y - Dataset.pointSet.get(i).y)*(point.y - Dataset.pointSet.get(i).y)) <= radius ){
				idx.add(i);
			}
		}
	}
	
	public static void main(String[] args) throws NumberFormatException, IOException {
		
		
		//long time200 = 0,time400 = 0,time600=0,time800=0,time1000=0;
		//for (int loop = 0; loop < 1000; loop++){
			//System.out.println("No. "+loop+" Monte Carlo Simulation ");
				
		NaiveSatscan myNaiveSatscan = new NaiveSatscan();
		Dataset myDataset = new Dataset();
		SigCirMain mySigCirMain = new SigCirMain();
		
		ArrayList<double[]> finalLR = new ArrayList<double[]>();
		
		int ctot = myDataset.readDataset(myNaiveSatscan.dataset_path);
		double xmax = Dataset.activitySet.get(0).x;
		double xmin = xmax;
		double ymax = Dataset.activitySet.get(0).y;
		double ymin = ymax;
		
		for (int i = 0; i < ctot; i++){
			if (Dataset.activitySet.get(i).x > xmax){
				xmax = Dataset.activitySet.get(i).x;
			}
			if (Dataset.activitySet.get(i).x < xmin){
				xmin = Dataset.activitySet.get(i).x;
			}
			if (Dataset.activitySet.get(i).y > ymax){
				ymax = Dataset.activitySet.get(i).y;
			}
			if (Dataset.activitySet.get(i).y < ymin){
				ymin = Dataset.activitySet.get(i).y;
			}
		}
		double studyArea = (xmax - xmin) * (ymax - ymin);
				
		System.out.println("Activity Set |A| = " + ctot);
		System.out.println("Study Area S = " + studyArea);
		
		double[] RSmax = new double[2];
		
		if ( (xmax - xmin) <= (ymax - ymin) ) {
			RSmax[0] = (xmax - xmin) / 2;
			RSmax[1] = (ymax - ymin) / 2;
			xmax = xmin + ymax - ymin;
			ymax = ymin + ymax - ymin; //????
		}
		else {
			RSmax[0] = (ymax - ymin) / 2;
			RSmax[1] = (xmax - xmin) / 2;
			ymax = ymin + xmax - xmin;
			xmax = xmin + xmax - xmin; //???
		}		
		
		// !!! in this Java code, "pointSet" is equal to "new ActivitySet" in the Matlab code
		int counter = 0;
		boolean flag = true;
		int newCounter = Dataset.activitySet.size(); // varying
		
		long startTime1 = System.currentTimeMillis();

		while (flag == true){
			double[][] tmpLRList = new double[Dataset.pointSet.size()][7];
			for (int i = 0; i < Dataset.pointSet.size(); i++){
				double likelihoodRatioPrevious = 0;
				double likelihoodRatio = 0;
				
				double[] distance = new double[Dataset.pointSet.size()];
				myNaiveSatscan.knnsearch(i, distance);
				
				for (int j = 0; j < Dataset.pointSet.size(); j++){
					int c = j+1;
					double radius = distance[j];
					double area = Math.PI * radius * radius;
					double B = area * ctot / studyArea;
					likelihoodRatio = mySigCirMain.LikelihoodRatio(ctot, c, B);
					if ( (likelihoodRatio > likelihoodRatioPrevious) && (c > 1) ){
						likelihoodRatioPrevious = likelihoodRatio;
						tmpLRList[i][0] = Dataset.pointSet.get(i).x;
						tmpLRList[i][1] = Dataset.pointSet.get(i).y;
						tmpLRList[i][2] = radius;
						tmpLRList[i][3] = area;
						tmpLRList[i][4] = B;
						tmpLRList[i][5] = (double)c;
						tmpLRList[i][6] = likelihoodRatio;
					}
				}
			}
			//System.out.println("newCounter: "+newCounter);

			for(int i = 0; i < Dataset.pointSet.size(); i++) {
				if (tmpLRList[i][6] < myNaiveSatscan.theta){
					tmpLRList[i][0] = 0;
					tmpLRList[i][1] = 0;
					tmpLRList[i][2] = 0;
					tmpLRList[i][3] = 0;
					tmpLRList[i][4] = 0;
					tmpLRList[i][5] = 0;
					tmpLRList[i][6] = 0;
				}
			}
			
			Arrays.sort(tmpLRList, new Comparator<double[]>() {
				public int compare(double[] a, double[] b) {
					return -Double.compare(a[6], b[6]);
				}
			});
			newCounter = mySigCirMain.unique(tmpLRList, Dataset.pointSet.size(), 7);
			
			if (newCounter > 1){
				ArrayList<Integer> idx = new ArrayList<Integer>();
				myNaiveSatscan.rangesearch(new SSTD_Point(tmpLRList[0][0], tmpLRList[0][1]), tmpLRList[0][2], idx);
				for (int i = idx.size()-1; i >= 0; i--){
					Dataset.pointSet.remove((int)idx.get(i)); // this loop used as "setdiff"
				}
				finalLR.add(tmpLRList[counter]);
			}
			//System.out.println("size: "+Dataset.pointSet.size());
			if ( (Dataset.pointSet.size() <= 0) || (newCounter <= 1) ) {
				flag = false;
			}
			//System.out.println("counter: "+counter);
			counter++;
		}
		
		long startTime2 = System.currentTimeMillis();

		System.out.println("SatScan: "+(startTime2-startTime1)+" milliseconds");

		
		/*
		if (loop == 199) {time200 = System.currentTimeMillis() - startTime1;}
		else if (loop == 399) {time400 = System.currentTimeMillis() - startTime1;}
		else if (loop == 599) {time600 = System.currentTimeMillis() - startTime1;}
		else if (loop == 799) {time800 = System.currentTimeMillis() - startTime1;}
		else if (loop == 999) {time1000 = System.currentTimeMillis() - startTime1;}
		}
		*/
		
		/*
		long startTime1000 = System.currentTimeMillis();
		System.out.println("200: "+time200+" milliseconds");
		System.out.println("400: "+time400+" milliseconds");	
		System.out.println("600: "+time600+" milliseconds");	
		System.out.println("800: "+time800+" milliseconds");
		System.out.println("1000: "+time1000+" milliseconds");*/

	}
}
