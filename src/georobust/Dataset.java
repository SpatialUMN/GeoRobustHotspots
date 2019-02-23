package georobust;
import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

public class Dataset {
	
	static ArrayList<SSTD_Point> activitySet = new ArrayList<SSTD_Point>(); // activitySet
	static ArrayList<SSTD_Point> pointSet = new ArrayList<SSTD_Point>(); // pointSet
	
	int readDataset(String path) throws NumberFormatException, IOException{
		// read CSV dataset into two arraylists, returns number of activities
		BufferedReader br = null;
		String splitBy = ","; // numbers split by ","
		String line ="";
		
		br = new BufferedReader(new FileReader(path));
		br.readLine(); // skip first line of title
			
		activitySet.clear();
		pointSet.clear();
		while ((line = br.readLine()) != null) { // read dataset line by line
			String [] rowContent = line.split(splitBy);
			activitySet.add(new SSTD_Point(Double.parseDouble(rowContent[1]), Double.parseDouble(rowContent[2]))); // second and third element stored into point_x
		}
		br.close();
		
		for (int i = 0; i < Dataset.activitySet.size(); i++){ // duplicate the activity set
			Dataset.pointSet.add(Dataset.activitySet.get(i));
		}
		
		return Dataset.activitySet.size(); // return size of dataset
	}
	
	public static void main(String[] args) throws NumberFormatException, IOException {
		Dataset test = new Dataset();
		int size = test.readDataset("/Users/tangxun/Desktop/SSTD_2015/circletoy150.csv");
		System.out.println("size: "+size);
		for (int i=0; i < 130; i++){
			System.out.println(i+": "+Dataset.activitySet.get(i).x+" , "+Dataset.pointSet.get(i).x);
			//System.out.println(i+": "+Dataset.activitySet_x.get(i)+", "+", "+Dataset.activitySet_y.get(i));
		}
	
	}
}
