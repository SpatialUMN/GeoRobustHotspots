package georobust;
import java.io.IOException;

public class RunGeo {
	static String path = "C:/Users/JingWe/Desktop/capstone/activity.csv";
	static int cellSize = 15;
	static double theta = 100;
	
	public static void main(String[] args) throws NumberFormatException, IOException {
		
		
		SigCirMain mySigCirMain = new SigCirMain();

		
		mySigCirMain.circle(path);
		
		
	}

}
