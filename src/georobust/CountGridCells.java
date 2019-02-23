package georobust;
import java.util.ArrayList;

// this class is equivalent to "countGridCells" plus "countGrid" in Matlab code
public class CountGridCells {

	public ArrayList<SSTD_Point> points; // a cell contains a set of points
	public int countGrid; // number of points
	
	public CountGridCells(){
		this.points = new ArrayList<SSTD_Point>();
		this.countGrid = 0;
	}
}
