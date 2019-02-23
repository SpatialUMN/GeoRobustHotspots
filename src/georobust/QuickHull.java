package georobust;

//This is a java program to find a points in convex hull using quick hull method
import java.util.ArrayList;
import java.util.Scanner;

//this methods will change the input arraylist (i.e., points)

public class QuickHull
{
  public ArrayList<SSTD_Point> quickHull(ArrayList<SSTD_Point> points)
  {
      ArrayList<SSTD_Point> convexHull = new ArrayList<SSTD_Point>();
      
      //if (points.size() < 3)
          //return (ArrayList) points.clone();

      int minPoint = -1, maxPoint = -1;
      double minX = Double.MAX_VALUE;
      double maxX = Double.MIN_VALUE;
      for (int i = 0; i < points.size(); i++)
      {
          if (points.get(i).x < minX)
          {
              minX = points.get(i).x;
              minPoint = i;
          }
          if (points.get(i).x > maxX)
          {
              maxX = points.get(i).x;
              maxPoint = i;
          }
      }
      SSTD_Point A = new SSTD_Point(points.get(minPoint).x, points.get(minPoint).y);
      SSTD_Point B = new SSTD_Point(points.get(maxPoint).x, points.get(maxPoint).y);

      convexHull.add(A);
      convexHull.add(B);
      points.remove(A);
      points.remove(B);

      ArrayList<SSTD_Point> leftSet = new ArrayList<SSTD_Point>();
      ArrayList<SSTD_Point> rightSet = new ArrayList<SSTD_Point>();

      for (int i = 0; i < points.size(); i++)
      {
          SSTD_Point p = points.get(i);
          if (pointLocation(A, B, p) == -1)
              leftSet.add(p);
          else if (pointLocation(A, B, p) == 1)
              rightSet.add(p);
      }
      hullSet(A, B, rightSet, convexHull);
      hullSet(B, A, leftSet, convexHull);

      return convexHull;
  }

  public double distance(SSTD_Point A, SSTD_Point B, SSTD_Point C)
  {
      double ABx = B.x - A.x;
      double ABy = B.y - A.y;
      double num = ABx * (A.y - C.y) - ABy * (A.x - C.x);
      if (num < 0)
          num = -num;
      return num;
  }

  public void hullSet(SSTD_Point A, SSTD_Point B, ArrayList<SSTD_Point> set,
          ArrayList<SSTD_Point> hull)
  {
      int insertPosition = hull.indexOf(B);
      if (set.size() == 0)
          return;
      if (set.size() == 1)
      {
          SSTD_Point p = set.get(0);
          set.remove(p);
          hull.add(insertPosition, p);
          return;
      }
      double dist = Double.MIN_VALUE;
      int furthestPoint = -1;
      for (int i = 0; i < set.size(); i++)
      {
          SSTD_Point p = set.get(i);
          double distance = distance(A, B, p);
          if (distance > dist)
          {
              dist = distance;
              furthestPoint = i;
          }
      }
      SSTD_Point P = set.get(furthestPoint);
      set.remove(furthestPoint);
      hull.add(insertPosition, P);

      // Determine who's to the left of AP
      ArrayList<SSTD_Point> leftSetAP = new ArrayList<SSTD_Point>();
      for (int i = 0; i < set.size(); i++)
      {
          SSTD_Point M = set.get(i);
          if (pointLocation(A, P, M) == 1)
          {
              leftSetAP.add(M);
          }
      }

      // Determine who's to the left of PB
      ArrayList<SSTD_Point> leftSetPB = new ArrayList<SSTD_Point>();
      for (int i = 0; i < set.size(); i++)
      {
          SSTD_Point M = set.get(i);
          if (pointLocation(P, B, M) == 1)
          {
              leftSetPB.add(M);
          }
      }
      hullSet(A, P, leftSetAP, hull);
      hullSet(P, B, leftSetPB, hull);

  }

  public int pointLocation(SSTD_Point A, SSTD_Point B, SSTD_Point P)
  {
      double cp1 = (B.x - A.x) * (P.y - A.y) - (B.y - A.y) * (P.x - A.x);
      if (cp1 > 0)
          return 1;
      else if (cp1 == 0)
          return 0;
      else
          return -1;
  }

  public static void main(String args[])
  {
  	/*
      System.out.println("Quick Hull Test");
      Scanner sc = new Scanner(System.in);
      System.out.println("Enter the number of points");
      int N = sc.nextInt();

      ArrayList<SSTD_Point> points = new ArrayList<SSTD_Point>();
      System.out.println("Enter the coordinates of each points: <x> <y>");
      for (int i = 0; i < N; i++)
      {
          int x = sc.nextInt();
          int y = sc.nextInt();
          SSTD_Point e = new SSTD_Point(x, y);
          points.add(i, e);
      }

      QuickHull qh = new QuickHull();
      ArrayList<SSTD_Point> p = qh.quickHull(points);
      System.out
              .println("The points in the Convex hull using Quick Hull are: ");
      for (int i = 0; i < p.size(); i++)
          System.out.println("(" + p.get(i).x + ", " + p.get(i).y + ")");
      sc.close();
      */
  	
  	ArrayList<SSTD_Point> points = new ArrayList<SSTD_Point>();
  	points.add(new SSTD_Point(0, 0));
  	points.add(new SSTD_Point(1, 0.1));
  	points.add(new SSTD_Point(2, 2));
  	points.add(new SSTD_Point(1, 1));
  	points.add(new SSTD_Point(0, 1));

  	QuickHull q = new QuickHull();
  	points = q.quickHull(points);
  }
}