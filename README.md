# Geographically Robust Hotspots Detection  
* [What Can You Get](https://github.com/SpatialUMN/GeoRobustHotspots/blob/master/README.md#what-can-you-get)  
* [Usage](https://github.com/SpatialUMN/GeoRobustHotspots/blob/master/README.md#usage)   
  * [Data Format](https://github.com/SpatialUMN/GeoRobustHotspots/blob/master/README.md#data-format)  
  * [Download and Run](https://github.com/SpatialUMN/GeoRobustHotspots/blob/master/README.md#Download-and-Run)  
    * [How to import a GitHub project into Eclipse](https://github.com/collab-uniba/socialcde4eclipse/wiki/How-to-import-a-GitHub-project-into-Eclipse)  
    * [Set Variables](https://github.com/SpatialUMN/GeoRobustHotspots/blob/master/README.md#set-variables) 
    * [Output Result](https://github.com/SpatialUMN/GeoRobustHotspots/blob/master/README.md#output-result)
* [Code Explanation - Java Diagram](https://github.com/SpatialUMN/GeoRobustHotspots/wiki/Java-Class-Diagram) 
* [Case Study](https://github.com/SpatialUMN/GeoRobustHotspots/wiki/Case-Study)  
* [Bug Report](https://github.com/SpatialUMN/GeoRobustHotspots/issues)  
* [Link to Paper](https://ieeexplore.ieee.org/abstract/document/7395840)
  

# What Can You Get
Geographically Robust Hotspot Detection (GRHD) finds hotspot areas where the concentration of points inside is significantly high.  
![G1](https://github.com/SpatialUMN/GeoRobustHotspots/blob/master/image/G1.PNG)  
See [application domain](https://github.com/SpatialUMN/GeoRobustHotspots/wiki/Application-Domain) to see where you can use GRHD.   
Basic [concepts](https://github.com/SpatialUMN/GeoRobustHotspots/wiki/Basic-Concepts) can help better understanding the problem.  


# Usage  
## Data Format  
We need 1 input file `activity`. It has 3 attributes:  
`ID` is the activity id.   
`X` is the x-axis value of the activity.  
`Y ` is the Y-axis value of the activity.  

## Download and Run  
### [How to import a GitHub project into Eclipse](https://github.com/collab-uniba/socialcde4eclipse/wiki/How-to-import-a-GitHub-project-into-Eclipse)  
### Set Variables   
Open [`RunGeo.java`](https://github.com/SpatialUMN/GeoRobustHotspots/blob/master/src/georobust/RunGeo.java) file, change line 5, 6, and 7.    
`path` is the path to your activity file.  
`cellSize` divide the input data into multiple cellsize x cellsize squares.  
`theta` is the log likelihood ratio threshold.  

### Output Results   
The output will contain the dataset information, how does the algorithm participate the data, the hotspot result information and the running time.  
[Here](https://github.com/SpatialUMN/GeoRobustHotspots/blob/master/Output) is an outcome example you might see.
