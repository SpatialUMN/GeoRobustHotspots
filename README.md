# Geographically Robust Hotspots Detection  
* [What Can You Get](https://github.com/SpatialUMN/GeoRobustHotspots/blob/master/README.md#what-can-you-get)  
* [Usage](https://github.com/SpatialUMN/GeoRobustHotspots/blob/master/README.md#usage)   
  * [Data Format](https://github.com/SpatialUMN/GeoRobustHotspots/blob/master/README.md#data-format)  
  * [Download and Run](https://github.com/SpatialUMN/GeoRobustHotspots/blob/master/README.md#Download-and-Run)  
    * [How to import a GitHub project into Eclipse](https://github.com/collab-uniba/socialcde4eclipse/wiki/How-to-import-a-GitHub-project-into-Eclipse)  
    * [Set Variables](https://github.com/SpatialUMN/GeoRobustHotspots/blob/master/README.md#set-variables) 
* [Code Explanation - Java Diagram](https://github.com/SpatialUMN/GeoRobustHotspots/wiki/Java-Class-Diagram) 
* [Case Study]()  
* [Bug Report](https://github.com/SpatialUMN/GeoRobustHotspots/issues)  
* [Link to Paper]()
  

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
Open [`RunGeo.java`](https://github.com/SpatialUMN/EllipticalHotspots/blob/master/src/elliptical/RunElliptic.java) file, change line 6, 7, and 8.  
`dataset_path` is the path to your activity file.  
`Method` set Method = 1 if you want to use grid method, set Method = 0 if you need naive method.  
`step_size` only has effect if you choose naive method. It is the step length used on denominator.   
