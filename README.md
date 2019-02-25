# Geographically Robust Hotspots Detection  
* [What Can You Get From Geographically Robust Hotspots Detection](https://github.com/SpatialUMN/GeoRobustHotspots/blob/master/README.md#what-can-you-get)  
* [Usage](https://github.com/SpatialUMN/GeoRobustHotspots/blob/master/README.md#usage)   
  * [Input Data Format](https://github.com/SpatialUMN/GeoRobustHotspots/blob/master/README.md#input-data-format)  
  * [Download and Run](https://github.com/SpatialUMN/GeoRobustHotspots/blob/master/README.md#Download-and-Run)  
    * [How to import a GitHub project into Eclipse](https://github.com/collab-uniba/socialcde4eclipse/wiki/How-to-import-a-GitHub-project-into-Eclipse)  
    * [Set Variables](https://github.com/SpatialUMN/GeoRobustHotspots/blob/master/README.md#set-variables) 
* [Code Explanation (Java Diagram)]() 
* [Case Study]()  
* [Bug Report](https://github.com/SpatialUMN/GeoRobustHotspots/issues)  
  

# What Can You Get
Elliptical Hotspot Detection (EHD) finds ellipse shaped hotspot areas where the concentration of activities inside is significantly higher
than the concentration of activities outside.   
![E1b](https://github.com/SpatialUMN/GeoRobustHotspots/blob/master/images/E1b.PNG)  
See [application domain](https://github.com/SpatialUMN/GeoRobustHotspots/wiki/Application-Domain) to see where you can use elliptical hotspot detection

# Usage  
## Input Data Format  
We need 1 input file `activity`. It has 3 attributes:  
`ID` is the activity id.   
`X` is the x-axis value of the activity.  
`Y ` is the Y-axis value of the activity.  

## Download and Run  
### [How to import a GitHub project into Eclipse](https://github.com/collab-uniba/socialcde4eclipse/wiki/How-to-import-a-GitHub-project-into-Eclipse)  
### Set Variables   
Open [`RunElliptic.java`](https://github.com/SpatialUMN/EllipticalHotspots/blob/master/src/elliptical/RunElliptic.java) file, change line 6, 7, and 8.  
`dataset_path` is the path to your activity file.  
`Method` set Method = 1 if you want to use grid method, set Method = 0 if you need naive method.  
`step_size` only has effect if you choose naive method. It is the step length used on denominator.   
