# coherentImaging
Simple simulations to understand effects of coherent imaging.

## Overview
A point object emits spherical waves. A lens intercepts a portion of this spherical wave, and refocuses it onto a blurred point in the image plane. The point spread function (PSF) describes the response of an imaging system to a point source or point object [1]. For a single lens, an on-axis point source in the object plane produces an Airy disc PSF in the image plane. Below are some simulations of coherent absorption imaging of a ring shaped atomic cloud. The probe used is gaussian in shape.

<table>
<tr>
<td> <div align="left"> NA = 0.1 </div></td>
<td><img src="/images/ImageBlur_ring_Probe_gaussInt_NA10E-2.png" alt="Drawing" width="400"/> </td>
</tr>
<tr>
<td> <div align="left"> NA = 0.25 </div></td>
<td><img src="/images/ImageBlur_ring_Probe_gaussInt_NA25E-2.png" alt="Drawing" width="400"/> </td>
</tr>
<tr>
<td> <div align="left"> NA = 0.40 </div></td>
<td><img src="/images/ImageBlur_ring_Probe_gaussInt_NA40E-2.png" alt="Drawing" width="400"/> </td>
</tr>
</table>

## Contents
This repo contains the source code, some use case examples and sample data to play with.
- *../src/* contains the source code with all modules.
- *../test/* contains the examples.

## Getting Started

- Clone this repository. 
```git clone https://github.com/sbanik1/coherentImaging <lcl_dir>```
- Modify the startup code *../test/startUp.py* as follows. 
  - Change the  path directory  for the command ```sys.path.insert(1, '../coherentImaging/src')``` to the path for the directory *../src/* on your local machine.
  - Change the variable *saveDir* and set it to the path to directory *../test/out/* on your local machine.
<table>
<tr>
<td><img src="/images/Installation.png" alt="Drawing" width="700"/> </td>
</tr>
</table>

- Run the examples in the directory *../test/*, by either running them on an IDE or via the command line.
  - Change directory ```cd "../coherentImaging/test"```
  - Run the example ```ImageBlurring.py``
  
  
 ## References
<a id="1">[1]</a> 
“Point spread function.” Wikipedia, Wikimedia Foundation, 7 Dec. 2020, https://en.wikipedia.org/wiki/Point_spread_function.






