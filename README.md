# FGPG2

_Fine Involute Gear Profile Generator 2 - with python3_and uv package manager


![](./img/Screenshot_GUI.png)


## UI Buttons

* __Working Directory__ : Output files are saved here.
* __Load__ : Load Inputs.dat and every parameters are set.
* __Run__ : Calculate, Plot, Save plotted images, Save dxf, Save csv
* __Toggle__ : Change output image
* __Exit__ : Finish



## Input parameters

* __m__ : Module
* __z__ : Teeth number
* __alpha__ : Pressure angle
* __x__ : Offset factor
* __b__ : Backlash factor
* __a__ : Addendum factor
* __d__ : Dedendum factor
* __c__ : Radius factor of edge round of hob (root)
* __e__ : Radius factor of edge round of tooth (edge)
* __x0__ : Center position
* __y0__ : Center position
* __seg circle__ : Number of contol points for pitch, offset and base circles
* __seg involute__ : Number of contol points for involute curve
* __seg edge r__ : Number of contol points for edge trocoid rounding
* __seg root r__ : Number of contol points for root trocoid rounding
* __seg outer__ : Number of contol points for outer arc
* __seg root__ : Number of contol points for root arc
* __scale__ : Scale factor for one tooth plot



## Output files

* __Inputs.dat__ : Parameters data for UI
* __Result.csv__ : Gear spec data

![](./img/Screenshot_csv.jpg)

* __Result.dxf__ : dxf CAD file for one tooth

![](./img/Screenshot_dxf.jpg)

* __Result.png__ : Gear geometry for whole teeth

![](./img/Result.png)

* __Result2.png__ : Gear geometry for one tooth

![](./img/Result2.png)

## Ref.

* [Calculate Gear Spec](https://tro.kr/11) : Excel Sheet
* [Calculate Gear Spec](https://tro.kr/49) : Web



## Thank you!
