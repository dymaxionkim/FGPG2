# FGPG2

_Fine Involute Gear Profile Generator 2_

_with python3_

![./Result/Result.png]

![./Result/Result2.png]


## Pre-requisites

```
pip install numpy
pip install matplotlib
pip install ezdxf
pip install pysimplegui
```

## Build exe in MS Windows 10

* Development version of pyinstaller is needed because of matplotlib matching

```
pip install https://github.com/pyinstaller/pyinstaller/archive/develop.zip
pyinstaller -w -F FGPG2.py
```

## Using on Android by Pydroid3

* Install Pydroid3
* Install some libraries by Pip menu in Pydroid (numpy, matplotlib, ezdxf, pysimplegui)
* Install mgit
* clone this repository
* Load FGPG2.py and run

![./img/Screenshot_Pydroid3.png]

![./img/Screenshot_Android.png]


