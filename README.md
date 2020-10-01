# comparatorPloting

## Requirements

PYTHON:
  * numpy: `pip3 install numpy`
  * matplotlib: `pip3 install matplotlib`
  * matplotlib_venn: `pip3 install matplotlib-venn`

## Usage
The command `python3 comparatorPloting.py -h` return:

```
usage: comparatorPloting.py [-h] [-f FPR] directory N [N ...]

positional arguments:
  directory          directory with data
  N                  list of models names used in analisys

optional arguments:
  -h, --help         show this help message and exit
  -f FPR, --FPR FPR  FPR, def=1.9*10^(-4)
```

Example run:

```
python3 comparatorPloting.py \
/dir/with/results \
pwm1 pwm2 pwm3
```

All plots will be writen in `/dir/with/results/plots`
