import numpy as np
from solar_mod import *

phi= 45.5
ym = 26.0
zm = 45.0
x1 = -8.0
x2 = 52.0
tan1 = ym/zm
gams1 = arctand(x1/zm)
gams3 = arctand(x2/zm)
als1 = arctand(tan1*cosd(gams1))
als3 = arctand(tan1*cosd(gams3))
f_out= open('fichier_ombrage.txt','w');
line = '   001     ! for each orientation (or opening)'  + '\n'
f_out.write(line)
line = '   045     ! Slope for each orientation' + '\n'
f_out.write(line)
line = '   000     ! Azimuth for each orientation'  + '\n'
f_out.write(line)
step = 5
azi = np.arange(-180,180,step)
n = len(azi)  # Ce nombre doit être <= 80 pour TRNSYS
for i in range(0,n):
    strr = '%.1f' % azi[i]
    f_out.write(strr)
    f_out.write(' ')
f_out.write('\n')
strr = '! Obstruction height for orientation 1, angle  '
for i in range(0,n):
    str2 = strr + '%d'% i + ' (' + '%.2f' % azi[i] +  'to ' + '%.2f' % (azi[i] + step) + ')'
    ga = azi[i] + step/2.0
    if ga >= gams1 and ga <= gams3 :
        al = arctand(ym/zm*cosd(ga))
        str3 =  '%.2f' % al   + str2 + '\n'
    else:
        str3 = '%.2f' % 0   + str2 + '\n'
    f_out.write(str3)
f_out.close()

