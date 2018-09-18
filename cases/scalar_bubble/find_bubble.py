# include microHH tools
execfile('../../python/microhh_tools.py')
filetoread='myscalar.0000020'

# read binary file
res=read_restart_file(filetoread,64,64,64,endian='little')

# find bubble
import numpy as np
bubble=np.where(res!=0)
print('Current bubble position (t=20):',bubble,'\n')
print('Current bubble magnitude (t=20):',res[bubble],'\n')
