# include microHH tools
execfile('../../python/microhh_tools.py')
filetoread='myscalar.0000010'

# read binary file
res=read_restart_file(filetoread,64,64,64,endian='little')

# write bubble
res[32][32][32]=1.0
print('Setting scalar=1 at 32,32,32. Total scalar=',sum(sum(sum(res))),'\n')

# rewrite binary file
write_restart_file(res,64,64,64,filetoread,'per_slice=True',endian='little')
res=read_restart_file(filetoread,64,64,64,endian='little')
print('Quick check: total scalar written=',sum(sum(sum(res))),'\n')

# modify restart.ini file
replace_namelist_value('starttime',10,'scalar_bubble_restart.ini')
replace_namelist_value('endtime',20,'scalar_bubble_restart.ini')
print('Restart namelist modified \n')
