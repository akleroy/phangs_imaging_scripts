# Staging script for NGC 5128

os.chdir('../ngc5128/')

os.system('rm -rf ngc5128_803_7m.ms.contsub')
uvcontsub(vis='ngc5128_803_7m.ms',
          field='',
          fitspw='0:100~1000;0:3000~3900',
          solint='int',
          fitorder=0,
          spw='0')

os.system('rm -rf ngc5128_803_12m.ms.contsub')
uvcontsub(vis='ngc5128_803_12m.ms',
          field='',
          fitspw='0:100~1000;0:3000~3900',
          solint='int',
          fitorder=0,
          spw='0')

line_ext = '.contsub'

cont_ext = ''
