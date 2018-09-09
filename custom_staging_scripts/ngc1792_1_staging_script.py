# Staging script for NGC 1792

os.chdir('../ngc1792/')

os.system('rm -rf ngc1792_1_886_7m_1.ms.contsub')
uvcontsub(vis='ngc1792_1_886_7m_1.ms',
          field='',
          fitspw='16:100~1000;16:3000~3900',
          solint='int',
          fitorder=0,
          spw='16')

os.system('rm -rf ngc1792_1_886_7m_2.ms.contsub')
uvcontsub(vis='ngc1792_1_886_7m_2.ms',
          field='',
          fitspw='16:100~1000;16:3000~3900',
          solint='int',
          fitorder=0,
          spw='16')

os.system('rm -rf ngc1792_1_886_7m_3.ms.contsub')
uvcontsub(vis='ngc1792_1_886_7m_3.ms',
          field='',
          fitspw='16:100~1000;16:3000~3900',
          solint='int',
          fitorder=0,
          spw='16')

os.system('rm -rf ngc1792_1_886_7m_4.ms.contsub')
uvcontsub(vis='ngc1792_1_886_7m_4.ms',
          field='',
          fitspw='16:100~1000;16:3000~3900',
          solint='int',
          fitorder=0,
          spw='16')

#os.system('rm -rf ngc1792_1_886_12m_1.ms.contsub')
#uvcontsub(vis='ngc1792_1_886_12m_1.ms',
#          field='',
#          fitspw='25:100~1000;25:3000~3900',
#          solint='int',
#          fitorder=0,
#          spw='25')

#os.system('rm -rf ngc1792_1_886_12m_2.ms.contsub')
#uvcontsub(vis='ngc1792_1_886_12m_2.ms',
#          field='',
#          fitspw='25:100~1000;25:3000~3900',
#          solint='int',
#          fitorder=0,
#          spw='25')

line_ext = '.contsub'

cont_ext = ''
