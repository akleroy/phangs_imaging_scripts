import glob

fname_list = [
    'ic5332_7m_co21.ms',
    'ic5332_12m_co21.ms',
    'ngc0628_7m_co21.ms',
    'ngc0628_12m_co21.ms',
    'ngc1087_7m_co21.ms',
    'ngc1300_7m_co21.ms',
    'ngc1385_7m_co21.ms',
    'ngc1433_7m_co21.ms',
    'ngc1512_7m_co21.ms',
    'ngc1566_7m_co21.ms',
    'ngc1672_7m_co21.ms',
    'ngc1672_12m_co21.ms',
    'ngc2835_7m_co21.ms',
    'ngc2835_12m_co21.ms',
    'ngc3351_7m_co21.ms',
    'ngc3351_12m_co21.ms',
    'ngc3627north_7m_co21.ms',
    'ngc3627north_12m_co21.ms',
    'ngc3627south_7m_co21.ms',
    'ngc3627south_12m_co21.ms',
    'ngc4254north_7m_co21.ms',
    'ngc4254north_12m_co21.ms',
    'ngc4254south_7m_co21.ms',
    'ngc4254south_12m_co21.ms',
    'ngc4303_7m_co21.ms',
    'ngc4303_12m_co21.ms',
    'ngc4321north_7m_co21.ms',
    'ngc4321north_12m_co21.ms',
    'ngc4321south_7m_co21.ms',
    'ngc4321south_12m_co21.ms',
    'ngc4535_7m_co21.ms',
    'ngc4535_12m_co21.ms',
    'ngc5068north_7m_co21.ms',
    'ngc5068north_12m_co21.ms',
    'ngc5068south_7m_co21.ms',
    'ngc5068south_12m_co21.ms',
    'ngc6744north_7m_co21.ms',
    'ngc6744north_12m_co21.ms',
    'ngc6744south_7m_co21.ms',
    'ngc6744south_12m_co21.ms',
    ]

f = open('baseline_stats.txt','w')
f.write('# column 1: data set')
f.write('# column 2: uv distace - minimum\n')
f.write('# column 3: uv distace - 20\n')
f.write('# column 4: uv distace - 30\n')
f.write('# column 5: uv distace - 50\n')
f.write('# column 6: uv distace - 75\n')
f.write('# column 7: uv distace - maximum\n')

for fname in fname_list:
    this_file = glob.glob('../*/'+fname)
    nbase, min_uv, max_uv, median_uv, mean_uv, rms_uv \
        , perc20_uv, perc25_uv, perc30_uv, perc75_uv, perc90_uv \
        = au.getBaselineStats(this_file[0])
    line = fname+' '+\
        str(nbase)+' '+\
        str(min_uv)+' '+\
        str(perc20_uv)+' '+\
        str(perc30_uv)+' '+\
        str(median_uv)+' '+\
        str(perc75_uv)+' '+\
        str(perc90_uv)+' '+\
        str(max_uv)
    f.write(line+'\n')
    
f.close()
