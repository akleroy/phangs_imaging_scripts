# This is the line list.

# Drawn from Splatalogue at http://www.cv.nrao.edu/php/splat/

# Lists that combine multiple transitions for a single line

line_families = {
    'co':['co10','co21','co32','co43','co54','co65'],
    '13co':['13co10','13co21','13co32','13co43','13co54','13co65'],
    'c18o':['c18o10','c18o21','c18o32','c18o43','c18o54','c18o65'],
    'hcn':['hcn10','hcn21','hcn32','hcn43','hcn54','hcn65','hcn76'],
    'h13cn':['h13cn10','h13cn21','h13cn32','h13cn43','h13cn54','h13cn65','h13cn76'],
    'hnc':['hnc10','hnc21','hnc32','hnc43','hnc54','hnc65','hnc76'],
    'hn13c':['hn13c10','hn13c21','hn13c32','hn13c43','hn13c54','hn13c65','hn13c76'],
    'hcop':['hcop10','hcop21','hcop32','hcop43','hcop54','hcop65','hcop76'],
    'h13cop':['h13cop10','h13cop21','h13cop32','h13cop43','h13cop54','h13cop65','h13cop76'],
    'cs':['cs10','cs21','cs32','cs43','cs54','cs65','cs76','cs87','cs98','cs109','cs1110','cs1211','cs1312','cs1413'],
    '13cs':['13cs10','13cs21','13cs32','13cs43','13cs54','13cs65','13cs76','13cs87','13cs98','13cs109','13cs1110','13cs1211','13cs1312','13cs1413'],
    'sio':['sio10','sio21','sio32','sio43','sio54','sio65','sio76','sio87','sio98','sio109','sio1110','sio1211','sio1312','sio1413','sio1514','sio1615'],
    'hi':['hi21cm'],
    'ci':['ci10','ci21'],
    }

# The line list dictionary

line_list = {
    'co65':691.47308,
    'co54':576.26793,
    'co43':461.04077,
    'co32':345.79599,
    'co21':230.53800,
    'co10':115.27120,
    '13co65':661.06728,
    '13co54':550.92629,
    '13co43':440.76517,
    '13co32':330.58797,
    '13co21':220.39868,
    '13co10':110.20135,
    'c18o65':658.55328,
    'c18o54':548.83101,
    'c18o43':439.08877,
    'c18o32':329.33055,
    'c18o21':219.56035,
    'c18o10':109.78217,
    'hcn10':88.63185, # J=1-0, F=2-1
    'hcn21':177.26111, # J=2-1, F=2-1
    'hcn32':265.88618,
    'hcn43':354.50548,
    'hcn54':443.11616,
    'hcn65':531.71639,
    'hcn76':620.30410,
    'h13cn10':86.33992140,
    'h13cn21':172.67785120,
    'h13cn32':259.01179760,
    'h13cn43':345.33976930,
    'h13cn54':431.65977480,
    'h13cn65':517.96982100,
    'h13cn76':604.26791400,
    'cs10':48.99095,
    'cs21':97.98095,
    'cs32':146.96903,
    'cs43':195.95421,
    'cs54':244.93556,
    'cs65':293.91209,
    'cs76':342.88285,
    'cs87':391.84689,
    'cs98':440.80323,
    'cs109':489.75092,
    'cs1110':538.68900,
    'cs1211':587.61649,
    'cs1312':636.53246,
    'cs1413':685.43592,
    '13cs10':46.24756320,
    '13cs21':92.49430800,
    '13cs32':138.73933500,
    '13cs43':184.98177200,
    '13cs54':231.22068520,
    '13cs65':277.45540500,
    '13cs76':323.68497300,
    '13cs87':369.90855050,
    '13cs98':416.12527510,
    '13cs109':462.33429010,
    '13cs1110':508.53473910,
    '13cs1211':554.72576570,
    '13cs1312':600.90648000,
    '13cs1413':647.07615000,
    'hcop10':89.18852,
    'hcop21':178.37506,
    'hcop32':267.55763,
    'hcop43':356.73422,
    'hcop54':445.90287,
    'hcop65':535.06158,
    'hcop76':624.20836,
    'h13cop10':86.75428840,
    'h13cop21':173.50670030,
    'h13cop32':260.25533900,
    'h13cop43':346.99834400,
    'h13cop54':433.73383270,
    'h13cop65':520.45988430,
    'h13cop76':607.17464560,
    'hnc10':90.66357,
    'hnc21':181.32476,
    'hnc32':271.98114,
    'hnc43':362.63030,
    'hnc54':453.26992,
    'hnc65':543.89755,
    'hnc76':634.51083,
    'hn13c10':87.09085000,
    'hn13c21':174.17940800,
    'hn13c32':261.26331010,
    'hn13c43':348.34026950,
    'hn13c54':435.40796260,
    'hn13c65':522.46407300,
    'hn13c76':609.50628400,
    'ci10':492.16065, # 3P1-3P0
    'ci21':809.34197, # 3P2-3P1
    'sio10':43.42376,
    'sio21':86.84696,
    'sio32':130.26861,
    'sio43':173.68831,
    'sio54':217.10498,
    'sio65':260.51802,
    'sio76':303.92696,
    'sio87':347.33063,
    'sio98':390.72845,
    'sio109':434.11955,
    'sio1110':477.50310,
    'sio1211':520.87820,
    'sio1312':564.24396,
    'sio1413':607.59942,
    'sio1514':650.94359,
    'sio1615':694.27543,
    'hi21cm':1.420405751,
    }

# Run some consistency checks

def run_checks():
    """
    """
    all_okay = True
    
    for family in line_families:
        this_list = line_families[family]
        for this_line in this_list:
            if this_line not in line_list.keys():
                print("Line missing from line list but in line families: "+this_line)
                all_okay = False
    
    if all_okay:
        print("All lines in line families present in line list.")

    no_repeats = True

    for this_line in line_list:
        for other_line in line_list:
            if this_line == other_line:
                continue
            if line_list[this_line] == line_list[other_line]:
                print("Duplicate frequencies for: "+this_line+" and "+other_line+" . Check for typos.")
                no_repeats = False

    if no_repeats:
        print("No repeat frequencies in list.")
