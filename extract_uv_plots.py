# Extract u vs v and uv distance vs amp. from the central channels of
# the CO 2-1 line. This is both klugy and pretty specific to the
# PHANGS data.

try:
    input_vis
except NameError:
    print "Please specify an input measurement set as input_vis."

delta = 20

# Work out number of channels
vm = au.ValueMapping(input_vis)
nchan = vm.spwInfo[0]['numChannels']
middle_channel = round(nchan/2)

# Call plotms
plotms(vis=input_vis,
       xaxis='u',
       yaxis='v',
       spw='0:'+str(middle_channel),
       plotfile=input_vis+'.uv.txt',
       expformat='txt',
       showgui=False)
       
plotms(vis=input_vis,
       xaxis='uvdist',
       yaxis='amp',
       spw='0:'+str(middle_channel-delta/2)+'~'+str(middle_channel+delta/2),
       avgchannel=10000,
       plotfile=input_vis+'.amp_vs_uv.txt',
       expformat='txt',
       showgui=False)
