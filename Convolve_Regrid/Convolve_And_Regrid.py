# -*- coding: utf-8 -*-
"""
File meant to do the convolve and regrid at various stages.
"""
from Convolve import master_convolve
from regrid import resample

# File to Convolve 
file = 'Sky_Remove/Median_Removed/LinearAvgd/Subtraction/AdjustMask_Left/3count/LinearSub_15_24micronV2_avg.fits'
f2 = 'Sky_Remove/Median_Removed/LinearAvgd/Subtraction/AdjustMask_Left/3count/LinearSub_15_70micron_new_avg.fits'
f3 = 'Sky_Remove/Median_Removed/Compared/Subtraction/sub_100micron_using_final70bkgd.fits'
# Kernels to Convolve with
kern70_160 = '../Kernels/Kernel_LoRes_PACS_70_to_PACS_160.fits'
kern24_160 = '../Kernels/Kernel_LoRes_MIPS_24_to_PACS_160.fits'
kern100_160 = '../Kernels/Kernel_LoRes_PACS_100_to_PACS_160.fits'

# Save Name for Convolution
sname_c24_160 = 'Sky_Remove/Median_Removed/Regrid_Convolved/Convolved_Only/final_sub_24_convto160.fits'
sname_c70_160 = 'Sky_Remove/Median_Removed/Regrid_Convolved/Convolved_Only/final_sub_70_convto160.fits'
sname_c100_160 = 'Sky_Remove/Median_Removed/Regrid_Convolved/Convolved_Only/final_sub_100_convto160.fits'
# Convolve
#c70_160 = master_convolve(kern70_160,f2,sname_c70_160)
#c24_160 = master_convolve(kern24_160,file,sname_c24_160)
c100_160 = master_convolve(kern100_160,f3,sname_c100_160)
print("Convolution Complete")
# Files to regrid to
f100 = 'Sky_Remove/Median_Removed/100um_medianRemoved.fits'
f160 = '../Original_Files/e0102_pacs160_new.fits'
# Save Name for Regrid

sname_r24_160 = 'Final/24um_subtraction_skyRemoved_linearAvgd_3count_leftExtend_15.fits'
sname_r70_160 = 'Final/24um_subtraction_skyRemoved_linearAvgd_3count_leftExtend_15.fits'
sname_r100_160 = 'Final/24um_subtraction_skyRemoved_modeled_with_70um.fits'
# Regrid
print("Regridding...")
#r24_100 = resample(sname_c24_100, f100, sname_r24_100)
#print('Regrid 100 file done')
#r24_160 = resample(sname_c24_160, f160, sname_r24_160)
#r70_160 = resample(sname_c70_160, f160, sname_r70_160)
r100_160 = resample(sname_c100_160, f160, sname_r100_160)
print('Regrid 160 file done')
print("Complete")
"""
# File to Convolve 
file70 = 'Sky_Remove/Median_Removed/LinearAvgd/Background/LinearBkgd_70micron_new_avg.fits'
# Kernels to Convolve with
kern70_100 = '../Kernels/Kernel_LoRes_PACS_70_to_PACS_100.fits'
kern70_160 = '../Kernels/Kernel_LoRes_PACS_70_to_PACS_160.fits'
# Save Name for Convolution
sname_c70_100 = 'Sky_Remove/Median_Removed/Regrid_Convolved/Convolved_Only/70um_bkgd_medianRemoved_convto100.fits'
sname_c70_160 = 'Sky_Remove/Median_Removed/Regrid_Convolved/Convolved_Only/70um_bkgd_medianRemoved_convto160.fits'
# Convolve
c24_100 = master_convolve(kern70_100,file70,sname_c70_100)
c24_160 = master_convolve(kern70_160,file70,sname_c70_160)
# Files to regrid to
f100 = 'Sky_Remove/Median_Removed/100um_medianRemoved.fits'
f160 = 'Sky_Remove/Median_Removed/160um_medianRemoved.fits'
# Save Name for Regrid
sname_r70_100 = 'Sky_Remove/Median_Removed/Regrid_Convolved/70um_bkgd_MedianRemoved_Regridded_Convolved_to_100.fits'
sname_r70_160 = 'Sky_Remove/Median_Removed/Regrid_Convolved/70um_bkgd_MedianRemoved_Regridded_Convolved_to_160.fits'
# Regrid
r70_100 = regrid(sname_r70_100, f100, sname_r70_100)
r70_160 = regrid(sname_r70_160, f160, sname_r70_160)



f = 'Sky_Remove/Comparison/Background/backgroundmod70_100.fits'
k = '../Kernels/Kernel_LoRes_PACS_100_to_PACS_160.fits'
cn = 'Sky_Remove/Regrid_Convolved/Convolved_Only/convolve100scaledto160.fits'
conv = master_convolve(k,f,cn)
orig_160 = 'Sky_Remove/Const_rm_160.fits'
rn = 'Sky_Remove/Regrid_Convolved/100_Convolved_to_160_Regrid.fits'
r = regrid(cn, orig_160, rn)
# Background Files to Convolve and Regrid
f24 = 'Sky_Remove/LinInterp_BKG_24micron_avgd_avg.fits'
f70 = 'Sky_Remove/LinInterp_BKG_70micron_new_avgd_avg.fits'
# Necessary Kernels
kern24_160 = '../Kernels/Kernel_LoRes_MIPS_24_to_PACS_160.fits'
kern24_100 = '../Kernels/Kernel_LoRes_MIPS_24_to_PACS_100.fits'
kern70_160 = '../Kernels/Kernel_LoRes_PACS_70_to_PACS_160.fits'
kern70_100 = '../Kernels/Kernel_LoRes_PACS_70_to_PACS_100.fits'
kern70_24 = '../Kernels/Kernel_LoRes_PACS_70_to_MIPS_24.fits'
# Save Name For Convolutions crm = constant removed
conv24_100 = 'Sky_Remove/Regrid_Convolved/Convolved_Only/convolve24crmto100.fits'
conv24_160 = 'Sky_Remove/Regrid_Convolved/Convolved_Only/convolve24crmto160.fits'
conv70_100 = 'Sky_Remove/Regrid_Convolved/Convolved_Only/convolve70crmto100.fits'
conv70_160 = 'Sky_Remove/Regrid_Convolved/Convolved_Only/convolve70crmto160.fits'
conv70_24 = 'Sky_Remove/Regrid_Convolved/Convolved_Only/convolve70crmto24.fits'
print("Files Imported")
# Convolution
c24_100 = master_convolve(kern24_100, f24, conv24_100)
c24_160 = master_convolve(kern24_160, f24, conv24_160)
c70_100 = master_convolve(kern70_100, f70, conv70_100)
c70_160 = master_convolve(kern70_160, f70, conv70_160)
c70_24 = master_convolve(kern70_24,f70,conv70_24)
print("Convolution Complete")
# Original Files to regrid to, with constant removed. 
orig_24 = 'Sky_Remove/Const_rm_24.fits'
orig_70 = 'Sky_Remove/Const_rm_70.fits'
orig_100 = 'Sky_Remove/Const_rm_100.fits'
orig_160 = 'Sky_Remove/Const_rm_160.fits'
# Save names for Regrids
fsave24_100 = 'Sky_Remove/Regrid_Convolved/24_Convolved_to_100_Regrid.fits'
fsave24_160 = 'Sky_Remove/Regrid_Convolved/24_Convolved_to_160_Regrid.fits'
fsave70_24 = 'Sky_Remove/Regrid_Convolved/70_Convolved_to_24_Regrid.fits'
fsave70_100 = 'Sky_Remove/Regrid_Convolved/70_Convolved_to_100_Regrid.fits'
fsave70_160 = 'Sky_Remove/Regrid_Convolved/70_Convolved_to_160_Regrid.fits'
# Regrid
r24_100 = regrid(conv24_100, orig_100, fsave24_100)
r24_160 = regrid(conv24_160, orig_160, fsave24_160)
r70_100 = regrid(conv70_100, orig_100, fsave70_100)
r70_160 = regrid(conv70_160, orig_160, fsave70_160)
r70_24 = regrid(conv70_24,orig_24,fsave70_24)
print("Regrids Complete")


f24 = 'Background/LinSpace/Background/24micron_avgd_avg.fits'
f70 = 'Background/LinSpace/Background/70micron_new_avgd_avg.fits'
kern24_160 = '../Kernels/Kernel_LoRes_MIPS_24_to_PACS_160.fits'
kern24_100 = '../Kernels/Kernel_LoRes_MIPS_24_to_PACS_100.fits'
kern70_24 = '../Kernels/Kernel_LoRes_PACS_70_to_MIPS_24.fits'

conv24_100 = 'Convolutions/newconvolve24to100.fits'
conv24_160 = 'Convolutions/newconvolve24to160.fits'
conv70_24 = 'Convolutions/convolve70to24.fits'

c24_100 = master_convolve(kern24_100, f24, conv24_100)
c24_160 = master_convolve(kern24_160, f24, conv24_160)
c70_24 = master_convolve(kern70_24,f70,conv70_24)

orig_24 = '../Original_Files/mosaic_24_n76.fits'
orig_100 = '../Original_Files/e0102_pacs100_new.fits'
orig_160 = '../Original_Files/e0102_pacs160_new.fits'

fsave_24 = '70convolvedto24_Regrid.fits'
fsave_100 = 'New_24convolvedto100_Regrid.fits'
fsave_160 = 'New_24convolvedto160_Regrid.fits'

r24_100 = regrid(conv24_100, orig_100, fsave_100)
r24_160 = regrid(conv24_160, orig_160, fsave_160)
r70_24 = regrid(conv70_24,orig_24,fsave_24)
"""