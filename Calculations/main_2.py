import numpy as np
from mymain import my_main_1, preparation

eventfiles, background, scale, x_range = preparation()
eventfile_czt = eventfiles[0]
scale_czt = scale[0]
x_range_czt = x_range[0]

eventfile_coating = eventfiles[1]
scale_coating = scale[1]
x_range_coating = x_range[1]

eventfile_plastic = eventfiles[2]
scale_plastic = scale[2]
x_range_plastic = x_range[2]


all_contrib_at116Cd_coating, all_contrib_at130Te_coating = my_main_2(eventfile_coating, scale_coating, x_range_coating, background[0], Q_val_ret=True)
my_main_2(eventfile_czt, scale_czt, x_range_czt, background[1])
sum_glyp_at116Cd = all_contrib_at116Cd_coating[0] + all_contrib_at116Cd_coating[2] + all_contrib_at116Cd_coating[4]
sum_epox_at116Cd = all_contrib_at116Cd_coating[1] + all_contrib_at116Cd_coating[3] + all_contrib_at116Cd_coating[5]
sum_glyp_at130Te = all_contrib_at130Te_coating[0] + all_contrib_at130Te_coating[2] + all_contrib_at130Te_coating[4]
sum_epox_at130Te = all_contrib_at130Te_coating[1] + all_contrib_at130Te_coating[3] + all_contrib_at130Te_coating[5]
sum_coating_at_116Cd = sum_glyp_at116Cd + sum_epox_at116Cd
sum_coating_at_130Te= sum_glyp_at130Te + sum_epox_at130Te
new_data = np.array([['glyptal', 'epoxy', 'sum'], [sum_epox_at116Cd, sum_glyp_at116Cd, sum_coating_at_116Cd], [sum_epox_at130Te, sum_glyp_at130Te, sum_coating_at_130Te]])
descriptions = ['Contributions at Qvalues', 'N in 1/kg/keV/yr ']
var_names = ['Isotope', 'N_at116Cd', 'N_at130Te']
write_detailed_txtfile(np.transpose(new_data), var_names, descriptions, './calc_solutions/', 'events_at_Qvalues_allcoating.txt' )
