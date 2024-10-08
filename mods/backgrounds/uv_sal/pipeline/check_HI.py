"""
Recreate the Trident ion fraction field calculation for HI,
which already exists in the dataset 
and so Trident won't calculate it.
"""
import yt
import numpy as np
from os import path
from trident.ion_balance import IonBalanceTable
from yt.utilities.linear_interpolators import \
    TrilinearFieldInterpolator
from yt.frontends.ytdata.utilities import \
    save_as_dataset

# Define field functions needed for interpolation
def _log_nH(field, data):
    return np.log10(data["gas", "H_nuclei_density"])

def _log_T(field, data):
    return np.log10(data["gas","temperature"])

def _redshift(field, data):
    try:
        current_redshift = data.ds.current_redshift
    except AttributeError:
        current_redshift = 0.
    return current_redshift * \
        np.ones(data["gas", "density"].shape, 
                dtype=data["gas", "density"].dtype)

def spoof_ray_HI(saved_ray_filename, ion_table, out_path, n):

    # Load ion balance table
    tab = IonBalanceTable(filename=ion_table, atom="H")

    # Prepare dataset for interpolation
    ds = yt.load(saved_ray_filename)

    if ("gas", "log_nH") not in ds.derived_field_list:
        ds.add_field(("gas", "log_nH"), function=_log_nH, units="",
                    sampling_type="cell")
        
    if ("gas", "redshift") not in ds.derived_field_list:
        ds.add_field(("gas", "redshift"), function=_redshift, units="",
                    sampling_type="cell")

    if ("gas", "log_T") not in ds.derived_field_list:
        ds.add_field(("gas", "log_T"), function=_log_T, units="",
                    sampling_type="cell")

    # Interpolate
    ionFraction = tab.ion_fraction[0] # neutral
    n_param = tab.parameters[0]
    z_param = tab.parameters[1]
    t_param = tab.parameters[2]
    bds = [n_param.astype("=f8"), z_param.astype("=f8"), t_param.astype("=f8")]

    interp = TrilinearFieldInterpolator(ionFraction, bds,
                                        [("gas", "log_nH"),
                                        ("gas", "redshift"),
                                        ("gas", "log_T")],
                                        truncate=True)

    ad = ds.all_data()
    fraction = np.power(10, interp(ad))
    H_p0_number_density = fraction * ad["gas","H_nuclei_density"]

    # Save a copy of the ray with new HI field
    data = {key:ad[key] for key in ds.field_list if key[0]=='grid'}
    data[('grid','H_p0_number_density')] = H_p0_number_density

    field_types = dict([(field, "grid") for field in data.keys()])

    # doing some weird path-name stuff

    ray_basename = saved_ray_filename[-(6+n):]
    ray_num = ray_basename[3:3+n]

    yt.save_as_dataset(ds,
                       filename=out_path+"/ray"+ray_num+".h5",
                       data=data,
                       field_types=field_types,
                       extra_attrs={"data_type":"yt_light_ray"})
    
# if __name__ == "__main__":

#     tab_name = path.expanduser("~/.trident/hm2012_ss_hr.h5")
#     dataset = path.expanduser("~/research/ray.h5")

#     spoof_ray_HI(dataset, tab_name)