{
 "cells": [
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "388e0ea8-7736-4f25-ad20-ccb05bd27bca",
   "metadata": {},
   "outputs": [],
   "source": [
    "from hpcc_sal_the_snake import *\n",
    "import numpy as np\n",
    "\n",
    "print(\"let's do some math, kids\")\n",
    "\n",
    "path = '~/test_sal/test1/'\n",
    "\n",
    "def visualize(ds_file, center_list, ray_dir, ray_num, ion='O VI', name='example_multiplot', num_dense_min = 1e-11, num_dense_max=1e-5, **vis_args):\n",
    "\n",
    "    \"\"\"\n",
    "    Uses salsa.AbsorberPlotter() to generate multiplots of data produced by salsa.AbsorberExtractor()\n",
    "\n",
    "    :ds_file: Halo dataset used in generating data\n",
    "\n",
    "    :center_list: x, y, and z coordinates of center of galaxy on plot\n",
    "\n",
    "    :ray_dir: Directory where ray.h5 files are stored (same one that's passed to sal())\n",
    "\n",
    "    :ray_num: Same as ray_num from run_sal(). Here, it is used to index ray.h5 file for information\n",
    "\n",
    "    :name: String used to name multiplot file\n",
    "\n",
    "    :num_dense_min: To be passed to salsa.AbsorberPlotter(). Default 1e-11\n",
    "\n",
    "    :num_dense_max: To be passed to salsa.AbsorberPlotter(). Default 1e-5\n",
    "\n",
    "    :vis_args: Mainly used by run_sal() to pass a multiplot name other than \"example_multiplot.png\". Not meant to be manually passed by user. \n",
    "    \"\"\"\n",
    "\n",
    "    if len(str(ray_num)) == 1:\n",
    "        new_ray_num = f'0{ray_num}'\n",
    "    else:\n",
    "        new_ray_num=ray_num\n",
    "\n",
    "    newname = f'{name}_ray{ray_num}.png'\n",
    "\n",
    "    ray = yt.load(f'{ray_dir}/ray{new_ray_num}.h5')\n",
    "    plotter = salsa.AbsorberPlotter(ds_file, ray, ion, center_gal=center_list, use_spectacle=False, plot_spectacle=False, plot_spice=True, num_dense_max=num_dense_max, num_dense_min=num_dense_min)\n",
    "\n",
    "    fig, axes = plotter.create_multi_plot(outfname=newname)\n",
    "\n",
    "\n",
    "\n",
    "def update_ionlist(list_name, change_index, ion):\n",
    "\n",
    "    \"\"\"\n",
    "    Used to more efficiently generate lists of ions to be passed to run_sal() (later passed to sal())\n",
    "\n",
    "    :list_name: Original list to be updated\n",
    "\n",
    "    :change_index: Integer -- index of ion to be replaced\n",
    "\n",
    "    :ion: Ion that will replace the old ion\n",
    "    \"\"\"\n",
    "\n",
    "    new_list = list_name.copy()\n",
    "    new_list[change_index] = ion\n",
    "\n",
    "    return new_list\n",
    "\n",
    "\n",
    "\n",
    "def generate_names(length, add='_'):\n",
    "\n",
    "    \"\"\"\n",
    "    Returns a list of generic names for multiplots anda list of generic names for raw data. These lists are typically passed to run_sal()\n",
    "\n",
    "    :length: length of lists to be generated. Do not account for indexing starting at zero.\n",
    "\n",
    "    :add: Additional relevant information for keeping track of multiplots and data later on. Default add=''\n",
    "    \"\"\"\n",
    "\n",
    "    vis_name_list = []\n",
    "    saved_filename_list = []\n",
    "\n",
    "    for i in range(length):\n",
    "        vis_name_list.append(f'multiplot_row{i}{add}')\n",
    "        saved_filename_list.append(f'data_row{i}{add}')\n",
    "\n",
    "    return vis_name_list, saved_filename_list\n",
    "\n",
    "def run_sal(vis_name, saved_filename, vis_tf, ray_dir, ray_num, path, n_rays, vis_add='_', saved_add = '_', **kwargs):\n",
    "\n",
    "    \"\"\"\n",
    "    Calls sal() and visualize() functions and saves data. \n",
    "\n",
    "    :vis_name: Indexed from vis_name_list and used to name generated multiplot of data according to iteration number (and other naming characteristics which are passed through vis_add)\n",
    "\n",
    "    :saved_filename: Indexed from saved_filename_list and used to name saved data according to iteration number (and other naming characteristics which are pass through saved_add)\n",
    "\n",
    "    :vis_tf: If True, uses visualize() to generate multiplots to accompany data\n",
    "\n",
    "    :ray_num: Used to keep track of iteration value from for loop in run_sal() to aid in naming multiplots and data\n",
    "\n",
    "    :path: Path to directories where data is stored. Default path=''\n",
    "\n",
    "    :vis_add: Used to add relevant information to the names of generated multiplots. Should be passed through **kwargs.\n",
    "\n",
    "    :saved_add: Used to add relevant information to the names of generated data. Should be passed through **kwargs.\n",
    "\n",
    "    :kwargs: See description from run_sal()\n",
    "    \"\"\"\n",
    "    print(f'KWARGS: {kwargs}')\n",
    "\n",
    "\n",
    "\n",
    "    catalog = sal(ray_dir=ray_dir, ray_num=ray_num, n_rays = n_rays, df_type='multiple', **kwargs)\n",
    "\n",
    "    new_saved_filename = saved_add+saved_filename\n",
    "    catalog.to_csv(f'{path}{new_saved_filename}.txt', sep = ' ')\n",
    "    catalog.to_csv(f'{path}{new_saved_filename}.csv', sep = ' ')\n",
    "\n",
    "    if vis_tf == True:\n",
    "        new_vis_name = vis_add+vis_name\n",
    "        vis_args = dict(name = f'{path}{new_vis_name}')\n",
    "        for r in range(n_rays):\n",
    "            visualize(ds_file='/mnt/research//galaxies-REU/sims/FOGGIE/halo_002392/nref11c_nref9f/RD0020', center_list=[0.53, 0.53, 0.53], ray_dir=ray_dir, ray_num=r, **vis_args)\n",
    "\n",
    "    print('go look at your data!')\n",
    "\n",
    "\n",
    "\n",
    "list1 = ['Ne VIII', 'Mg X', 'O VI', 'S IV', 'Si III', 'C II', 'N I']\n",
    "#list2 = update_ionlist(list1, 5, 'Fe II')\n",
    "#list3 = update_ionlist(list2, 3, 'Si IV')\n",
    "#list4 = update_ionlist(list1, 3, 'Si IV')\n",
    "\n",
    "#ions = np.array([[list1], [list2], [list3], [list4]])\n",
    "\n",
    "kwargs = dict(reading_func_args=dict(ray_dir=f'{path}rays', ion_list = list1)\n",
    "\n",
    "vis, saved = generate_names(25)\n",
    "\n",
    "selr = 0\n",
    "\n",
    "for i in range(25):\n",
    "    #kwargs['ion_list'] = list(ions[i][0])\n",
    "#     kwargs['reading_func_args']['select_row'] = selr\n",
    "    run_sal(vis[i], saved[i], vis_tf=False, ray_num=i, path=path, n_rays=10, saved_add=f'data/', **kwargs)\n",
    "#     selr += 1"
   ]
  }
 ],
 "metadata": {
  "kernelspec": {
   "display_name": "Python 3",
   "language": "python",
   "name": "python3"
  },
  "language_info": {
   "codemirror_mode": {
    "name": "ipython",
    "version": 3
   },
   "file_extension": ".py",
   "mimetype": "text/x-python",
   "name": "python",
   "nbconvert_exporter": "python",
   "pygments_lexer": "ipython3",
   "version": "3.8.8"
  }
 },
 "nbformat": 4,
 "nbformat_minor": 5
}
