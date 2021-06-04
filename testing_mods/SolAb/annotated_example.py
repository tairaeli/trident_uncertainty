def test_SolAb(Data_Filename, proj_filepath, spec_filepath1, spec_filepath2, data_path = '~/git_env/research/oshea/trident_modifications/testing_mods/SolAb/enzo_cosmology_plus/RD0009/RD0009', line_list = ['H', 'C', 'O', 'Fe', 'Mg', 'Si'], **reading_func_args):
    import yt
    import trident
    fn = data_path
    print('LOADING COSMOLOGY PLUS')
    ds = yt.load(fn)
    print('SETTING RAY START AND END')
    ray_start = ds.domain_left_edge
    ray_end = ds.domain_right_edge
    if reading_func_args:
    	print('SETTING FIELD PARAMETERS')
    	field_params = {'reading_func_args':reading_func_args}
    	ray = trident.make_simple_ray(ds, start_position=ray_start, end_position=ray_end, 		data_filename = Data_Filename, lines = line_list, ftype = 'gas', field_parameters = 	field_params)
    else:
    	print('NO FIELD PARAMS')
    	ray = trident.make_simple_ray(ds, start_position=ray_start, end_position=ray_end, 		data_filename = Data_Filename, lines = line_list, ftype = 'gas')
    
    print('MAKING PROJECTION')	
    p = yt.ProjectionPlot(ds, 'x', 'density')
    print('ANNOTATING....')
    p.annotate_ray(ray, arrow=True)
    print('SAVING...')
    p.save(proj_filepath)
    print('ESTABLISHING SPECTRUM GENERATOR')
    sg = trident.SpectrumGenerator('COS-G130M')
    print('MAKING SPECTRUM')
    if reading_func_args:
    	sg.make_spectrum(ray, lines = line_list, abundance_table_args = reading_func_args)
    else:
    	sg.make_spectrum(ray, lines = line_list)
    print('SAVING SPECTRUM...')
    sg.save_spectrum(spec_filepath1)
    print('SPECTRUM COMPLETED')
    sg.plot_spectrum(spec_filepath2)
    
def setup(file = '~/git_env/research/oshea/trident_modifications/testing_mods/SolAb/cgm_abundances_2eb.txt', names = ['ray_cgm64.h5', 'proj_cgm64.png', 'spec_raw_cgm64.txt', 'spec_raw_cgm64.png']):
    import pandas as pd

    reading_func_args = {'filename': file, 'select_row': 0, 'ratios': False, 'kwargs': {'delim_whitespace': True}}
	
    df = pd.read_csv(file, delim_whitespace=True)
	
    cols = list(df.columns)
	
    test_SolAb(names[0], names[1], names[2], names[3], line_list = cols, **reading_func_args)
