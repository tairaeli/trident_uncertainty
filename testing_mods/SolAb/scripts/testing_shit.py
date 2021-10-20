const = ['this/shit/andstuff/']
const.append(input('give me a string\n'))
bitch = ''.join(const)
print(bitch)

preliminary_dic_str = input("provide a dictionary (or type 'None' to enter manually) with the directory where the data and rays will be stored as 'path', the number of lightrays to be generated as 'nrays', the halo data to be uused as 'ds_file', and the abundances ('None' if using SolAb) to be used as 'file_path'.")
if preliminary_dic_str == 'None':
	path = input("Enter path to directory where the data and rays will be stored (must end in /): \n")
	print(f"PATH: {path}")
	nrays = input("Enter the number of lightrays to be generated: \n")
	nrays = int(nrays)
	ds_file = input("Enter the path to the halo data: \n")
	print(f"TYPE NRAYS: {type(nrays)}")
	file_path = input("Enter the path to the abundance file ('None' if using solar abundances): \n")
	if file_path == 'None':
		nrows = nrays
		litty = 'False'
	else:
		df = pd.read_csv(file_path, delim_whitespace=True)
		nrows = len(df)
		litty = 'True'
		
	print(f"NROWS: {nrows}")
else:
	preliminary_dic = json.loads(preliminary_dic_str)
	path = preliminary_dic['path']
	nrays = preliminary_dic['nrays']
	ds_file = preliminary_dic['ds_file']
	file_path = preliminary_dic['file_path']
