import configparser, os, argparse

parser = argparse.ArgumentParser()

parser.add_argument('--uv_back_filename', nargs = '?', action = 'store', dest = 'uv_back_filename', default = "hm2012_hr.h5")
parser.add_argument('--uv_back_dir', nargs = '?', action = 'store', dest = 'uv_back_dir', default = os.path.expanduser("~/.trident"))
parser.add_argument('--config_loc', nargs = '?', action = 'store', dest = 'config_loc', default = os.path.expanduser("~/.trident/config.tri"))

args = parser.parse_args()

config = configparser.ConfigParser()

config.read(args.config_loc)

config.set("Trident","ion_table_file",args.uv_back_filename)

config.set("Trident","ion_table_dir",args.uv_back_dir)

with open(args.config_loc, 'w') as new_config:
    config.write(new_config)

#####################################
config = configparser.ConfigParser()

config.read(args.config_loc)

print(config["Trident"]["ion_table_file"])