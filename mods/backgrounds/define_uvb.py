import configparser, os, argparse

parser = argparse.ArgumentParser()

parser.add_argument('--background_file', nargs = '?', action = 'store', dest = 'background_file', default = os.path.expanduser("~/.trident/hm2012_hr.h5"))

parser.add_argument('--config_loc', nargs = '?', action = 'store', dest = 'config_loc', default = os.path.expanduser("~/.trident/config.tri"))

args = parser.parse_args()

config = configparser.ConfigParser()

config.read(args.config_loc)

print(config.sections())

config.set("Trident","ion_table_file",args.background_file)

with open(args.config_loc, 'w') as new_config:
    config.write(new_config)

#####################################
config = configparser.ConfigParser()

config.read(args.config_loc)

print(config["Trident"]["ion_table_file"])