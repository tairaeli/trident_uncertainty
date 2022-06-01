import configparser, os, argparse

parser = argparse.ArgumentParser()

parser.add_argument('background_file', nargs = '?', action = 'store', default = os.path.expanduser("~/.trident/hm2012_ss_hr.h5"))

parser.add_argument('config_loc', nargs = '?', action = 'store', default = os.path.expanduser("~/.trident/config.tri"))

args = parser.parse_args()

config = configparser.ConfigParser()

config.read(args.config_loc)

config.set("Trident","ion_table_file",args.background_file)

