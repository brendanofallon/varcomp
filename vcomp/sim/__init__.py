import sys
import os

parent_dir = os.sep.join( os.path.dirname(__file__).split(os.sep)[0:-1])
sys.path.append( os.sep.join([parent_dir, "util"]) )

