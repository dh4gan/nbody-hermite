# Written 30/11/16 by dh4gan
# Reads in data from nbody_hermite (individual body files)
# Animates the result

import io_nbody_hermite as io


prefix = raw_input('What is the run prefix?')
prefix = prefix+'*'
io.animate_bodies(prefix)
