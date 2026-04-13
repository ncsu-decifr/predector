#!/usr/bin/env python3
import sys
import subprocess
import os

if __name__ == "__main__":
    effector_bin = os.path.join(os.path.dirname(os.path.realpath(__file__)), "EffectorP3", "EffectorP.py")
    sys.exit(subprocess.call([sys.executable, effector_bin] + sys.argv[1:]))
