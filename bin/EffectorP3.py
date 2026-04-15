#!/usr/bin/env python3
import sys
import subprocess
import os

if __name__ == "__main__":
    # Absolute path to the EffectorP3 installation directory
    EFFECTORP3_DIR = "/icarbon_pipp/predector/bin/EffectorP3"
    EFFECTORP_BIN = os.path.join(EFFECTORP3_DIR, "EffectorP.py")
    
    # Ensure the directory is in PYTHONPATH so it can find its companion functions.py
    env = os.environ.copy()
    env["PYTHONPATH"] = EFFECTORP3_DIR + os.pathsep + env.get("PYTHONPATH", "")
    
    # Execute the underlying script using the same Python interpreter
    sys.exit(subprocess.call([sys.executable, EFFECTORP_BIN] + sys.argv[1:], env=env))
