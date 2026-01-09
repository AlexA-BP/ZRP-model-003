import h5py as h5
from pathlib import Path
import numpy as np

def main():
    with h5.File(Path("data/test.h5"), "r") as fid:
        dset = np.array(fid["mocksim/particles"])

        print(dset)
    

if __name__ == "__main__":
    main()