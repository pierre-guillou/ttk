#!/usr/bin/env pvpython

import argparse

import paraview.simple
import paraview.smstate


def main(stateFile):
    paraview.simple.LoadState(stateFile)

    pyfilename = f'{stateFile.split(".")[0]}.py'

    with open(pyfilename, "w") as dst:
        dst.write(paraview.smstate.get_state())


if __name__ == "__main__":
    desc = "Converts a ParaView state file into a Python script"
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("stateFile", type=str, help="Path to ParaView state file")
    args = parser.parse_args()

    main(args.stateFile)
