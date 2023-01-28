# !/usr/bin/env python3
import os

# This is the regression test base module used for all test cases.

class SphinxsysRegressionTest:

    def __init__(self, casename, bodyname, parametername):
        self.sphinxsys_exec_path = os.path.abspath(os.path.join(os.getcwd()))
        self.sphinxsys_case_path = os.path.abspath(os.path.join(self.sphinxsys_exec_path))
        self.sphinxsys_case_name = casename
        self.sphinxsys_body_name = bodyname
        self.sphinxsys_parameter_name = parametername
        self.enter_sphinxsys_exec_folder = f"cd {self.sphinxsys_exec_path};"
        self.enter_sphinxsys_case_folder = f"cd {self.sphinxsys_case_path};"
        self.input_file_path = os.path.join(self.sphinxsys_exec_path, "input")
        self.condition_file_path = os.path.join(self.input_file_path, f"{bodyname}_{parametername}_runtimes.dat")

    def compile_case(self) -> None:
        print('Start compiling test case....')
        command = "make -j8"
        os.system(self.enter_sphinxsys_case_folder)
        os.system(command)
        print('Compiling test case is finished...')

    def run_particle_relaxation(self) -> None:
        print('Start particle relaxation for the simulation...')
        command = f".{os.sep}{self.sphinxsys_case_name} --r=true"
        os.system(self.enter_sphinxsys_exec_folder)
        os.system(command)
        print('Simulating case is finished...')

    def run_case(self) -> None:
        print('Start case simulation...')
        print(self.enter_sphinxsys_exec_folder)
        command = f".{os.sep}{self.sphinxsys_case_name} --r=false --i=true --rt=true"
        os.system(self.enter_sphinxsys_exec_folder)
        os.system(command)
        print('Simulating case is finished...')

    def read_dat_file(self):
        file = open(self.condition_file_path)
        ifconverged = file.readline(4)
        file.close()
        return ifconverged
    