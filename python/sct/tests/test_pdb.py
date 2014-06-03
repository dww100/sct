
correct_lines = [
"ATOM     21  CB  LEU X   3      19.573  10.448-163.020  1.00  0.00           C\n"]

faulty_lines = [
"ATOMM     21  CB  LEU X   3      19.573  10.448-163.020  1.00  0.00           C\n",
"ATOM     21  CB  LEUUUU X   3      19.573  10.448-163.020  1.00  0.00           C\n",
"ATOM    21  CB  LEU X   3      19.573  10.448-163.020  1.00  0.00           C\n",
"ATOM     21  CB  LEU X   3          19.573  10.448-163.020  1.00  0.00           C\n"]

from sct.pdb import *

def test_pdb_res_line_parse_correct_lines():
    for l in correct_lines:
        try:
            pdb_res_line_parse(l,"TEST")
        except:
            print l
            assert False
        assert True


def test_pdb_res_line_parse_faulty_lines():
    for l in faulty_lines:
        try:
            pdb_res_line_parse(l,"TEST")
        except IOError:
            continue
        else:
            print l
            assert False
    assert True
        
