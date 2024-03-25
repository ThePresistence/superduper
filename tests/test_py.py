import sys, os 
sys.path.append(
    os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir)))

from ignore_files.main import add_bro 

def test_bro_done(): 
    assert add_bro(7, 7) == 14 


