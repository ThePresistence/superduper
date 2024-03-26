import sys, os 
sys.path.append(
    os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir)))

from ignore_files.main import add_bro, minus_bro 
import pytest 


## parametrized testing 
@pytest.mark.parametrize( 
        argnames="a, b, out", 
        argvalues=[
            (7,7, 14), 
            (8,8, 16), 
            (9,9, 18)
        ], 
    
)
def test__bro_done(a:int, b:int, out:int): 
    assert add_bro(a, b) == out

# markers 
@pytest.mark.slow 
def test__bro_minus(): 
    assert minus_bro(7, 8) == 1

## If you dont' want to exectue this you can use this command 
## pytest -m "not slow"
## If you want to only run the tests are slow, you the bellow command 
## pytest -m "slow"

@pytest.mark.fast 
def test___bro_super_duper(): 
    assert minus_bro(7, 8) == 1 

