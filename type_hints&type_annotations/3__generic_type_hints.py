## Generics 
## - It's basically a alias type in python by storing them as a variable.

## Example 
Tsting = str 
x: Tsting = "hi"  # This is generics ( storing the alias type in variable )

## - This is mainly used for stroring the complex types like bellow 
from typing import Dict, List, Union 

TStudentDict = Dict[str, Union[str, int, float]]
super_dict: TStudentDict = ...  # more complex example :)  


### Story Time 1.1
def add(a, b): 
    return a+b  
## This function is the one of the most complicated thing for type hinting because any datatype can go inside like "list","str", "float" 
## LEt's say you're storing like this: 
st = Union[str, float, list, int, tuple]

### Story Time 1.2 
def add(a:st, b:st): 
    return a + b   ## This function is more complicated we don't know what comes inside to solve thi we can use "Typevar types" 

## Typevar: what if we could make an alias that was a stand in for a type and gave us the ability to do this restriction.It makes sure the the data type are same :) 

from typing import TypeVar 
st_var = TypeVar("st_var")  ## give the variable name inside 

def adds(a:st_var, b:st_var) -> List[st_var]: 
    return [a, b]  ## now we are making sure both the types are same type 


adds(a="bro", b="super")











