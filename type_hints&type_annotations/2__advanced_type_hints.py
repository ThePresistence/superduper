# Nested types: Dict 

from typing import Dict, Union, List 
student: Dict[Union[str], Union[str, int]]= {
    "name": "Marcy", 
    "age": 25
}

## The last type annotatoions is really right, let's re-itrate the "General principle of pipeline", it means type as specific as you possibly can with your type annotaions. 
## -- because the more specific you are, the more usefull type annotations are going to be in terms of "auto-completion". 
## Let's introdude the new tool to type more specific, the tool name is called "TypedDict" 

from typing import TypedDict 

class Student(TypedDict):  ## this is custom type 
    name: str
    age: int 

student: Student = {
    "name": "Marcy", 
    "age": 25
}

## you can also use the TypedDict dictionary to initiate the dictionary from scratch. 
other_student = Student(name="age", age="name")
print(other_student)

## There is a very similar concept is called "NamedTuple"
from collections import namedtuple

Point = namedtuple("Point", ["x", "y"])
point2d = Point(1, 2)
point2d.x 
point2d.y 

from typing import NamedTuple

class SuperTuple(NamedTuple):  ## This is custom type 
    x: int 
    y: int 

point2 = SuperTuple(3, 4)
point2.x
point2.y 
point2[0] 

## what are all the things we seen everything is a complex type: complex types are type which is a composition of similar types. 
## Example of complex type bellow: 

class Student(TypedDict):  ## this is custom type 
    """
    This is a composition of primitive types and the complex types SuperTuple 
    """
    name: str
    age: int 
    point = SuperTuple 
    


## This is the example of compelx types, also you can specify the "recursive types" where types are self refreshing. 

class Students(TypedDict):  ## this is custom type 
    """
    Recursive type and self referencive type example  
    """
    name: str
    age: int 
    point: SuperTuple 
    friends: List[Students]  ## Here it's giving error, because it doesn't know what the type it is !!! 
    ## So will change this to "self-referential types"
    friends: List["Students"]  ## Now it's getting, and this is the example of self referential types
    ## we used quotes here to reference a type that technically hasn't been defined by the way. We can do the similar things before we even run the class. 
    ## by importing like this you can do the self-referential without quotes 

    from __future__ import annotations 
    ## it's a special type of import in python, don't import until unless you need this specifically :) 
    ## it must be begining of the python file :) 

## Data Class :) 
    ## In general, it's really phenomenal that we can take dictionaries and we can assign types to them. 
    ## And suddenly get autocompletion on both the keys and the keys and the values being correct type in the dictionary. 
    ## But actually, if you don't need to use a dictionary, I would recommend using a data class instead. 

from dataclasses import dataclass

@dataclass
class NewDataClass: 
    name: str 
    age: int 
    bro: float 
    position: Point


bro = NewDataClass(**{"name": "aravind", "age": "superage"})

print(bro)









    









