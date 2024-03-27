### - Simple tyeps: These are types that you can use to annotate your code without importing anything special.
###     - Primitive types: int, float, bool, str, bytes, object
### - Any type: If you can't find a good type for some value, you can fall back to Any type. It's not a builtin type
### - Generic types (compelx types that are composed of simple types): list, tuple, dict, Iterable, Mapping, type

## 0. Simple types


def consume_many_types(
    num: int, decimal: float, boolean: bool, string: str, binary: bytes, obj: object
) -> None:
    ...

    ## try to type the parameter and ., you sill get all the suggestions methods available
    # obj.  ## un comment and check


## after writing the function name, press ctrl + space for the "auto completion".
# consume_many_types(

# )

### 1. Any type
from typing import Any

any_value: Any = (
    0  # you don't need to use this, you won't get any chance to use this also.
)


## 2. Generic types
nums: list[int] = [1, 2, 3, 4]  ## this is pronounced as list of ints
### sometimes `list[int]` this kind of format not supported in older version of the python, so solution for this is to import those.
from typing import List

nums: List[int] = [1, 2, 3, 4]


from typing import Tuple

three_dimensional_vector: Tuple[int, float, str] = (
    1,
    2.0,
    "3",
)  # now you will get the proper auto completion.
n_dimensional_vector: Tuple[int, ...] = (
    1,
    2,
    3,
    4,
    5,
)
# ... means syntax for defining the tuple in any way

from typing import (
    Dict,
    Mapping,
)

students_to_agents: Dict[str, int] = {
    "bobby": 25,
    "murph": 27,
    "alice": 21,
}

students_to_agents: Mapping[str, int] = (
    {  # mapping is used to repreent the dict in borader way
        "bobby": 25,
        "murph": 27,
        "alice": 21,
    }
)

## recommadable way is use the narrower way like Dict don't use broader way Mapping.

from typing import Sequence  # sequence implies that the data type is ordered
from typing import (
    Iterable,
    Set,
)

fruits: Set[str] = {
    "apple"
}  # you can get the auto completion easily for the Set perspective.
fruits: Iterable[str] = {
    "apple"
}  # you can get the auto completion easily for the Iterable perspective

## 3. Union (complex type)

from typing import (
    Type,
    Union,
)

miscellaneous_values: List[Union[int, float, str, Type]] = [1, 1.0, "hi", object]

## union means, whatever the value is going to be assummed by a variable of the union type, it allowed to be any of the types inside the Union list.
## If you specified Union[int, flaot, str], the value should contains any of the values from the list you have created inside the Union.
x: Union[int, float, str, Type] = [2.0, "super"]

## this is the example of complex type because it contias the composition of the types.
## Recent version of the python contains fancier version of the Union types by denoting pipe symbol
x: List[int | float | str | Type] = [2.0, "super"]  # it is lot cleaner compare to union


## 4.Optional
## - Optional keyword is essentially syntax sugar, or in other words, an abbreviated way of writing the type union, some type like int and then more.

from typing import Optional

x: Optional[int] = 1  # it might be int or None.


def greet(name: Optional[str] = None) -> str:  # pretty redable way
    return name


## another way of representing a optional is pipe symbol
x: str | None = "bro"
