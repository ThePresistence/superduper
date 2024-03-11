# Type hints

## 1. History

* In python 3.5, a new feature was added to Python, that is nothing but "type annotations".

    ```python
    variable: int=25  # type annotations
    ```

* By using type annotations in your software project, the code that you write is going to be a lot more maintainable and it'll be far more likely people will actually want to adopt software libraries that you write.
* It's the largest impact of the code quality system, it proivides redability, maintainability and extendability.
* **What is typing?** Python is a dynamically typed language, this means that the python interpreter does type checking only as code runs, and the type of variable is allowded to change over its lifetime. It means dynamic typing we can change the type of the variable anytime, however **static typing** does not allow us to change the type of the variable. Here the **process of type checking is called static type** checking where types are actually determined before a program is executed.
* **FastAPI** is the new productionized web framework for flask. It has a very good auto completion.
* One great advantage of **type annotations**, it can save you apparently 15 percent of debugging time at least. It will give you the red squiggly line, if you're typing wrong types.
* **Auto completion** is the one of the biggest problem in **Python**.
* Whenver you're writing a code you should follow the principle **the principle of least surprise**. This states that, if there is a problem, there's probably multiple implementation you could chosse you should choose one that is more obivious like the soultion most expereienced developers choose.
* This is the first [docmentation](https://peps.python.org/pep-0484/) for the type hints by the python.
* One of the advantage of adding type annotations in your code, it allows you to do **static analysis** meaning we can analyze the type annotations in our python code without even running the program and predict errors before they happen.

## 2. Mypy

* Mypy is used to do the **static code analysis**, it's a another cli tool that allows us to statically analyze our code and in this case the function of mypy is that it analyzes our type annotations and evaluvate whether or we are consistently and correctly using variables according to the types that we assing them.
* Mypy is the default type checker for **python programming**.

    ```python
    pip install mypy
    mypy <file-name>
    ```

* They also having a **vs code extensions**. It will give you the red squiggly line after installing.
* There are different kinds of type checking tools **Wasabi** from facebook (for larger codebase), **Pyre** from google, and **pyrite** from microsoft.

## 3. Basic typing

* Type inference is the idea that youc can guess the type of a variable by watching what values or what types of values the variable has taken on over the course of a program in various in your code.
* You can have two different type declaration for a variable and that can be **implicit type** and **explicit type**. 
* Explicit type: ```nums: List[int] = [2, 1]```. 
* Inferred type: ```nums: [2, 1]```, if we remove the explicit type, that is implicit types.