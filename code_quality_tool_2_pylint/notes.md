# Write clean and efficient code :) 

## 1. Pylint

* Black is a auto formatter, it helps to format space and identation however **pylint** instead of looking space and identation it will give some **qualitative recommedation**.
* Pylint is a static code analyzer, which means it's a tool that's able to treat your code as a bunch of text files and detect quality issues with your code **without actually executing your code**.
* Incase your code is slow, it doesn't affect **pytlint** ability to check your code.
* Pylint looks for erros and deviations from coding standards like pep8, and it also it catches the **code smells**.
* **Code smell** is a stylistic issue with your code that doesn't provide your code from executing or working properly, may indicate the serious problem with the structure of your code. To fix the code smell we need to do refractoring, pylint will give you the suggestions for the refactoring.
* In simple words, <i>coding patterns that indicate that something is wrong with the design of a programme</i>.

    ```bash
    pytlint <file-name.py>
    ```

## 2. Flake8

* Pylint: Linting tool is a CLI tool that can do some level of against your code. Linting tools usually have a set of rules that your code is checked against and any time you violate one of those rules.
* Pytlint is not only a linting tool in python space, there is a another popular tool called **Flake8**.
* **Flake8** is yet another CLI tool that has its own rules set that we can check our code against. PYlint and Flake8 has somewhat differing feautre set in the rules to enforce.
* Flake8 supports a rich plugin system, and because of the large number of community contributed plugins have been created for FLAKE8 and those plugins enhance **Flake8** functionality in all kinds of differnet ways.
* You can check the **flake8 plugins** in google.

    ```python
    flake8 <file-name.py> or <folder-name>
    ```

## 3. isort

* It's a another code quality tool to organize our import statements.
* If we are not keeping your import statements organized, in the best case, the redability of all our import statements is not going to be great. And in worst case, we will have a bug in our code.
* isort sorts your import statements.
* isort basically sorts the import statmenet based on the types of modules. Let's say if it's a standard package it will import the standard package first then the third party package like "numpy", "pandas" and so more. 
* There are lot of config you can do, you can check that in the documentation.

    ```python
    isort <file-name.py>
    ```

## 4. Radon

* It's a another set of command line tool, that validate the quality of the code and in this case, the quality check that we're going to be applying about **code complexity**.
* The important definition of the clean code is code that's easy to change to extend, because if your code isn't fast or secure or redable in some cases, as long as your code is easy to change, you can iterate towards making your code be those things. Something that makes code hard to change if the code is highly complex. If the code is higly complex, that will make it difficult to read and therefore reason about.
* We've a heuristic principle for clean code and it's called **the principle of least surprise**.
* It would be amazing if there were a perfect code metric that could tell us exactly how maintainable or redable a particular block of code is, because if that were the case, we would almost never need to have humans in the loop during the code review process because our automated tools could just do the peer review for us, unfortunately there is no perfect metric and not matter the code analysis we use, at least we need humans.
* **Radon** this is used for **cyclomatic complexity**.
* We want to able to reduce the number of decisions that a block of code contains.

    ```python
    pip install radon
    radon cc <file-name.py> -s ## cc -> cyclomatic complexity
    radon raw <file-name.py> ## raw metrics 
    

    # you can check more about this in documentation.
    ```

