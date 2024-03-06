# Write clean and efficient code :) 

## 1. Clean code ✨

* You should write the code in **idiomatic** way. Idiomatic means in software is whether or not your code conforms to standard style patterns for whatever the language or framework is that you're working on. Like proper variable name case for function, classes and constonts.
* In python, we have a special name of **idiomatic** code, that is called **pythonic** way. !!
* In software engineering interviews, some times people will check the way you written the code. So you should write the code in pythonic way !!
* Clean code means 
    - Easy to change, understand, safe/secure, relevant, and highly subjective.
    - Redability is the one of the key characteristics of clean code.
* You should use this principle called **refactoring** for better readability in code.  -> which is basically means incremental changes to make code more redable than it was when you first started. If you want to read a book for this you can read this book: <i>Refactoring improving the design of existing code</i>.
* If you're having stronger desire to write cleaner code, it will take years to achieve that. ONce you're able to write clean code, you will spend less time in debugging, better software, better at interviews and this is the organized way of designing.
* In few years, it will be possible python won't be there, people might use some other programming language but Principiles of clean code will never change. This is the crux of the programming language.

## 2. Python style guides

* Python has something called [**PEP8**](https://peps.python.org/pep-0008/) style guide. Basically python language is maintained by community and the community has a group called **steering council**. You can go to pep8 org for checking the python sytle guides and the proposals. 
* Again Python has some other styling guide by google called [**Google python style guide**](https://google.github.io/styleguide/pyguide.html). It's has become super popular. (most recommandable way!!)

## 3. Refractoring the code

* which is basically means incremental changes to make code more redable than it was when you first started.
* Refractoring is a skill it needs repeated practice to master. 
* Please read this [😌<-click](https://testdriven.io/blog/clean-code-python/) for how to write a clean python code. 

## 3 Auto Formatters 

* It will be automatically format the code what you have written.
* The first tool is [**BLACK**](https://github.com/psf/black)

### 3.1 Black

* Black is a python package it works in CLI.
* Black modifies the files in place until unless you're specifically telling to the black.
* This is one of the best tool to use in the teams for formatting a code for production use case. It contains lot of functionalities for formatting the code, best tool for pythonic way of formatting. 

### 3.2 VS code Format document

* You can type "ctrl+shift+p" and type "Format Document" to do the formatting. If you're doing first time you need to set up your formatter.
* You can install you're favourite formatter like black and type "ctrl+shift+p" and type "format document". now it will work. 

### 3.3 yapf

* Same as black formatter.

### 3.4 auto pep8 

* same as a black formatter. 

## 4.Rulers

* For single line, you should write only 79 characters according to the pep 8 standards.
* To see visually you can activate the ruler in vs code by go the setting and type for ruler and click the "edit in settings" and paste this dicionary inside the "editor.rulers": [ {"color": "lightblue", column: 79}]
* Sometimes it's ok to increase the 79 to 99 characters or sometimes you can go to 119 characters in single line.
* You might notice that all the numbers are ending at "9" this is because everything will be divisible by 4. 

---

## [Resources](https://github.com/phitoduck/python-software-development-course)

1. [clean python code](https://testdriven.io/blog/clean-code-python/)
