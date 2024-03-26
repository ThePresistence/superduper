# pytest

Whenever you're writing a test case, create a folder called "test". It's a recommendable way to do that.

## Testing theory

The test should not tide the implementation of the function, they should only test the behaviour of the function and not just function too. This kind of testing is called **black box testing**. The idea of the black box test is you cannot see the internals.

The **opposite** of the **black box** is **glass box tests**, it does exposes the internals of the function. We want our **test** can should be **black box tests** the reason should be we just need to test the input and output of the function not the internals.

By the way, the word people use the describe the tests that depends on implementation details is **brittle**, developers don't want the test case to be **brittle**.

**Test case universe**, it means the set of all possible test cases or all possible inputs that we could run through this function. Adding more test cases will prove that the written code is super robust. The trade of off increasing the number of test in **test case universe** is that the **more the test you have the slower the code**.

## Parameterized test cases

Often times you can divide the test case universe into different categories or sections. Let's say you have a function to add, instead of having the test case like infinte numbers you can divide the test cases like this positive numbers, negative numbers and zero. IF the function is working for this numbers, it means it's working good. Here we're are dividing the test case into some sections. In simple words we are adding parameters to the test cases.

## Markers

This is the another way of testing the code. Makers are basically a way for us to assign lables or tags so that we can place our various tests into different categories. Basically we can tag the function like **fast** and **slow**. If you want to run the **fast** functions first and **slow** function second, we can use the markers to do that.

There is an another term called **unit test**, it's a type of test that tests on very specific thing or set differently it tests a single unit of your code base that could be a single unit of your app's functionality.

## Test coverage

Test coverage is the percentage of lines of your program that actually got executed. So If you have full est coverage, it means that every single line of code in your program was executed by atleast one test. To do that we have a tool called **coverage.py**. This is a **separate** from the pytest. If you want to integerate with **pytest** we can include the **plugin** in pytest called **pytest-cov**. To add pluging simply install the plugin.

```code
pip install pytest-cov
```

And to see the test coverage, you can include this parameters in the pytest cli.

```code
pytest --cov --cov-report html # that's all very simple
```

After that it will save the test coverage results in htmlcov folder. To host the application, run the command bellow.

```python
python -m http.server -d ./htmlcov
```

























