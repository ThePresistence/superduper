# PRE-Commit Framework

We're going to integerate all the linting tools and code quality tools we have seen in the previous folders will bring them into the **workflow** that will make it as easy as possible to use them as the part of the continuous integeration.

There are lot of **branching strategy** we are gonna learn about the **trunk based development** means you can commit directly to the main or you create very short feature branches that get ideally merged back into the main in a same day.

Let's say you're following this branching strategy, now you created a PR. In this area you can have some **automated testing** for the linting, code-quality tests and more. Once the test is passed, you can approve the PR because testing is the necesarry part of the continuous integeration. If the test not passed,  you should review the code properly.

Let's see the approach of integerating all the tools we have seen before. There are lot of approach we can do this. Let's see one by one

## Approach 1

Writing a bash script and calling that in CI. It's the simplest and easiest way to do it. Please check the file `check-code-quality.sh`

## Approach 2

One problem you will face, if you're using the previous approach, the problem is all contributors should run this test before giving a PR. By making this operationalization takes lot of time. To make this automated **Git Hoooks** comes into the picture. **Git hooks** basically scripts that are run whenever you do certain git cli commands. Let's say we have a `git commit` command, after user writing this command, we can trigger a file that checks all the test, Isn't cool!




