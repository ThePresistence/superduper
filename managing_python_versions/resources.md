### SEMANTIC VERSIONING 

**MAJOR.MINOR.PATCH** 

* **MAJOR** -> breaking change 
* **MINOR** -> Adding new functionality, feature, and backward compatible. 
* **PATCH** -> just a bug fix.

Official documentation for the semantic version: [click](https://semver.org/)


### [Managing Python environemnt](https://ericriddoch.notion.site/Resources-for-installing-pyenv-1cd45005cc6a4eef906f5165f259db94)


### Notes; 
* When-ever if you are specifying the pacakges in python, you should consider specifying by the package with **semantic versoning**. 
* This is because if you are specifying the package normally like `pip install numpy` sometimes they will upgrade the package that may contians some breakable change, that will affect your whole code. So instead you can specify like this `pip install numpy<3.0.2`, it will make sure it don't affect any of your exisiting package in your application or you can also mention like this `pip install numpy==3.2.0`. 
* The only drawback of this your package won't upgrade until you're upgrading manually. depends on the situation decide what you can use. 


#### Specify the version properly 

##### **0.0.0**

* If you're specifying the semantic version `0.0.0` client may think like your application in "alpha" or "unstable" phase. It's for getting feedback from the people. It can also communicate "this package can break anytime", so use with caution. 
* It's not generalized, it's published only for the feedback or showcasing. People don't consider your package is stable. 


##### **Deprecation warning** 

* Deprecation is whenever you remove an existing functionality or peice of existing API. 
* So whenver you're going to remove anything you should give the **Deprecation** warning. 
* 