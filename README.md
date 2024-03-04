# taking python to production learning 


### [All resources](https://ericriddoch.notion.site/Taking-Python-to-Production-A-Professional-Onboarding-Guide-799409731bf14c78a531ac779f1bd76d)


# 1. Install ohMyZsh on codespace ✨
```
rm -rf /home/codespace/.oh-my-zsh
sh -c "$(wget https://raw.githubusercontent.com/ohmyzsh/ohmyzsh/master/tools/install.sh -O -)"
```


### 1.1 Some more customization 
```
# To see the plugins 
cd ~ && cd .oh-my-zsh/plugins 

# Go to the .zshrc folder (step1: do this)
vim ~/.zshrc 
## change ZSH_THEME="bira"
```
### 1.2 Enable and disable plugins 
* Go to the .zshrc folder 
```
cd ~
vim .zshrc 
```
* Search all the plugins plugins folder and add those in the plugins 
```
cd ~
vim .zshrc 
plugins=(git web-search python pyenv virtualenv pip zsh-autosuggestions) 
```
