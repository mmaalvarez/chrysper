# REQUIRED - git global setup
git config --global user.name "NAME SURNAME" # your gitlab name here

git config --global user.email "NAME.SURNAME@irbbarcelona.org" # your e-mail here


-----------------------------------------------------------------------------------


## update existing repository

### go to existing local repo directory
cd PATH/TO/REPO/

### add all files with changes (can specify them one by one, as well)
git add --all

### call commit, with info of which are the changes
git commit -m "I did this and that in there"

### pull possible remote repo updates into local
git pull -u origin master

### now push the commit from local to remote
git push -u origin master

### add commit list to changelog
git log --pretty="- %s" > CHANGELOG


------------------------------------------------------------------------------------


## clone existing repo into local
cd PATH/WHERE/CLONE_REPO/

git clone gitlab@fsupeksvr.irbbarcelona.pcb.ub.es:NAME_AUTHOR/NAME_PROJECT.git

### or clone specific repo branch
git clone --single-branch --branch NAME_BRANCH gitlab@fsupeksvr.irbbarcelona.pcb.ub.es:NAME_AUTHOR/NAME_PROJECT.git
