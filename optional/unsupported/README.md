### Getting started
##### Clone:
###### HTTPS (github recommends)
```
git clone https://github.com/gismo/gismo.git
```
###### or SSH
```
git clone git@github.com:gismo/gismo.git
```
###### submodules
Submodules should be downloades by cmake automaticly. The are forced with `GISMO_<SUBMODULE>=ON` commands.
A example for this is `cmake . -DGISMO_UNSUPPORTED=ON`, as it is ussed in this submodule.

###### How to remove unused submodules (example)
```
git submodule deinit extensions/gsElasticity extensions/gsExaStencils extensions/motor
```
###### "Break up" detached state for submodule unsupported
```
cd gismo/extensions/unsupported  
git checkout master
```
###### or do it for all submodules - as example, don't do it
```
git submodule foreach git checkout master
```
##### Pull (Update to latest version)
```
git pull
git submodule update --remote
```
##### Build unsupported
```
cmake . -DGISMO_UNSUPPORTED=ON
make gismo_dev
```

This README contains only brief information. For more details, use [public wiki](http://gs.jku.at/gismo) and [internal wiki](https://github.com/gismo/internal/wiki).

[Information about Project Structure and Build Targets](https://github.com/gismo/internal/wiki/Git-Structure)

## Clone and use on your computer

This is a **submodule** of ***Gismo***. To use it, you have to clone Gismo with it's submodules.

### SSH Key

To do so, you have to use a SSH key - found in your [account settings](https://github.com/settings/keys).
_Please **read the manuals** there_. If you have still questions, ask us or your local system administrator.
> Note: There is also a way to set a **token** for a https connection - this is very long password, that you can than use with your account name.
> But, **we strongly recommend to use a SSH key!**
> You may also use a GPG key to sign your commits.

#### GPG Key
You can use a **GPG** key to sign your commit. Just follow the instructions above.

> Note: You have to use the same e-mail address as you use for commiting on github - therefore your "keep my email address private" address found under https://github.com/settings/emails

To ensure it works with your local git, add following lines to your global git config - type the commands into console with your email address and gpg key id.
```
git config --global user.email [your email]
git config --global user.name [your name]
git config --global user.signingkey [gpg key id]
git config --global commit.gpgsign true
```

You can look your settings with `git config --global -l`, there also exists `--system` and `--local` (only that project your currently in) instead of `--global`.

### Clone **gismo** with *submodules*:

    git clone --recurse-submodules git@github.com:gismo/gismo.git

The submodule folders are now checked out, but in a state called **detached** - this is a reference to a specific commit of each submodule.

#### unneeded submodules
You should *remove* submodules that you don't need.

First, look what submodules exist:

    git submodule
    
Remove unnedded submodules with:
    
    git submodule deinit <submodule>

To ***break up*** this state, go to the submodule you want to edit (for unsupported: `extensions/unsupported`) and call `git checkout master`

Now your submodule folder refers to the latest (called `HEAD`) commit.

### Pull (Update from remote)

To pull your ***superproject*** (this is a term git people use for the most outer project, for us, it is gismo) you can use `git pull`

To ***pull the subprojects***, call `git submodule update --remote`

### Push (Update your local changes to remote)
##### Commit local changes
First you have to ***save*** your changes in your local repository. This is also called **commit**:
`git commit`
##### Push to remote
Call `git push` to add your commits to remote server.

### Branch
Most times you want to do your changes in your own branch. To create one, go to the submodule and add a new branch with `git checkout -b [new_Branchname]`

When you ***push*** your branch for the first time, git wants you to add a upstream branch - at this moment, your branch is only local and not on github. Follow the output of git, it is somethink like `git push --set-upstream origin [new_Branchname]`.

After the first set of the upstream-branch, a normal `git push` will do it.

## Build unsupported

Create a build folder as usual - or use the one your IDE creates automaticly.

Add `GISMO_UNSUPPORTED=ON` to your CMake configuration. This is done by eighter `cmake . -DGISMO_UNSUPPORTED=ON` in your console or adding `-DGISMO_UNSUPPORTED=ON` to your IDEs settings.

