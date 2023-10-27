### TL;DR
##### Clone:
```
git clone --recurse-submodules git@github.com:gismo/gismo.git
git submodule foreach git checkout master
```
##### Pull (Update to latest version)
```
git pull
git submodule update --remote
```
##### Build motor
```
cmake . -DGISMO_UNSUPPORTED=ON -DGISMO_MOTOR=ON
make motor
```

This README contains only brief informations. For more details, use [public wiki](http://gs.jku.at/gismo) and [internal wiki](https://github.com/gismo/internal/wiki).
