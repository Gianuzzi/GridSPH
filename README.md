# GridSPH
Repository with code to grid SPH particles. 
Uses  [Röttgers & Alexander, 2018](https://arxiv.org/abs/1803.03652): [https://arxiv.org/pdf/1803.03652.pdf](https://arxiv.org/pdf/1803.03652.pdf)

## Requerimientos
Compilar y ejecutar:
- ```gfortran```

Recompilar y ejecutar (en caso de realizar modificaciones al código que cambien las dependencias):
- ```pip``` 
- ```python```
  - [fortdepend](https://github.com/ZedThree/fort_depend.py) (paquete de python utilzado para generar nuevo archivo de dependencias)

## Modo de uso
Si se usa el código tal y como está, solo es necesario compilarlo y ejecutarlo.
En caso de realizarle modificaciones que alteren las dependencias, es necesario instalar [fortdepend](https://github.com/ZedThree/fort_depend.py).
### Crear y ejecutar
```console
$ make
$ ./main
```
En caso de que no exista el archivo [main.dep](./main.dep), es necesario tener instalado [fortdepend](https://github.com/ZedThree/fort_depend.py) para generarlo.
### Instalar _fortdepend_ y generar nuevo archivo de dependencias
```console
$ make install
$ make deps
```
### Borrar archivos de compilación para luego recompilar
```console
$ make clean
```
En este caso, también se borra el archivo de dependencias [main.dep](./main.dep), por lo que luego será necesario regenerarlo.

## Licencia
Distribucído bajo la Licencia MIT.  (Archivo [LICENSE](./LICENSE))

## Autor
Emmanuel Gianuzzi - egianuzzi@mi.unc.edu.ar