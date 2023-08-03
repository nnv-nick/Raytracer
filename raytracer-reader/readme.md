# Raytracer P2

Это вторая подзадача [домашнего задания про рейтрейсер](../raytracer). В рамках этой задачи вам нужно реализовать ридер .obj файлов, в которых задается сцена для рейтрейсера.

Описание формата можно посмотреть на вики: https://en.wikipedia.org/wiki/Wavefront_.obj_file. В нашем случае можно отметить следующее:

* Необходимо поддерживать только строчки с v, vt, vn, f, mtllib и usemtl, остальное нужно игнорировать.
* Нумерация сущностей в f глобальная, т.е. нужно игнорировать различные группирующие модификаторы вроде g или o.
* В f может быть задано произвольное число вершин, но вы можете сразу нарезать такой многоугольник на треугольники $`(P_0, P_i, P_{i+1})`$, где $`i`$ пробегает по вершинам. Необходимо поддерживать все возможные варианты
задания вершин в f (с индексом нормали и/или текстуры). При этом саму поддержку текстурирования делать не нужно.
* Гарантируется, что файл, определенный в mtllib, находится в той же директории, что и .obj файл.
* Для удобства в нашем задании также нужно уметь обрабатывать строчки вида `S x y z r`. Такая строка задает сферу с центром в `(x, y, z)`
радиуса `r`.
* Помимо этого нужно обрабатывать строки вида `P x y z r g b`. Такая строка задает точечный источник света с координатами `(x, y, z)`, имеющий
интенсивность `(r, g, b)`. Обратите внимание, что `(r, g, b)` необязательно лежат в диапазоне `[0,1].`

Касательно .mtl файлов описание приведено там же на вики, есть следующие нюансы:

* Необходимо поддерживать строки с newmtl, Ka (поле `ambient_color` в `Material`), Kd (`diffuse_color`), Ks (`specular_color`), Ke (`intensity`), Ns (`specular_exponent`), Ni (`refraction_index`). Все остальное нужно игнорировать.
* Модификатор Ke нестандартный, но часто встречается в реальных .mtl файлах. В нашем задании его нужно обрабатывать так же, как и Ka, т.е.
как аддитивную прибавку к суммарному освещению объекта.
* В нашем задании нужно поддерживать следующий дополнительный модификатор `al a b c` (`albedo`), задающий светоотражающие характеристики материала. Как они учитываются в вычислении освещения описано в основном задании. Если модификатор не задан, то по умолчанию значения должны быть `al 1 0 0`.

Все необходимые для реализации функции приведены в `scene.h`. Функции-геттеры должны возвращать объекты в том порядке, в котором они заданы в файле. В Read и GetMaterials ключом мапа является имя материала.