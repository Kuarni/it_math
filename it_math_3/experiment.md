# Эксперимент

---
Были взяты следующие алгоритмы:

- [numpy.linalg.svd](https://numpy.org/doc/stable/reference/generated/numpy.linalg.svd.html) - `npSVD`
- [power method SVD](https://www.jeremykun.com/2016/05/16/singular-value-decomposition-part-2-theorem-proof-algorithm/) - `pwmSVD`
- [block power method SVD](https://www.emis.de/journals/ASUO/mathematics/anale2015vol2/Bentbib_A.H.__Kanber_A..pdf) - `bpmSVD`

---

## Время выполнения

Реализация из numpy является самой быстрой, в то время как блочный метод - самый медленный.

## Сравнение результатов

### параметры:

- изображение сжимается в 2 раза
- для power и block метода используется `tolerance = 10e-2`

### черно-белое изображение

| исходное                                              | numpy                                                                  | power                                                                   | block power                                                             |
|-------------------------------------------------------|------------------------------------------------------------------------|-------------------------------------------------------------------------|-------------------------------------------------------------------------|
| ![black-white.bmp](img_src/the_end_24/the_end_24.bmp) | ![black-white.bmp](img_src/the_end_24/the_end_24_npSVD_decompress.bmp) | ![black-white.bmp](img_src/the_end_24/the_end_24_pwmSVD_decompress.bmp) | ![black-white.bmp](img_src/the_end_24/the_end_24_bpmSVD_decompress.bmp) 

тут явно проигрывает power метод, из-за того, что он добавил цветные шумы, хотя изначально картинка была черно-белой.
Numpy
и block power показали примерно одинаковый результат. Оба алгоритма добавили примерно одинаковое количество шумов.

### простое цветное изображение

| исходное                                   | numpy                                                       | power                                                        | block power                                                  |
|--------------------------------------------|-------------------------------------------------------------|--------------------------------------------------------------|--------------------------------------------------------------|
| ![simple-color.bmp](img_src/4k24/4k24.bmp) | ![simple-color.bmp](img_src/4k24/4k24_npSVD_decompress.bmp) | ![simple-color.bmp](img_src/4k24/4k24_pwmSVD_decompress.bmp) | ![simple-color.bmp](img_src/4k24/4k24_bpmSVD_decompress.bmp) |

power метод снова проигрывает - у него получилось более шумное изображение. У numpy и block power примерно одинаковые
результаты.

### изображение в низком разрешении с большим количеством деталей

| исходное                                            | numpy                                                                | power                                                                 | block power                                                           |
|-----------------------------------------------------|----------------------------------------------------------------------|-----------------------------------------------------------------------|-----------------------------------------------------------------------|
| ![nagibator.bmp](img_src/Kirilenich/Kirilenich.bmp) | ![nagibator.bmp](img_src/Kirilenich/Kirilenich_npSVD_decompress.bmp) | ![nagibator.bmp](img_src/Kirilenich/Kirilenich_pwmSVD_decompress.bmp) | ![nagibator.bmp](img_src/Kirilenich/Kirilenich_bpmSVD_decompress.bmp) |

тут результаты одинаковы

### большое изображение с большим кол-вом деталей

| исходное                              | numpy                                                  | power                                                   | block power                                             |
|---------------------------------------|--------------------------------------------------------|---------------------------------------------------------|---------------------------------------------------------|
| ![big.bmp](img_src/Hazbik/Hazbik.bmp) | ![big.bmp](img_src/Hazbik/Hazbik_npSVD_decompress.bmp) | ![big.bmp](img_src/Hazbik/Hazbik_pwmSVD_decompress.bmp) | ![big.bmp](img_src/Hazbik/Hazbik_bpmSVD_decompress.bmp) |

На данном изображении разница между алгоритмами невидна.

Попробуем сжать не в 2, а в 4 раза

| исходное                              | numpy                                                    | power                                                     | block power                                               |
|---------------------------------------|----------------------------------------------------------|-----------------------------------------------------------|-----------------------------------------------------------|
| ![big.bmp](img_src/Hazbik/Hazbik.bmp) | ![big.bmp](img_src/Hazbik/Hazbik_npSVD_4_decompress.bmp) | ![big.bmp](img_src/Hazbik/Hazbik_pwmSVD_4_decompress.bmp) | ![big.bmp](img_src/Hazbik/Hazbik_bpmSVD_4_decompress.bmp) |

Тут тоже не видно.

### простые фигуры

| исходное                                       | numpy                                                           | power                                                            | block power                                                      |
|------------------------------------------------|-----------------------------------------------------------------|------------------------------------------------------------------|------------------------------------------------------------------|
| ![simp_fig.bmp](img_src/simp_fig/simp_fig.bmp) | ![simp_fig.bmp](img_src/simp_fig/simp_fig_npSVD_decompress.bmp) | ![simp_fig.bmp](img_src/simp_fig/simp_fig_pwmSVD_decompress.bmp) | ![simp_fig.bmp](img_src/simp_fig/simp_fig_bpmSVD_decompress.bmp) |

Тут видно, что на power методе больше шумов. Результаты numpy и block power одинаковы.

Посмотрим результат на сжатии в 10 раз:

| исходное                                       | numpy                                                             | power                                                              | block power                                                        |
|------------------------------------------------|-------------------------------------------------------------------|--------------------------------------------------------------------|--------------------------------------------------------------------|
| ![simp_fig.bmp](img_src/simp_fig/simp_fig.bmp) | ![simp_fig.bmp](img_src/simp_fig/simp_fig_npSVD_4_decompress.bmp) | ![simp_fig.bmp](img_src/simp_fig/simp_fig_pwmSVD_4_decompress.bmp) | ![simp_fig.bmp](img_src/simp_fig/simp_fig_bpmSVD_4_decompress.bmp) |

Тут разницы между алгоритмами особо нет. Но можно сделать вывод, что SVD алгоритмы лучше всего работают с прямоугольными фигурами.

---

# Вывод

Алгоритмы можно ранжировать так:

1. numpy (из-за того, что он сильно быстрее, чем block power)
2. block power method
3. power method (проигрывает по качеству)
