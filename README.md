# Вариант 4: Движение твердого тела в жидкости под действием силы архимеда, силы тяжести и сопротивления среды
Численное интегрирование выполняется с использованием **метода Рунге-Кутты 4-го порядка**
### Доп задание: 
реализовать симуляцию при помощи метода Эйлера и метода средней точки 
## Запуск  
https://github.com/BelyaevMaks/computational-mathematics.git
g++ euler.cpp -o euler_simulation -lGL -lGLU -lglut  
./euler_simulation
g++ midpoint.cpp -o midpoint_simulation -lGL -lGLU -lglut  
./midpoint_simulation
g++ runge_kutta_4.cpp -o runge_kutta_4_simulation -lGL -lGLU -lglut  
./runge_kutta_4
