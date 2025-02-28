using Plots
using Printf

# Аппроксимация первой производной
function approximate_first_derivative(f, x, h=1e-6)
    return (f(x + h) - f(x - h)) / (2h)
end

# Аппроксимация второй производной
function approximate_second_derivative(f, x, h=1e-6)
    return (f(x + h) - 2f(x) + f(x - h)) / h^2
end

# Функция для нахождения интервала унимодальности
function find_unimodal_interval(f, a, b, h=1e-5)
    # Проверим первую производную на интервале [a, b]
    x_vals = a:h:b
    first_derivative = [approximate_first_derivative(f, xi, h) for xi in x_vals]
    
    # Найдем точку, где производная изменяет знак
    sign_changes = 0
    start_index = 1
    end_index = length(x_vals)
    
    # Ищем первую точку, где меняется знак производной
    for i in 2:length(first_derivative)
        if first_derivative[i] * first_derivative[i-1] < 0
            sign_changes += 1
            end_index = i
            break
        end
    end
    
    # Если не было изменений знака на интервале, то это унимодальная функция
    if sign_changes == 0
        return a, b
    else
        # Возвращаем интервал, на котором производная меняет знак
        return x_vals[start_index], x_vals[end_index]
    end
end

# Функция для нахождения интервала унимодальности по второй производной
function find_unimodal_interval_second_derivative(f, a, b, h=1e-5)
    # Проверим вторую производную на интервале [a, b]
    x_vals = a:h:b
    second_derivative = [approximate_second_derivative(f, xi, h) for xi in x_vals]
    
    # Найдем точку, где вторая производная изменяет знак
    sign_changes = 0
    start_index = 1
    end_index = length(x_vals)
    
    # Ищем первую точку, где меняется знак второй производной
    for i in 2:length(second_derivative)
        if second_derivative[i] * second_derivative[i-1] < 0
            sign_changes += 1
            end_index = i
            break
        end
    end
    
    # Если не было изменений знака на интервале, то это унимодальная функция
    if sign_changes == 0
        return a, b
    else
        # Возвращаем интервал, на котором производная меняет знак
        return x_vals[start_index], x_vals[end_index]
    end
end

# Пример функции
f(x) = 1/3*(x^3) + 2*x^2 + 3*x 
a, b = -5, 5

interval = find_unimodal_interval(f, a, b)

println("Интервал, на котором функция является унимодальной: ", interval)

interval2 = find_unimodal_interval_second_derivative(f, a, b)

println("Интервал, на котором функция является унимодальной по второй производной: ", interval2)
# Метод бисекции
function bisection(f, a, b, eps)
    a = Float64(a)
    b = Float64(b)
    intervals = [(a, b)]
    iters = 0
    x_vals, y_vals = [a], [f(a)]  # Для графика
    while b - a > eps
        iters += 1
        m = (a + b) / 2
        if f(m - eps) < f(m + eps)
            b = m
        else
            a = m
        end
        push!(intervals, (a, b))
        push!(x_vals, m)  # Сохраняем x для графика
        push!(y_vals, f(m))  # Сохраняем y для графика
    end
    min = (a + b) / 2
    return min, f(min), iters, x_vals, y_vals
end

# Метод золотого сечения
function golden_section(f, a, b, eps)
    k = (sqrt(5) - 1) / 2
    x1 = a + (1 - k) * (b - a)
    x2 = a + k * (b - a)
    a = Float64(a)
    b = Float64(b)
    intervals = [(a, b)]
    iters = 0
    x_vals, y_vals = [a], [f(a)]  # Для графика
    while abs(x1 - x2) > eps
        iters += 1
        if f(x1) <= f(x2)
            b = x2
            x2 = x1
            x1 = a + b - x1
        else
            a = x1
            x1 = x2
            x2 = a + b - x2
        end
        push!(intervals, (x1, x2))
        push!(x_vals, (x1 + x2) / 2)  # Сохраняем x для графика
        push!(y_vals, f((x1 + x2) / 2))  # Сохраняем y для графика
    end
    min = (a + b) / 2
    return min, f(min), iters, x_vals, y_vals
end

# Метод Фибоначчи
function fibonacci(n)
    if n == 1 || n == 2
        return 1
    end
    return fibonacci(n - 1) + fibonacci(n - 2)
end

# Метод Фибоначчи
function fibonacci_search(f, a, b, eps)
    n = 10
    fib = [fibonacci(i) for i in 1:n]  # Генерация чисел Фибоначчи
    x1 = a + (fib[n-2] / fib[n]) * (b - a)
    x2 = a + (fib[n-1] / fib[n]) * (b - a)
    a = Float64(a)
    b = Float64(b)
    intervals = [(a, b)]
    iters = 0
    x_vals, y_vals = [a], [f(a)]  # Для графика
    for k in 1:(n-3)
        iters += 1
        if f(x1) > f(x2)
            a = x1
            x1 = x2
            x2 = a + (fib[n-k-1] / fib[n-k]) * (b - a)
        else
            b = x2
            x2 = x1
            x1 = a + (fib[n-k-2] / fib[n-k]) * (b - a)
        end
        push!(x_vals, (a + b) / 2)  # Сохраняем x для графика
        push!(y_vals, f((a + b) / 2))  # Сохраняем y для графика
    end
    min = (a + b) / 2
    return min, f(min), iters, x_vals, y_vals
end

# Поиск минимума методом бисекции с графиком
x_min_bisection, f_min_bisection, iters_bisection, x_vals_bisection, y_vals_bisection = bisection(f, a, b, 1e-5)
println("Метод бисекции:")
println("Приближенное значение x*: ", x_min_bisection)
println("Количество итераций: ", iters_bisection)
println(x_vals_bisection)
# Поиск минимума методом золотого сечения с графиком
x_min_golden, f_min_golden, iters_golden, x_vals_golden, y_vals_golden = golden_section(f, a, b, 1e-5)
println("Метод золотого сечения:")
println("Приближенное значение x*: ", x_min_golden)
println("Количество итераций: ", iters_golden)
println(x_vals_golden)
# Поиск минимума методом Фибоначчи с графиком
x_min_fibonacci, f_min_fibonacci, iters_fibonacci, x_vals_fibonacci, y_vals_fibonacci = fibonacci_search(f, a, b, 1e-5)
println("Метод Фибоначчи:")
println("Приближенное значение x*: ", x_min_fibonacci)
println("Количество итераций: ", iters_fibonacci)
println(x_vals_fibonacci)

# Построение графиков для каждого метода в разных подграфиках
p1 = plot(a:0.1:b, f.(a:0.1:b), label="f(x) = 1/3*(x^3) + 2*x^2 + 3*x ", xlabel="x", ylabel="f(x)", linewidth=2, color=:blue)
scatter!(p1, x_vals_bisection, y_vals_bisection, color=:red, label="Точки итерации (бисекция)", markersize=3)
for i in 1:length(x_vals_bisection)
    annotate!(p1, (x_vals_bisection[i], y_vals_bisection[i]), text("$i", 15, :center))
end

p2 = plot(a:0.1:b, f.(a:0.1:b), label="f(x) = 1/3*(x^3) + 2*x^2 + 3*x ", xlabel="x", ylabel="f(x)", linewidth=2, color=:blue)
scatter!(p2, x_vals_golden, y_vals_golden, color=:green, label="Точки итерации (золотое сечение)", markersize=3)
for i in 1:length(x_vals_golden)
    annotate!(p2, (x_vals_golden[i], y_vals_golden[i]), text("$i", 8, :center))
end

p3 = plot(a:0.1:b, f.(a:0.1:b), label="f(x) = 1/3*(x^3) + 2*x^2 + 3*x ", xlabel="x", ylabel="f(x)", linewidth=2, color=:blue)
scatter!(p3, x_vals_fibonacci, y_vals_fibonacci, color=:orange, label="Точки итерации (Фибоначчи)", markersize=3)
for i in 1:length(x_vals_fibonacci)
    annotate!(p3, (x_vals_fibonacci[i], y_vals_fibonacci[i]), text("$i", 8, :center))
end

# Отображение всех 3 графиков
plot(p1, p2, p3, layout=(3,1), size=(800, 600))


gui()
readline()

