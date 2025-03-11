using Plots

function norm(p)
    return sqrt(sum(p .* p))
end

# Функция для вычисления градиента методом конечных разностей
function numerical_gradient(f, x, h=1e-5)
    grad = zeros(length(x))
    for i in eachindex(x)
        x_step = copy(x)
        x_step[i] += h
        grad[i] = (f(x_step) - f(x)) / h
    end
    return grad
end

# Метод золотого сечения для одномерной оптимизации
function golden_section_search(f, a, b, tol=1e-5)
    φ = (sqrt(5) - 1) / 2  
    c = b - (b - a) * φ
    d = a + (b - a) * φ
    while abs(b - a) > tol
        if f(c) < f(d)
            b = d
        else
            a = c
        end
        c = b - (b - a) * φ
        d = a + (b - a) * φ
    end
    return (a + b) / 2  
end

# Метод Милена и Контрелла
function milen_contrell_descent(f, x0; tol=1e-6, maxiter=1000)
    x = x0
    x_prev = x0
    history = [x]
    
    for k in 1:maxiter
        g = numerical_gradient(f, x)

        if norm(g) < tol
            println("Метод Милена-Контрелла сошелся на $k итерации")
            break
        end

        # Вычисляем разницу между шагами
        delta_x = x - x_prev

        # Оптимизируем λ0 (шаг по направлению градиента)
        λ0_func(λ) = f(x - λ * g)
        λ0_opt = golden_section_search(λ0_func, 1e-4, 0.5)

        # Оптимизируем λ1 (шаг по направлению предыдущего изменения)
        λ1_func(λ) = f(x + λ * delta_x)
        λ1_opt = golden_section_search(λ1_func, -0.5, 0.5)  

        # Обновляем x с учетом двух шагов
        x_new = x - λ0_opt * g + λ1_opt * delta_x

        if norm(x_new - x) < tol
            println("Метод сошелся по изменению координат на $k итерации")
            break
        end

        x_prev = x
        x = x_new
        push!(history, x)
    end

    return x, history
end

# Функции тестирования
# 1. Функция Равина (простая квадратичная)
function target_ravine(x)
    return sum(x .^ 2)
end

# 2. Функция Растригина (сложная, многоэкстремальная)
function target_rastrigin(x)
    A = 10
    return A * length(x) + sum(x_i^2 - A * cos(2π * x_i) for x_i in x)
end

# 3. Функция Швефеля (глобальный минимум далеко от нуля)
function target_schwefel(x)
    A = 418.9829
    return A * length(x) - sum(x_i * sin(sqrt(abs(x_i))) for x_i in x)
end

# 4. Функция Розенброка (использовалась ранее)
function target_rosenbrock(x)
    return (1 - x[1])^2 + 100 * (x[2] - x[1]^2)^2
end

# Визуализация траектории
function plot_descent(f, history; title="")
    n = 4
    x = range(-n, n, length=100)
    y = range(-n, n, length=100)
    Z = [f([xi, yi]) for xi in x, yi in y]

    plt = contour(x, y, Z, levels=30, title=title, aspect_ratio=:equal)
    traj_x = [p[1] for p in history]
    traj_y = [p[2] for p in history]

    plot!(plt, traj_x, traj_y, marker=:circle, linewidth=2, label=title)
    return plt
end

# Запуск метода и визуализация
function optimize_func_milen(f, title="")
    if !isempty(title)
        println("====== $title ======")
    end

    x0 = [4.0, 4.0]  

    (final_x, history) = milen_contrell_descent(f, x0)

    println("Milen-Conterell Optimization:   k = $(length(history))  x_min = $(final_x)  f(x_min) = ", f(final_x))
    println("==================================================\n")

    p1 = plot_descent(f, history, title="Milen-Conterell Optimization")
    display(p1)
end

# Запуск оптимизации на разных функциях
println("Минимизация функции Равина:")
optimize_func_milen(target_ravine, "Минимизация функции Равина")
readline()
# println("Минимизация функции Растригина:")
# optimize_func_milen(target_rastrigin, "Минимизация функции Растригина")

# println("Минимизация функции Швефеля:")
# optimize_func_milen(target_schwefel, "Минимизация функции Швефеля")

# println("Минимизация функции Розенброка:")
# optimize_func_milen(target_rosenbrock, "Минимизация функции Розенброка")
