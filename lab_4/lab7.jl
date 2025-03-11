using Plots, LinearAlgebra

function norm(p)
    return sqrt(sum(p .* p))
end

# --- Определение целевых функций ---
function rosenbrock(x::Vector)
    return 100.0 * (x[2] - x[1]^2)^2 + (1 - x[1])^2
end

function ravine(x::Vector)
    return sum(x .^ 2)
end

function rastrigin(x::Vector)
    A = 10
    return A * length(x) + sum(x_i^2 - A * cos(2π * x_i) for x_i in x)
end

function schwefel(x::Vector)
    A = 418.9829
    return A * length(x) - sum(x_i * sin(sqrt(abs(x_i))) for x_i in x)
end

# --- Функция вычисления градиента ---
function numerical_gradient(f, x::Vector; h=1e-8)
    n = length(x)
    grad = zeros(n)
    for i in 1:n
        x_forward = copy(x)
        x_backward = copy(x)
        for i in 1:n
            x_forward[i] += h
            x_backward[i] -= h
            grad[i] = (f(x_forward) - f(x_backward)) / (2h)
        end
    end
    return grad
end

# --- Метод золотого сечения для одномерной оптимизации ---
function golden_search(f, a, b; tol=1e-5)
    phi = (sqrt(5) - 1) / 2
    c = b - phi * (b - a)
    d = a + phi * (b - a)
    while abs(b - a) > tol
        if f(c) < f(d)
            b = d
        else
            a = c
        end
        c = b - phi * (b - a)
        d = a + phi * (b - a)
    end
    return (a + b) / 2
end

# --- Оптимизация параметров λ0 и λ1 ---
function optimize_lambdas(g; init=[0.0, 0.0], tol=1e-5, max_iter=20)
    lambda = copy(init)
    for iter in 1:max_iter
        f1(lambda0) = g([lambda0, lambda[2]])
        lambda_new0 = golden_search(f1, 0.0, 1.0; tol=tol)
        f2(lambda1) = g([lambda_new0, lambda1])
        lambda_new1 = golden_search(f2, -1.0, 1.0; tol=tol)

        if norm([lambda_new0, lambda_new1] .- lambda) < tol
            lambda = [lambda_new0, lambda_new1]
            break
        end
        lambda = [lambda_new0, lambda_new1]
    end
    return lambda
end

# --- Основной метод многокритериального поиска ---
function multicriteria_search(f, x0::Vector; eps=1e-2, max_iter=10000)
    x = copy(x0)
    dx_prev = zeros(length(x))
    path = [copy(x)]
    
    for k in 1:max_iter
        grad = numerical_gradient(f, x)
        if norm(grad) < eps
            println("Критерий останова достигнут на итерации $k")
            break
        end

        g(lambdas) = f(x .- lambdas[1] .* grad .+ lambdas[2] .* dx_prev)
        
        lambda_opt = optimize_lambdas(g, init=[0.0, 0.0])
        lambda0, lambda1 = lambda_opt
        
        x_new = x .- lambda0 .* grad .+ lambda1 .* dx_prev

        dx_prev = x_new .- x
        x = copy(x_new)
        push!(path, copy(x))
        
        if mod(k, 2) == 0
            dx_prev .= 0.0
        end
    end
    return x, path
end

# --- Визуализация 2D траектории оптимизации ---
function plot_descent(f, history; title="")
    n = 10
    x = range(-n, n, length=100)
    y = range(-n, n, length=100)
    Z = [f([xi, yi]) for xi in x, yi in y]
    
    plt = contour(x, y, Z, levels=30, title=title, aspect_ratio=:equal)
    traj_x = [p[1] for p in history]
    traj_y = [p[2] for p in history]

    plot!(plt, traj_x, traj_y, marker=:circle, markersize=4, color=:red, linewidth=2, label="Trajectory")
    return plt
end

# --- Функция запуска оптимизации и вывода результатов ---
function optimize_func(f, title="")
    # x0 = [-1.2, 4.0]  # Начальная точка rosenbrock
    x0 = [-1.2, 1.0]  # Начальная точка schwefel. rastrigin
    # x0 = [4.0, 4.0]  # Начальная точка ravine
    (final_x, history) = multicriteria_search(f, x0)

    println("Multicriteria Optimization:   k = $(length(history))  x_min = $(final_x)  f(x_min) = ", f(final_x))
    println("==================================================\n")

    # Визуализация 2D траектории
    p1 = plot_descent(f, history, title="Multicriteria Optimization")
    display(p1)
end

# --- Запуск оптимизации на разных функциях ---
# println("Минимизация функции Равина:")
# optimize_func(ravine, "Минимизация функции Равина")

# println("Минимизация функции Растригина:")
# optimize_func(rastrigin, "Минимизация функции Растригина")

println("Минимизация функции Швефеля:")
optimize_func(schwefel, "Минимизация функции Швефеля")

# println("Минимизация функции Розенброка:")
# optimize_func(rosenbrock, "Минимизация функции Розенброка")

readline()
