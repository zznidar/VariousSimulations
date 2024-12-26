using Pkg
#Pkg.activate(".")
using Plots
using ShiftedArrays
using LinearAlgebra

# Make the plot interactive
plotlyjs()

# Plot settings

function setup(dimensionality=3)
    markersize = 0.5
    # Create 4 birds positioned as edges of a regular tetrahedron with side length a
    a = 1

    coordinates = Dict(
        2 => [1 sqrt(3) 0; 0 0 0; 2 0 0], # 2D with 3 birds (equilateral triangle)
        3 => [0 0 0; a 0 0; a/2 sqrt(3)/2*a 0; a/2 sqrt(3)/6*a sqrt(2/3)*a], # 3D with 4 birds (Regular tetrahedron)
        4 => [sqrt(3) sqrt(5) sqrt(10) sqrt(30); sqrt(3) sqrt(5) sqrt(10) -sqrt(30); sqrt(3) sqrt(5) -sqrt(40) 0; sqrt(3) -sqrt(45) 0 0; -4*sqrt(3) 0 0 0], # 4D with 5 birds (5-cell)
        5 => 4*[1/sqrt(15) 1/sqrt(10) 1/sqrt(6) 1/sqrt(3) 1; 1/sqrt(15) 1/sqrt(10) 1/sqrt(6) 1/sqrt(3) -1; 1/sqrt(15) 1/sqrt(10) 1/sqrt(6) -2/sqrt(3) 0; 1/sqrt(15) 1/sqrt(10) -sqrt(3)/sqrt(2) 0 0; 1/sqrt(15) -2*sqrt(2)/sqrt(5) 0 0 0; -sqrt(5)/sqrt(3) 0 0 0 0] # 5D with 6 birds (Regular hexateron)   
    )

    birds = coordinates[dimensionality]

    # Scatter the birds
    p = scatter(title="Birds simulation", titlefont = font(12,"Computer Modern"), birds[:, 1], birds[:, 2], birds[:, 3], markersize=markersize, legend = false, size=(1000,1000))

    return birds, p, markersize
end

function draw_birds(p, birds)
    ds = calculate_distances(birds)
    scatter!(p, title=ds, birds[:, 1], birds[:, 2], birds[:, 3], markersize=markersize)
end


function distance(p1, p2)
    return sqrt(sum((p1 - p2).^2))
end

function plot_distance(p, p1, p2)
    plot!(p, [p1[1], p2[1]], [p1[2], p2[2]], [p1[3], p2[3]], label = "Distance", line = :dash)
end

function calculate_distances(birds)
    n = size(birds, 1)
    distances = zeros(n, n)
    for i in 1:n
        for j in 1:n
            distances[i, j] = distance(birds[i, :], birds[j, :])
        end
    end
    return distances
end

function draw_distances(p, birds)
    # plot_distance for each pair
    for i in 1:size(birds, 1)
        for j in (i+1):size(birds, 1)
            plot_distance(p, birds[i, :], birds[j, :])
        end
    end
end

function update_positions(birds, ds=0.001)
    # move each bird towards the next one
    vectors = ShiftedArrays.circshift(birds, -1) - birds

    # Moving vectors (\vec{v}) should be normalised, i. e. always the same speed.
    normalize!(vectors)

    distance₀ = calculate_distances(birds)
    
    # Move the birds
    birds = birds + vectors*ds
    
    distance₁ = calculate_distances(birds)
    
    if(norm(distance₀) <= norm(distance₁))
        println("Birds have converged!")
        birds = birds - vectors*ds
        return birds, true
    end
    return birds, false
end


function simulate(dimensionality=3)
    stepSize = 0.01
    global birds, p, markersize = setup(dimensionality)
    draw_distances(p, birds)

    if dimensionality > 3
        dimFactor = 10
    else
        dimFactor = 1
    end

    global disthistory = []

    println(calculate_distances(birds)[1, :])

    for i in 1:5000
        birds, finished = update_positions(birds, stepSize)
        if(finished)
            println("Finished at iteration ", i)
            break
        end

        distances = calculate_distances(birds)
        push!(disthistory, distances[1, :])


        if(i % dimFactor == 0)
            markersize = 0.5 + (i*stepSize)/dimFactor
            draw_birds(p, birds)
        end
        #draw_birds(p, birds[1:1, :])
        if(i % (10*dimFactor) == 0)
            draw_distances(p, birds)
            println(calculate_distances(birds)[1, :])
        end
    end
    
    println(calculate_distances(birds)[1, :])
    display(p)
end


simulate(5) # Pass number of dimensions as param (2:5)

#Plots.html(p, "birds_simulation.html")

# Save the plot to a file
#savefig(p, "birds_simulation2.png")

plot(stack(disthistory)'[:, 2:end÷2+1])