Birds simulation for Problem 2 of Physics Cup: https://physicscup.ee/physics-cup-taltech-2025-problem-2/

## Usage
Once you have Julia and needed packages installed, run the `j.jl` file. Pass the number of dimensions as the argument of the `simulate(dimensionality)` function.

Supported dimensionalities:
* 2: [equilateral triangle](https://en.wikipedia.org/wiki/Equilateral_triangle)
* 3: [Regular tetrahedron](https://en.wikipedia.org/wiki/Tetrahedron#Cartesian_coordinates)
* 4: [5-cell](https://en.wikipedia.org/wiki/5-cell#Coordinates)
* 5: [Regular hexateron](https://en.wikipedia.org/wiki/5-simplex#Regular_hexateron_cartesian_coordinates)

Dimensionalities higher than 3D only plot 3 dimensions (4th, 5th coordinates of each point are disregarded). Still, coordinates and distances are calculated in full dimensions.