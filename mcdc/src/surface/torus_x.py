"""
Torus: Implicit equation of a torus radially symmetric about the z-axis in the cartesian plane

r^2 = f(x, y, z) = ( sqrt( (x - A)^2 + (y - B)^2) - R )^2 + (z - C)^2

Where R is the radius of the shape as a whole, and r is the radius of the circle that is revolved to create the donut

As this definition will inevitably include square roots, an upper and lower surface of the torus will need to be combined to make the entire shape

Solving for z:
z = C (+ or -) sqrt( r^2 - ( sqrt( (x - A)^2 + (y - B)^2) - R )^2 )
"""



import math

from numba import njit

from mcdc.constant import (
    COINCIDENCE_TOLERANCE,
    INF,
)


@njit
def evaluate(particle_container, surface):
    """
    Description:
    Checking to see if the particle is on the surface in question

    Returns: (float)
    If the return is 0, the particle occupies the exact space of the shape
    """
  
    particle = particle_container[0]
    # Particle parameters
    x = particle["x"]
    y = particle["y"]
    z = particle["z"]

    # Surface parameters
    R = surface["R"]
    r = surface["r"]
    A = surface["A"]
    B = surface["B"]
    C = surface["C"]

    #Check if the particle is above or below the centerline of the torus and use the appropriate sign for the square root
    return (
        C + sqrt( r**2 - ( sqrt( (x - A)**2 + (y - B)**2) - R )**2 ) - z
        if z >= C else
        C - sqrt( r**2 - ( sqrt( (x - A)**2 + (y - B)**2) - R )**2 ) - z
    )


@njit
def reflect(particle_container, surface):
    particle = particle_container[0]
    # Particle coordinate
    x = particle["x"]
    y = particle["y"]
    z = particle["z"]
    ux = particle["ux"]
    uy = particle["uy"]
    uz = particle["uz"]

    # Surface coefficients
    R = surface["R"]
    r = surface["r"]
    A = surface["A"]
    B = surface["B"]
    C = surface["C"]

    # Surface normal
    dx = 
    dy = 
    dz = 
    norm = (dx**2 + dy**2 + dz**2) ** 0.5
    nx = dx / norm
    ny = dy / norm
    nz = dz / norm

    # Reflect
    c = 2.0 * (nx * ux + ny * uy + nz * uz)
    particle["ux"] -= c * nx
    particle["uy"] -= c * ny
    particle["uz"] -= c * nz


@njit
def get_normal_component(particle_container, surface):
    particle = particle_container[0]
    # Particle coordinate
    x = particle["x"]
    y = particle["y"]
    z = particle["z"]
    ux = particle["ux"]
    uy = particle["uy"]
    uz = particle["uz"]

    # Surface coefficients
    R = surface["R"]
    r = surface["r"]
    A = surface["A"]
    B = surface["B"]
    C = surface["C"]

    # Surface normal
    dx = 
    dy = 
    dz = 
    norm = (dx**2 + dy**2 + dz**2) ** 0.5
    nx = dx / norm
    ny = dy / norm
    nz = dz / norm

    return nx * ux + ny * uy + nz * uz


@njit
def get_distance(particle_container, surface):
    particle = particle_container[0]
    # Particle coordinate
    x = particle["x"]
    y = particle["y"]
    z = particle["z"]
    ux = particle["ux"]
    uy = particle["uy"]
    uz = particle["uz"]

    # Surface coefficients
    A = surface["A"]
    B = surface["B"]
    C = surface["C"]
    D = surface["D"]
    E = surface["E"]
    F = surface["F"]
    G = surface["G"]
    H = surface["H"]
    I = surface["I"]

    # Coincident?
    f = evaluate(particle_container, surface)
    coincident = abs(f) < COINCIDENCE_TOLERANCE
    if coincident:
        # Moving away or tangent?
        if (
            get_normal_component(particle_container, surface)
            >= 0.0 - COINCIDENCE_TOLERANCE
        ):
            return INF

    # Quadratic equation constants
    a = (
        A * ux * ux
        + B * uy * uy
        + C * uz * uz
        + D * ux * uy
        + E * ux * uz
        + F * uy * uz
    )
    b = (
        2 * (A * x * ux + B * y * uy + C * z * uz)
        + D * (x * uy + y * ux)
        + E * (x * uz + z * ux)
        + F * (y * uz + z * uy)
        + G * ux
        + H * uy
        + I * uz
    )
    c = f

    determinant = b * b - 4.0 * a * c

    # Roots are complex : no intersection
    # Roots are identical: tangent
    # ==> return huge number
    if determinant <= 0.0:
        return INF
    else:
        # Get the roots
        denom = 2.0 * a
        sqrt = math.sqrt(determinant)
        root_1 = (-b + sqrt) / denom
        root_2 = (-b - sqrt) / denom

        # Coincident?
        if coincident:
            return max(root_1, root_2)

        # Negative roots, moving away from the surface
        if root_1 < 0.0:
            root_1 = INF
        if root_2 < 0.0:
            root_2 = INF

        # Return the smaller root
        return min(root_1, root_2)
