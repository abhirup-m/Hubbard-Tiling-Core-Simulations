# define the edges of the Brillouin zone
(@isdefined K_MIN) || const K_MIN = -pi
(@isdefined K_MAX) || const K_MAX = pi

# tolerance for identifying energy contours
(@isdefined TOLERANCE) || const TOLERANCE = 1e-8

# tolerance for identifying RG fixed points
(@isdefined RG_TOL) || const RG_TOL = 1e-4

# overlap integral 't' set to 1
(@isdefined HOP_T) || const HOP_T = 1.0

(@isdefined OMEGA_BY_t) || const OMEGA_BY_t = -2 * HOP_T

(@isdefined NODAL_POINTS) || const NODAL_POINTS = [(-π/2, -π/2), (-π/2, π/2), (π/2, π/2), (π/2, -π/2)]
(@isdefined ANTINODAL_POINTS) || const ANTINODAL_POINTS = [(-π, 0.), (0., π), (π, 0.), (0., -π)]
