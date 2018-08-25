# Model Order Reduction (MOR)

Download pod2 folder from https://uofi.app.box.com/s/42r27y4kn3zkrfq82p3d3xin72f6p8oh and move the contents to pod folder to get started

## Cases

### Lid-Driven Cavity

A 2D case used to reproduce [1].

Domain: x,y \in [-1,1]
BCs: u(x,y) = \delta (y-1) * (1-x^2)^2

The top boundary condition is not the standard uniform velocity for consistency with [1]

[1]: Fick, Maday, Patera, Taddei; "A Reduced Basis Technique for Long-Time Unsteady Turbulent Flows"
