# Robotic manipulations using SMT solver dReal

To generate Polygon Pushing Problem `.drh` files, customize the final few lines of `generate_pusher_problem.jl`, i.e.,
```julia
square = Polygon([[0,0], [1,0], [1,1], [0,1]])
triangle = Polygon([[0,0], [1,0], [0.5,sqrt(3)/2]])

open("dreach/push_test.drh", "w") do f
    write(f, pushing_problem(square, v_max=1.0, α_max = 10*pi/180, pin_θ=true))
end

open("dreach/push_test_tri.drh", "w") do f
    write(f, pushing_problem(triangle, v_max=1.0, α_max = 10*pi/180, pin_θ=true))
end
```
and run at the command line
```
> julia generate_pusher_problem.jl
```
The directory `dircol` contains the transcription of `push_test.drh` as a multi-phase direct collocation optimization problem. You will need [opty](https://github.com/csu-hmc/opty) and [IPOPT](https://projects.coin-or.org/Ipopt) in order to run it.
