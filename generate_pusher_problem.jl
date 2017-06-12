import QuadGK.quadgk
cross2d(x, y) = x[1]*y[2] - x[2]*y[1]

immutable Polygon
    points::Vector{Vector{Float64}}
    edges::Vector{Vector{Float64}}
    c::Float64

    function Polygon(points)
        n = length(points)
        A = sum(cross2d(points[i], points[mod1(i+1,n)]) for i in 1:n)/2
        cm = sum((points[i]+points[mod1(i+1,n)])*cross2d(points[i], points[mod1(i+1,n)]) for i in 1:n)/(6*A)
        points = [p - cm for p in points]
        edges = diff(points)
        push!(edges, points[1] - points[end])
        mmax = 0.
        for i in 1:n
            r1 = hypot(points[i]...)
            r2 = hypot(points[mod1(i+1,n)]...)
            θ = mod2pi(atan2(reverse(points[mod1(i+1,n)])...) - atan2(reverse(points[i])...))
            α = (r2*cos(θ) - r1)/(r2*sin(θ))
            mmax += quadgk(t -> (r1/(cos(t) - α*sin(t)))^3,0,θ)[1]/3
        end
        c = mmax/A
        new(points, edges, c)
    end
end

function pushing_problem(P::Polygon; push_points = P.points + (P.edges ./ 2),
                                     t_max = 10,
                                     p_bound = 5,
                                     pin_θ = false,
                                     v_max = 10,
                                     α_max = 30*pi/180,
                                     x_min = -p_bound, x_max = p_bound,
                                     y_min = -p_bound, y_max = p_bound,
                                     vpx_min = -v_max, vpx_max = v_max,
                                     vpy_min = -v_max, vpy_max = v_max)
    n = length(push_points)
    body_coords = join((join(("#define rox$(i)_b ($(xy[1]))",
                              "#define roy$(i)_b ($(xy[2]))",
                              "#define r$(i)sq ($(norm(xy)^2))"), "\n")
                        for (i, xy) in enumerate(push_points)), "\n")
    world_coords = join((join(("#define rox$(i)_w (cos(theta)*rox$(i)_b - sin(theta)*roy$(i)_b)",
                               "#define roy$(i)_w (sin(theta)*rox$(i)_b + cos(theta)*roy$(i)_b)"), "\n")
                         for i in 1:n), "\n")
    product_terms = join((join(("#define rox$(i)roy$(i) (rox$(i)_b*roy$(i)_b*cos(2*theta) + 0.5*(rox$(i)_b^2 - roy$(i)_b^2)*sin(2*theta))",
                                "#define rox$(i)sq (0.5*(r$(i)sq + (rox$(i)_b^2 - roy$(i)_b^2)*cos(2*theta) - 2*rox$(i)_b*roy$(i)_b*sin(2*theta)))",
                                "#define roy$(i)sq (0.5*(r$(i)sq + (roy$(i)_b^2 - rox$(i)_b^2)*cos(2*theta) + 2*rox$(i)_b*roy$(i)_b*sin(2*theta)))"), "\n")
                          for i in 1:n), "\n")
    jump_to_any = join(("        (true) ==> @$(i+1) (and (x' = x) (y' = y) (theta' = theta) (vpx' <= $(vpx_max)) (vpy' <= $(vpy_max)));" for i in 1:n), "\n")

    function pushing_mode(i)
"""
// pushing on point $(i) = $(push_points[i])
{ mode $(i+1);
    invt:
        (vpx*$(P.edges[i][1]) + vpy*$(P.edges[i][2]))^2 <= $(sin(α_max)^2)*(vpx^2 + vpy^2)*$(dot(P.edges[i], P.edges[i]));
        -vpx*$(P.edges[i][2]) + vpy*$(P.edges[i][1]) > 0;
    flow:
        d/dt[x] = ((c^2 + rox$(i)sq)*vpx + rox$(i)roy$(i)*vpy)/(c^2 + r$(i)sq);
        d/dt[y] = (rox$(i)roy$(i)*vpx + (c^2 + roy$(i)sq)*vpy)/(c^2 + r$(i)sq);
        d/dt[theta] = $(pin_θ ? 0 : 1)*(rox$(i)_w*vpy - roy$(i)_w*vpx)/(c^2 + r$(i)sq);
        d/dt[vpx] = 0;
        d/dt[vpy] = 0;
    jump:
        (true) ==> @1 (and (x' = x) (y' = y) (theta' = theta) (vpx' <= $(vpx_max)) (vpy' <= $(vpy_max)));
        //(true) ==> @$(i+1) (and (x' = x) (y' = y) (theta' = theta) (vpx' <= $(vpx_max)) (vpy' <= $(vpy_max)));
}
"""
    end

    pushing_modes = join((pushing_mode(i) for i in 1:n), "\n")

"""
#define c ($(P.c))
$body_coords
$world_coords
$product_terms

[0, $(t_max)] time;
// center of mass (x, y) bounds
[$(x_min), $(x_max)] x;
[$(y_min), $(y_max)] y;
[-20, 20] theta;
// quasistatic pushing velocity bounds
[$(vpx_min), $(vpx_max)] vpx;
[$(vpy_min), $(vpy_max)] vpy;

// stationary
{ mode 1;
    invt:
    flow:
        d/dt[x] = 0;
        d/dt[y] = 0;
        d/dt[theta] = 0;
        d/dt[vpx] = 0;
        d/dt[vpy] = 0;
    jump:
$jump_to_any
}

$pushing_modes

init:
@1  (and (x = 0) (y = 0) (theta = 0));
goal:
@1  (and (x > 3) (y > 1));
"""

end


# here be numerical tricks/hacks/dragons
function pushing_problem_stct_trick(P::Polygon; push_points = P.points + (P.edges ./ 2),
                                     t_max = 10,
                                     p_bound = 5,
                                     v_max = 10,
                                     α_max = 30*pi/180,
                                     x_min = -p_bound, x_max = p_bound,
                                     y_min = -p_bound, y_max = p_bound,
                                     vpx_min = -v_max, vpx_max = v_max,
                                     vpy_min = -v_max, vpy_max = v_max)
    n = length(push_points)
    body_coords = join((join(("#define rox$(i)_b ($(xy[1]))",
                              "#define roy$(i)_b ($(xy[2]))",
                              "#define r$(i)sq ($(norm(xy)^2))"), "\n")
                        for (i, xy) in enumerate(push_points)), "\n")
    world_coords = join((join(("#define rox$(i)_w (ct*rox$(i)_b - st*roy$(i)_b)",
                               "#define roy$(i)_w (st*rox$(i)_b + ct*roy$(i)_b)"), "\n")
                         for i in 1:n), "\n")
    jump_to_any = join(("        (true) ==> @$(i+1) (and (x' = x) (y' = y) (theta' = theta) (ct' = cos(theta)) " *
                        "(st' = sin(theta)) (vpx' <= $(vpx_max)) (vpy' <= $(vpy_max)));" for i in 1:n), "\n")

    function pushing_mode(i)
"""
// pushing on point $(i) = $(push_points[i])
{ mode $(i+1);
    invt:
        ct^2 + st^2 = 1;
    flow:
        d/dt[x] = ((c^2 + rox$(i)_w*rox$(i)_w)*vpx + rox$(i)_w*roy$(i)_w*vpy)/(c^2 + r$(i)sq);
        d/dt[y] = (rox$(i)_w*roy$(i)_w*vpx + (c^2 + roy$(i)_w*roy$(i)_w)*vpy)/(c^2 + r$(i)sq);
        d/dt[theta] = (rox$(i)_w*vpy - roy$(i)_w*vpx)/(c^2 + r$(i)sq);
        d/dt[ct] = -((rox$(i)_w*vpy - roy$(i)_w*vpx)/(c^2 + r$(i)sq))*st;
        d/dt[st] = ((rox$(i)_w*vpy - roy$(i)_w*vpx)/(c^2 + r$(i)sq))*ct;
        d/dt[vpx] = 0;
        d/dt[vpy] = 0;
    jump:
        (true) ==> @1 (and (x' = x) (y' = y) (theta' = theta) (ct' = cos(theta)) (st' = sin(theta)) (vpx' <= $(vpx_max)) (vpy' <= $(vpy_max)));
        (true) ==> @$(i+1) (and (x' = x) (y' = y) (theta' = theta) (ct' = cos(theta)) (st' = sin(theta)) (vpx' <= $(vpx_max)) (vpy' <= $(vpy_max)));
}
"""
    end

    pushing_modes = join((pushing_mode(i) for i in 1:n), "\n")

"""
#define c ($(P.c))
$body_coords
$world_coords

[0, $(t_max)] time;
// center of mass (x, y) bounds
[$(x_min), $(x_max)] x;
[$(y_min), $(y_max)] y;
[-20, 20] theta;
[-1, 1] ct;
[-1, 1] st;
// quasistatic pushing velocity bounds
[$(vpx_min), $(vpx_max)] vpx;
[$(vpy_min), $(vpy_max)] vpy;

// stationary
{ mode 1;
    invt:
        ct^2 + st^2 = 1;
    flow:
        d/dt[x] = 0;
        d/dt[y] = 0;
        d/dt[theta] = 0;
        d/dt[ct] = 0;
        d/dt[st] = 0;
        d/dt[vpx] = 0;
        d/dt[vpy] = 0;
    jump:
$jump_to_any
}

$pushing_modes

init:
@1  (and (x = 0) (y = 0) (theta = 0) (ct = 1) (st = 0));
goal:
@1  (and (x > 1));
"""

end

#### Actual problem generation

square = Polygon([[0,0], [1,0], [1,1], [0,1]])
triangle = Polygon([[0,0], [1,0], [0.5,sqrt(3)/2]])

open("dreach/push_test.drh", "w") do f
    write(f, pushing_problem(square, v_max=1.0, α_max = 10*pi/180, pin_θ=true))
end

open("dreach/push_test_tri.drh", "w") do f
    write(f, pushing_problem(triangle, v_max=1.0, α_max = 10*pi/180, pin_θ=true))
end