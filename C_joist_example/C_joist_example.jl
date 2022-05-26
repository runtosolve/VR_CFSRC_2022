using ThinWalledBeam, CrossSection, Plots, CUFSM


#generate a C joist cross-section
#https://sfia.memberclicks.net/assets/TechFiles/SFIA%20Tech%20Spec%202022%20%20%202.7.22%20Final.pdf

#800S200-54, from the SFIA catalog
t = 0.0566 #in.

#define outside dimensions
L = [0.625, 2.0, 8.0, 2.0, 0.625] #in.

#define cross-section element orientations
θ = [π/2, π, -π/2, 0.0, π/2]

#define outside bend radii
r = [0.0849 + t, 0.0849 + t, 0.0849 + t, 0.0849 + t]

#define cross-section discretization
n = [2, 2, 2, 2, 2]
n_r = [3, 3, 3, 3]

#outside surface of cross-section, discretized into nodes
cross_section = CrossSection.generate_thin_walled(L, θ, n, r, n_r)

#shift cross-section nodes so that bottom flange is at y=0, outside web is at x=0
X = [cross_section[i][1] for i in eachindex(cross_section)]
Y = [cross_section[i][2] for i in eachindex(cross_section)]
cross_section = [[cross_section[i][1].- minimum(X),cross_section[i][2].- minimum(Y)] for i in eachindex(cross_section)]
plot(X, Y, markershape = :o, aspect_ratio = :equal, seriestype = :scatter)

#calculate surface normals at nodes
unit_node_normals = CrossSection.Tools.calculate_cross_section_unit_node_normals(cross_section)

#get centerline cross-section coordinates
Δ = -t/2
centerline = CrossSection.Tools.get_coords_along_node_normals(cross_section, unit_node_normals, Δ)

#calculate section properties from centerline nodes
X = [centerline[i][1] for i in eachindex(centerline)]
Y = [centerline[i][2] for i in eachindex(centerline)]
coord = [X Y]
num_cross_section_nodes = size(coords, 1)
ends = [1:num_cross_section_nodes-1 2:num_cross_section_nodes t * ones(Float64, num_cross_section_nodes-1)]
section_props = CUFSM.cutwp_prop2(coord, ends)

#########
#define inputs for joist structural analysis

span_length = 120.0 #in.
num_beam_segments = 12
top_flange_width = 2.0 #in.
joist_depth = 8.0 #in.
z = 0.0:span_length/num_beam_segments:span_length
num_joist_nodes = length(z)

Ix = ones(Float64, num_joist_nodes) * section_props.Ixx
Iy = ones(Float64, num_joist_nodes) * section_props.Iyy
Ixy = ones(Float64, num_joist_nodes) * section_props.Ixy
J = ones(Float64, num_joist_nodes) * section_props.J
Cw = ones(Float64, num_joist_nodes) * section_props.Cw

E = ones(Float64, num_joist_nodes) * 29500.0 #ksi
ν = ones(Float64, num_joist_nodes) * 0.3
G = E./(2 .*(1 .+ ν))

#distance from joist shear center to applied load (load is applied at the middle of top flange)
ax = ones(Float64, num_joist_nodes) * (abs(section_props.xs) + top_flange_width/2)
ay = ones(Float64, num_joist_nodes) * (joist_depth - abs(section_props.ys))

#assume no joist lateral bracing
kx = ones(Float64, num_joist_nodes) * 0.0
ay_kx = ones(Float64, num_joist_nodes) * 0.0
kϕ = ones(Float64, num_joist_nodes) * 0.0

#define uniform load on the joist (will do point load once you get this running Benjy)
qx = ones(Float64, num_joist_nodes) * 0.0  #no lateral force applied to joist
qy = ones(Float64, num_joist_nodes) * 0.1 #kips/in., positive is downward

#u''=v''=ϕ''=0 (simply-supported), u'=v'=ϕ'=0  (fixed), u''=v''=ϕ''=u'''=v'''=ϕ'''=0 (free, e.g., a cantilever)
end_boundary_conditions = ["simply-supported", "simply-supported"]

#supports
#z, u, v, ϕ  fixed means fixed in translation
supports = [(0.0, "fixed", "fixed", "fixed"),
            (span_length, "fixed", "fixed", "fixed")]

#run the analysis
model = ThinWalledBeam.solve(z, Ix, Iy, Ixy, J, Cw, E, G, kx, kϕ, ay_kx, qx, qy, ax, ay, end_boundary_conditions, supports)


#plot vertical deflection
plot(model.inputs.z, -model.outputs.v)

#plot lateral deflection
plot(model.inputs.z, -model.outputs.u)

#plot twist
plot(model.inputs.z, model.outputs.ϕ)



