import Plots

function get_frame(frame_idx::Int64, trj::Trajectory)::Frame
    Frame(trj.coordinates[:,frame_idx], trj.attributes)
end

function get_atom(query::Int64, trj::Trajectory)::Vector{Atom}
    map(xyz -> Atom(xyz, trj.attributes[query]), trj.coordinates[query, :])
end

function get_atom(query::Int64, frame::Frame)::Atom
    Atom(frame.coordinates[query], frame.attributes[query])
end

function display_frame(frame::Frame; particle_size = 5)
    x_coords = map(xyz -> xyz[1], frame.coordinates)
    y_coords = map(xyz -> xyz[2], frame.coordinates)
    z_coords = map(xyz -> xyz[3], frame.coordinates)
    Plots.plotlyjs()
    Plots.plot(x_coords, y_coords, z_coords, zcolor = reverse(1:length(frame.coordinates)),
               marker = (particle_size, :coolwarm, Plots.stroke(0)), leg = false, w = 0)
end
