# Functions

## Index
```@index
Order = [:function]
```

## File IO
```@docs
read_dcd
write_dcd
read_pdb
write_pdb
read_xyz
write_xyz
```

## Trajectory handling
```@docs
get_frame
get_atom
clip_trajectory
center_of_mass
radius_of_gyration
pair_length_matrix
pair_length_matrix_pbc
contact_bool_matrix
contact_probability_matrix
```

## Utility
```@docs
atom_mass
residue_mass
distance_pbc
```

### Multi-Threading
```@docs
pair_length_matrix_parallel
contact_bool_matrix_parallel
contact_probability_matrix_parallel
```
