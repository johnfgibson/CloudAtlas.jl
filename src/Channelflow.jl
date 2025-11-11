module Channelflow

export projectfield, field_to_coeff, coeff_to_field, changegrid, findsoln, continuesoln

import Channelflow_jll
using DelimitedFiles: readdlm
# include("BasisFunctions.jl")

"""
    myreaddlm(filename, cc='%')

Read matrix or vector from a file, dropping comments marked with cc.
"""
function myreaddlm(filename, cc='%')
    X = readdlm(filename, comments=true, comment_char=cc)
    if size(X,2) == 1
        X = X[:,1]
    end
    X
end

"""
    verify_file(file)

Throws an error if `file` doesn't exist. Better to have Julia error instead of a downstream binary.
"""
function verify_file(file::AbstractString)
    if !ispath(file)
        throw(ArgumentError(file, "File does not exist!"))
    end
end

"""
    kwargs_to_flags(kwargs)

Converts a list of keyword arguments to a string of flags to pass into a Channelflow binary.
"""
function kwargs_to_flags(kwargs)
    # Convert keyword arguments into a vector of strings
    flags = String[]
    for (key, value) in kwargs
        # Use the key as the flag (e.g., :Nx -> "-Nx")
        flag_name = "-$(key)"

        if value isa Bool
            # Handle boolean flags: only add the flag if true
            value && push!(flags, flag_name)
        else
            # Handle key-value pairs (e.g., "-Nx", "64")
            push!(flags, flag_name, string(value))
        end
    end
    return flags 
end

"""
    field_to_coeff(ijklfile, source, output; kwargs...)

`field_to_coeff` takes a `source` (full velocity field) and projects it onto the
ODE basis defined by `ijklfile`, storing the resulting coefficient vector in `output`.

In addition to writing to the `output` file, this returns the vector of ODE coefficients.
This calls the Channelflow binary `projectfield` in its default mode.

# Arguments
- `ijklfile::AbstractString`: The file path to the file specifying the ijkl basis.
- `source::AbstractString`: The file path to the input velocity field to be projected.
- `output_base::AbstractString`: The file path to store the output coefficient vector. 
                                 `.asc` is appended to this if it doesn't already end with it.

# Keyword Arguments
- `nrm::Bool`: Normalize basis elements.
- `sb::Bool`: Output basis elements as flowfields.
- `sIP::Bool`: Output the inner product matrix.
- `pb::Bool`: Save slices of basis elements for plotting.
- `pn::Bool`: Print L2Norm(psi[n]), L2Norm(div(u)), etc.
- `o::String`: Directory for basis plots (default: "plots/").
- `A::String`: File path to load inner product matrix from.
"""
function field_to_coeff(ijklfile::AbstractString, source::AbstractString, output_base::AbstractString; kwargs...)
    verify_file(ijklfile)
    verify_file(source)
    output = occursin(".asc", output_base) ? output_base : output_base * ".asc"

    flags = kwargs_to_flags(kwargs)
    run(`$(Channelflow_jll.projectfield()) $flags $ijklfile $source $output`)
    mv("x" * output, output)
    return myreaddlm(output)
end

"""
    field_to_coeff(ijkl::AbstractMatrix, source, output; kwargs...)

Converts a field to the coefficient basis described by `ijkl`, where `ijkl` is a
Julia matrix rather than a file.
"""
function field_to_coeff(ijkl::AbstractMatrix, source::AbstractString, output::AbstractString; kwargs...)
    # This function uses temp files and calls the other method.
    ijkl2file(ijkl, "temp_ijkl.asc")
    coeffs = field_to_coeff("temp_ijkl.asc", source, output; kwargs...)
    rm("temp_ijkl.asc")
    return coeffs
end

"""
    coeff_to_field(source_coeffs, ijklfile, field_example, output; kwargs...)

Converts a file containing ODE coefficients (`source_coeffs`) back into a full
velocity field (`output`), using `field_example` as a grid template.

This calls the Channelflow binary `projectfield` with the `-x` flag.

# Arguments
- `source_coeffs::AbstractString`: File path to the vector of ODE coefficients.
- `ijklfile::AbstractString`: File path to the file specifying the ijkl basis.
- `field_example::AbstractString`: File path to a velocity field to use as a template
                                 (for grid size, box length, etc.).
- `output::AbstractString`: The file path to store the reconstructed output velocity field.

# Keyword Arguments
- All flags from `field_to_coeff` are also valid here, e.g., `nrm=true`
"""
function coeff_to_field(source_coeffs::AbstractString, ijklfile::AbstractString, field_example::AbstractString, output::AbstractString; kwargs...)
    verify_file(source_coeffs)
    verify_file(ijklfile)
    verify_file(field_example)

    flags = kwargs_to_flags(kwargs)

    run(`$(Channelflow_jll.projectfield()) $flags -x $source_coeffs $ijklfile $field_example $output`)
end

"""
    coeff_to_field(x, ijkl, field_example, output; kwargs...)

Converts a coefficient vector `x` and `ijkl` matrix directly, without needing to
store them in files. The result is stored in a file at `output`.
"""
function coeff_to_field(x::AbstractVector, ijkl::AbstractMatrix, field_example::AbstractString, output::AbstractString; kwargs...)
    # This function uses temp files and calls the other method.
    ijkl2file(ijkl, "temp_ijkl.asc")
    save(x, "temp_x.asc")

    coeff_to_field("temp_x.asc", "temp_ijkl.asc", field_example, output; kwargs=kwargs)

    rm("temp_ijkl.asc")
    rm("temp_x.asc")
end

"""
    projectfield(ijklfile, source, target; flags)

Legacy interface, defer to `field_to_coeff` ideally.
"""
function projectfield(ijklfile, source, target; kwargs...)
    field_to_coeff(ijklfile, source, target; kwargs=kwargs)
end

"""
    projectfield(source_coeffs, ijklfile, field_example, output; flags)

Legacy interface, defer to `coeff_to_field` ideally.
"""
function projectfield(source_coeffs, ijklfile, field_example, output; kwargs...)
    coeff_to_field(source_coeffs, ijklfile, field_example, output; kwargs=kwargs)
end


"""
    findsoln(guess_flowfield::AbstractString; kwargs...)

Takes the flowfield stored at filepath `guess_flowfield` and attempts to find a
solution (e.g., equilibrium, traveling wave, periodic orbit) using a Newton-Krylov solver.

This is a Julia wrapper for the Channelflow binary `findsoln`.

# Arguments
- `guess_flowfield::AbstractString`: The file path to the initial guess flow field.

# Keyword Arguments
All command-line flags for the `findsoln` binary are passed to this function as keywords.
- Boolean switches (e.g., `-eqb`) are passed as `eqb=true`.
- Value flags (e.g., `-R 350`) are passed as `R=350`.
- String flags (e.g., `-solver gmres`) are passed as `solver="gmres"`.

See the `findsoln --help` output for a complete list.

## Common Keyword Examples

### Search Type:
- `eqb::Bool`: Search for a fixed point (equilibrium) or traveling wave.
- `orb::Bool`: Search for a periodic or relative periodic orbit.
- `xrel::Bool`: Search for a relative solution (TW or RPO) with shift in x.
- `zrel::Bool`: Search for a relative solution (TW or RPO) with shift in z.

### System Parameters:
- `R::Real`: Pseudo-Reynolds number (default: 400).
- `nu::Real`: Kinematic viscosity (overrides `R`).
- `symms::String`: Path to a file listing symmetry generators.
- `sigma::String`: Path to a file for the `sigma` matrix (for TWs/RPOs).

### Numerical/Solver Parameters:
- `T::Real`: Final time of integration or period of map (default: 20).
- `dt::Real`: Timestep (default: 0.03125).
- `solver::String`: "gmres", "fgmres", "eigen", "bicgstab".
- `opt::String`: "hookstep", "linear", "none".
- `Nn::Int`: Max Newton steps (default: 20).
- `Ng::Int`: Max GMRES iterations (default: 500).
- `od::String`: Output directory (default: "./").

# Example Usage
```julia
# Find an equilibrium at Re=350 with symmetries from local file
findsoln("u_Re350.nc"; eqb=true, R=350, T=10, symms="./sxy_sz_txz.asc")

# Find a z-relative traveling wave
findsoln("u_guess.nc"; eqb=true, zrel=true, T=10, R=200, sigma="sigma_guess.asc", symms="sxytxz.asc")
```
"""
function findsoln(guess_flowfield::AbstractString; kwargs...)
    verify_file(guess_flowfield)
    flags = kwargs_to_flags(kwargs)
    run(`$(Channelflow_jll.findsoln()) $flags $guess_flowfield`)
end


"""
    changegrid(infield::AbstractString, outfield::AbstractString; kwargs...)

Interpolates a given flowfield onto a different grid.
This calls the Channelflow binary `changegrid`.

# Arguments
- `infield::AbstractString`: The file path to the input flow field.
- `outfield::AbstractString`: The file path to store the output (interpolated) flow field.

# Keyword Arguments
- `p::Bool` (default: true): Set padding modes to zero.
- `dv::Bool` (default: true): Fix divergence and Dirichlet BCs.
- `Nx::Int`: New number of x gridpoints.
- `Ny::Int`: New number of y gridpoints.
- `Nz::Int`: New number of z gridpoints.
- `a::Real`: New lower wall position.
- `b::Real`: New upper wall position.
- `al::Real` or `alpha::Real`: New alpha (2pi/Lx).
- `ga::Real` or `gamma::Real`: New gamma (2pi/Lz).
- `lx::Real`: New Lx = 2 pi lx.
- `lz::Real`: New Lz = 2 pi lz.
- `Lx::Real`: Streamwise (x) box length.
- `Lz::Real`: Spanwise (z) box length.
- `np0::Int` or `nproc0::Int`: Number of MPI-processes for transpose.
- `np1::Int` or `nproc1::Int`: Number of MPI-processes for one fft.
"""
function changegrid(infield::AbstractString, outfield::AbstractString; kwargs...)
    verify_file(infield)
    flags = kwargs_to_flags(kwargs) 
    run(`$(Channelflow_jll.changegrid()) $flags $infield $outfield`)
end

"""
    continuesoln(initial_flowfield::AbstractString; kwargs...)

Performs a parametric continuation of an invariant solution (e.g., equilibrium,
traveling wave, or periodic orbit) starting from an initial guess.

This is a Julia wrapper for the Channelflow binary `continuesoln`.

# Arguments
- `initial_flowfield::AbstractString`: The file path to the initial solution
  (e.g., "ueqd_Re400.nc") from which to start the continuation.

# Keyword Arguments
All command-line flags for the `continuesoln` binary are passed as keywords.
- Boolean switches (e.g., `-eqb`) are passed as `eqb=true`.
- Value flags (e.g., `-R 400`) are passed as `R=400`.
- String flags (e.g., `-cont Re`) are passed as `cont="Re"`.

See `continuesoln --help` for the full list of ~70 options.

## Key Keyword Arguments
- `cont::String`: **(Required)** The continuation parameter. One of:
  "Re", "P", "Ub", "Uw", "ReP", "Theta", "ThLx", "ThLz", "Lx", "Lz",
  "Aspect", "Diag", "Lt", "Vs", "ReVs", "H", "HVs", "Rot".
- `dmu::Real`: The initial relative increment for the continuation parameter
  (default: 0.0001).
- `eqb::Bool` or `orb::Bool`: Specify the solution type (equilibrium or orbit).
- `R::Real`: The initial pseudo-Reynolds number (default: 400).
- `symms::String`: Path to a file listing symmetry generators.
- `al::Bool`: Use arclength continuation (default: false).
- `ds0::Real`: Initial arclength increment (if `al=true`).
- `targ::Bool`: Abort when the target value `targMu` is reached.
- `targMu::Real`: The target value for the continuation parameter.
- `od::String`: Output directory (default: "./").

# Example Usage
```julia
# Continue an equilibrium in Reynolds number
continuesoln("ueqd_Re400.nc";
    cont="Re",
    eqb=true,
    R=400,
    T=19,
    dmu=-0.01,
    symms="sxy_sztx.asc"
)
```
"""
function continuesoln(initial_flowfield::AbstractString; kwargs...)
    verify_file(initial_flowfield)
    flags = kwargs_to_flags(kwargs)
    run(`$(Channelflow_jll.continuesoln()) $flags $initial_flowfield`)
end

end