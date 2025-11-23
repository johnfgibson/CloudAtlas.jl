module CloudAtlasVisualizationExt
# Visualization tools
using CairoMakie
using CloudAtlas
using Statistics
using DelimitedFiles

export velocity_fields, velocity_fields_dns, velocity_fields_comparison, PlotSettings, DNSData, VelocityField

"""
    myreaddlm(filename)

Read delimited file with comments (lines starting with %).
"""
myreaddlm(filename) = readdlm(filename, comments=true, comment_char='%')

"""
    VelocityField

Callable struct that evaluates velocity components from a CloudAtlas model solution.
"""
struct VelocityField{T<:Real, M}
    Ψ::Vector{CloudAtlas.BasisFunction{T}}
    x::Vector{T}
    model::M
    add_baseflow::Bool
    
    function VelocityField(model::Union{CloudAtlas.ODEModel{T}, CloudAtlas.TWModel{T}}, x::Vector{T}; add_baseflow::Bool=false) where T<:Real
        new{T, typeof(model)}(model.Ψ, x, model, add_baseflow)
    end
end

"""
    (vf::VelocityField)(component::Symbol, x, y, z)

Evaluate velocity component at a point.
component ∈ [:u, :v, :w] for streamwise, wall-normal, spanwise velocity.
"""
function (vf::VelocityField)(component::Symbol, x::Real, y::Real, z::Real)
    idx = component == :u ? 1 : (component == :v ? 2 : 3)
    perturbation = sum(vf.Ψ[i].u[idx](x, y, z) * vf.x[i] for i in eachindex(vf.x))
    
    # Add baseflow for streamwise velocity if requested
    if vf.add_baseflow && component == :u
        return perturbation + y  # Linear profile: u_base = y
    else
        return perturbation
    end
end

# Convenience methods for all three components at once
function (vf::VelocityField)(x::Real, y::Real, z::Real)
    u = vf(:u, x, y, z)
    v = vf(:v, x, y, z)
    w = vf(:w, x, y, z)
    return (u, v, w)
end

"""
    plot_xz_plane!(ax, vf::VelocityField, settings::CloudAtlas.PlotSettings; Lx=2π, Lz=π, y_slice=0.0)

Plot (u, w) velocity field in the xz-plane at fixed y.
"""
function plot_xz_plane!(ax, vf::VelocityField, settings::CloudAtlas.PlotSettings; 
                        Lx=2π, Lz=π, y_slice=0.0)
    xs = range(0, Lx, settings.num_points)
    zs = range(0, Lz, settings.num_points)
    
    points = [Point2f(x, z) for x in xs for z in zs]
    
    if y_slice == :mean
        ys = range(-1, 1, settings.num_points)
        arrows_data = [(
            mean(vf(:u, x, y, z) for y in ys) * settings.arrow_scale,
            mean(vf(:w, x, y, z) for y in ys) * settings.arrow_scale
        ) for x in xs for z in zs]
    else
        arrows_data = [(
            vf(:u, x, y_slice, z) * settings.arrow_scale,
            vf(:w, x, y_slice, z) * settings.arrow_scale
        ) for x in xs for z in zs]
    end
    
    # Uses shaftwidth, tiplength, tipwidth
    arrows2d!(ax, points, arrows_data, 
            lengthscale=settings.arrow_lengthscale,
            tiplength=settings.arrow_tiplength,
            tipwidth=settings.arrow_tipwidth,
            shaftwidth=settings.arrow_shaftwidth)
end

"""
    plot_xy_plane!(ax, vf::VelocityField, settings::CloudAtlas.PlotSettings; Lx=2π, z_slice=0.0)

Plot (u, v) velocity field in the xy-plane at fixed z.
"""
function plot_xy_plane!(ax, vf::VelocityField, settings::CloudAtlas.PlotSettings; 
                        Lx=2π, z_slice=0.0)
    xs = range(0, Lx, settings.num_points)
    ys = range(-1, 1, settings.num_points)
    
    points = [Point2f(x, y) for x in xs for y in ys]
    arrows_data = [(
        vf(:u, x, y, z_slice) * settings.arrow_scale,
        vf(:v, x, y, z_slice) * settings.arrow_scale
    ) for x in xs for y in ys]
    
    arrows2d!(ax, points, arrows_data,
            lengthscale=settings.arrow_lengthscale,
            tiplength=settings.arrow_tiplength,
            tipwidth=settings.arrow_tipwidth,
            shaftwidth=settings.arrow_shaftwidth)
end

"""
    plot_yz_plane!(ax, vf::VelocityField, settings::CloudAtlas.PlotSettings; Lz=π, x_slice=0.0)

Plot (v, w) velocity field in the yz-plane at fixed x, with u heatmap background.
"""
function plot_yz_plane!(ax, vf::VelocityField, settings::CloudAtlas.PlotSettings; 
                        Lz=π, x_slice=0.0)
    zs = range(0, Lz, settings.num_points)
    ys = range(-1, 1, settings.num_points)
    
    # Background heatmap of u velocity
    u_field = [vf(:u, x_slice, y, z) for z in zs, y in ys]
    heatmap!(ax, zs, ys, u_field, colormap=settings.colormap)
    
    # Overlay (v, w) arrows
    points = [Point2f(z, y) for z in zs for y in ys]
    arrows_data = [(
        vf(:w, x_slice, y, z) * settings.arrow_scale,
        vf(:v, x_slice, y, z) * settings.arrow_scale
    ) for z in zs for y in ys]
    
    arrows2d!(ax, points, arrows_data,
            lengthscale=settings.arrow_lengthscale,
            tiplength=settings.arrow_tiplength,
            tipwidth=settings.arrow_tipwidth,
            shaftwidth=settings.arrow_shaftwidth,
            color=:black)
end

"""
    velocity_fields(model::Union{ODEModel, TWModel}, x::Vector; 
                    settings=PlotSettings(), Lx=2π, Lz=π, save_path=nothing, add_baseflow=false)

Create a three-panel visualization of velocity field (xz, xy, yz planes).

# Arguments
- `model`: ODEModel or TWModel instance
- `x`: Solution vector
- `settings`: PlotSettings for customization
- `Lx`: Domain length in x-direction (default: 2π)
- `Lz`: Domain length in z-direction (default: π)
- `save_path`: Optional path to save figure
- `add_baseflow`: If true, adds laminar Couette baseflow u = y (default: false)
"""
function CloudAtlas.velocity_fields(model::Union{CloudAtlas.ODEModel, CloudAtlas.TWModel}, x::Vector;
                        settings::CloudAtlas.PlotSettings=CloudAtlas.PlotSettings(),
                        Lx::Real=2π, Lz::Real=π,
                        save_path::Union{String,Nothing}=nothing,
                        add_baseflow::Bool=false)
    
    vf = VelocityField(model, x; add_baseflow=add_baseflow)
    fig = Figure(size=settings.fig_size)
    
    # XZ plane: mean (u, w) averaged over y
    ax1 = Axis(fig[1, 1],
        title="Mean (u, w) in xz-plane",
        xlabel="X", ylabel="Z",
        aspect=DataAspect())
    plot_xz_plane!(ax1, vf, settings, Lx=Lx, Lz=Lz, y_slice=:mean)
    
    # XY plane: (u, v) at z = 0
    ax2 = Axis(fig[2, 1],
        title="(u, v) in xy-plane at z = 0",
        xlabel="X", ylabel="Y",
        aspect=DataAspect())
    plot_xy_plane!(ax2, vf, settings, Lx=Lx, z_slice=0.0)
    
    # YZ plane: (v, w) at x = 0 with u heatmap
    ax3 = Axis(fig[3, 1],
        title="(v, w) in yz-plane at x = 0",
        xlabel="Z", ylabel="Y",
        aspect=DataAspect())
    plot_yz_plane!(ax3, vf, settings, Lz=Lz, x_slice=0.0)
    
    Colorbar(fig[3, 2], 
             colormap=settings.colormap,
             limits=(-1, 1),
             label="Streamwise velocity u")
    
    if save_path !== nothing
        CairoMakie.save(save_path, fig)
        println("Saved figure to: $save_path")
    end
    
    return fig
end

"""
    DNSData

Container for DNS velocity field data.
"""
struct DNSData{T<:Real}
    u_xy::Matrix{T}
    v_xy::Matrix{T}
    w_xy::Matrix{T}
    u_xz::Matrix{T}
    v_xz::Matrix{T}
    w_xz::Matrix{T}
    u_yz::Matrix{T}
    v_yz::Matrix{T}
    w_yz::Matrix{T}
    
    function DNSData(root::String)
        u_xy = myreaddlm(root * "_u_xy.asc")
        v_xy = myreaddlm(root * "_v_xy.asc")
        w_xy = myreaddlm(root * "_w_xy.asc")
        u_xz = myreaddlm(root * "_u_xz.asc")
        v_xz = myreaddlm(root * "_v_xz.asc")
        w_xz = myreaddlm(root * "_w_xz.asc")
        u_yz = myreaddlm(root * "_u_yz.asc")
        v_yz = myreaddlm(root * "_v_yz.asc")
        w_yz = myreaddlm(root * "_w_yz.asc")
        
        T = promote_type(eltype(u_xy), eltype(v_xy), eltype(w_xy))
        new{T}(u_xy, v_xy, w_xy, u_xz, v_xz, w_xz, u_yz, v_yz, w_yz)
    end
end

"""
    plot_dns_xz_plane!(ax, dns::DNSData, settings::CloudAtlas.PlotSettings; Lx=2π, Lz=π)

Plot DNS (u, w) field in xz-plane.
"""
function plot_dns_xz_plane!(ax, dns::DNSData, settings::CloudAtlas.PlotSettings; 
                            Lx::Real=2π, Lz::Real=π)
    nz, nx = size(dns.u_xz)
    xs = range(0, Lx, length=nx)
    zs = range(0, Lz, length=nz)
    
    points = [Point2f(x, z) for x in xs, z in zs] |> vec
    arrows_data = [(
        dns.u_xz[j, i] * settings.arrow_scale,
        dns.w_xz[j, i] * settings.arrow_scale
    ) for i in 1:nx, j in 1:nz] |> vec
    
    arrows2d!(ax, points, arrows_data,
            lengthscale=settings.arrow_lengthscale,
            tiplength=settings.arrow_tiplength,
            tipwidth=settings.arrow_tipwidth,
            shaftwidth=settings.arrow_shaftwidth)
end

"""
    plot_dns_xy_plane!(ax, dns::DNSData, settings::CloudAtlas.PlotSettings; Lx=2π)

Plot DNS (u, v) field in xy-plane at z=0.
"""
function plot_dns_xy_plane!(ax, dns::DNSData, settings::CloudAtlas.PlotSettings; Lx::Real=2π)
    ny, nx = size(dns.u_xy)
    xs = range(0, Lx, length=nx)
    ys = range(-1, 1, length=ny)
    
    points = [Point2f(x, y) for y in ys, x in xs] |> vec
    arrows_data = [(
        dns.u_xy[ny-j+1, i] * settings.arrow_scale,
        dns.v_xy[ny-j+1, i] * settings.arrow_scale
    ) for j in 1:ny, i in 1:nx] |> vec
    
    arrows2d!(ax, points, arrows_data,
            lengthscale=settings.arrow_lengthscale,
            tiplength=settings.arrow_tiplength,
            tipwidth=settings.arrow_tipwidth,
            shaftwidth=settings.arrow_shaftwidth)
end

"""
    plot_dns_yz_plane!(ax, dns::DNSData, settings::CloudAtlas.PlotSettings; Lz=π)

Plot DNS (v, w) field in yz-plane at x=0 with u heatmap.
"""
function plot_dns_yz_plane!(ax, dns::DNSData, settings::CloudAtlas.PlotSettings; Lz::Real=π)
    ny, nz = size(dns.v_yz)
    zs = range(0, Lz, length=nz)
    ys = range(-1, 1, length=ny)
    
    # Background heatmap
    heatmap!(ax, zs, ys, reverse(transpose(dns.u_yz)), 
             colormap=settings.colormap)
    
    # Overlay arrows
    points = [Point2f(z, y) for z in zs, y in ys] |> vec
    arrows_data = [(
        dns.w_yz[ny-i+1, j] * settings.arrow_scale,
        dns.v_yz[ny-i+1, j] * settings.arrow_scale
    ) for j in 1:nz, i in 1:ny] |> vec
    
    arrows2d!(ax, points, arrows_data,
            lengthscale=settings.arrow_lengthscale,
            tiplength=settings.arrow_tiplength,
            tipwidth=settings.arrow_tipwidth,
            shaftwidth=settings.arrow_shaftwidth,
            color=:black)
end

"""
    velocity_fields_dns(root::String; 
                        settings=PlotSettings(), Lx=2π, Lz=π, save_path=nothing)

Visualize DNS velocity field data.
"""
function CloudAtlas.velocity_fields_dns(root::String;
                            settings::CloudAtlas.PlotSettings=CloudAtlas.PlotSettings(),
                            Lx::Real=2π, Lz::Real=π,
                            save_path::Union{String,Nothing}=nothing)
    
    dns = DNSData(root)
    fig = Figure(size=settings.fig_size)
    
    # XZ plane
    ax1 = Axis(fig[1, 1],
        title="DNS (u, w) in xz-plane",
        xlabel="X", ylabel="Z",
        aspect=DataAspect())
    plot_dns_xz_plane!(ax1, dns, settings, Lx=Lx, Lz=Lz)
    
    # XY plane
    ax2 = Axis(fig[2, 1],
        title="DNS (u, v) in xy-plane at z = 0",
        xlabel="X", ylabel="Y",
        aspect=DataAspect())
    plot_dns_xy_plane!(ax2, dns, settings, Lx=Lx)
    
    # YZ plane
    ax3 = Axis(fig[3, 1],
        title="DNS (v, w) in yz-plane at x = 0",
        xlabel="Z", ylabel="Y",
        aspect=DataAspect())
    plot_dns_yz_plane!(ax3, dns, settings, Lz=Lz)
    
    Colorbar(fig[3, 2],
             colormap=settings.colormap,
             limits=(-1, 1),
             label="Streamwise velocity u")
    
    ny, nx = size(dns.u_xy)
    nz = size(dns.u_xz, 1)
    println("DNS grid: Nx × Ny × Nz = $nx × $ny × $nz")
    
    if save_path !== nothing
        CairoMakie.save(save_path, fig)
        println("Saved figure to: $save_path")
    end
    
    return fig
end

"""
    velocity_fields_comparison(model::Union{ODEModel, TWModel}, x::Vector, root::String;
                               settings=PlotSettings(), Lx=2π, Lz=π, save_path=nothing, add_baseflow=false)

Create side-by-side comparison plots of model solution and DNS data.
Saves three separate figures (xz, xy, yz) if save_path is provided.

# Arguments
- `model`: ODEModel or TWModel instance
- `x`: Solution vector
- `root`: Root path for DNS data files (without extension)
- `settings`: PlotSettings for customization
- `Lx`: Domain length in x-direction (default: 2π)
- `Lz`: Domain length in z-direction (default: π)
- `save_path`: Optional directory path to save figures
- `add_baseflow`: If true, adds laminar Couette baseflow u = y to model (default: false)
"""
function CloudAtlas.velocity_fields_comparison(model::Union{CloudAtlas.ODEModel, CloudAtlas.TWModel}, 
                                    x::Vector,
                                    root::String;
                                    settings::CloudAtlas.PlotSettings=CloudAtlas.PlotSettings(),
                                    Lx::Real=2π, Lz::Real=π,
                                    save_path::Union{String,Nothing}=nothing,
                                    add_baseflow::Bool=false)
    
    vf = VelocityField(model, x; add_baseflow=add_baseflow)
    dns = DNSData(root)
    
    # Use explicit keyword arguments to avoid positional errors
    comparison_settings = CloudAtlas.PlotSettings(
        num_points = settings.num_points,
        arrow_scale = settings.arrow_scale,
        arrow_lengthscale = settings.arrow_lengthscale,
        arrow_tiplength = settings.arrow_tiplength,
        arrow_tipwidth = settings.arrow_tipwidth,
        arrow_shaftwidth = settings.arrow_shaftwidth,
        colormap = settings.colormap,
        fig_size = (1920, 1080)  # Wider format for side-by-side
    )

    # Use explicit keyword arguments to avoid positional errors
    comparison_settings_vertical = CloudAtlas.PlotSettings(
        num_points = settings.num_points,
        arrow_scale = settings.arrow_scale,
        arrow_lengthscale = settings.arrow_lengthscale,
        arrow_tiplength = settings.arrow_tiplength,
        arrow_tipwidth = settings.arrow_tipwidth,
        arrow_shaftwidth = settings.arrow_shaftwidth,
        colormap = settings.colormap,
        fig_size = (1080, 1080)  # Vertical format for stacked layout
    )
    
    
    # === XZ PLANE ===
    fig_xz = Figure(size=comparison_settings.fig_size)
    ax_xz_ode = Axis(fig_xz[1, 1],
        title="Model (u, w) in xz-plane",
        xlabel="X", ylabel="Z",
        aspect=DataAspect())
    ax_xz_dns = Axis(fig_xz[1, 2],
        title="DNS (u, w) in xz-plane",
        xlabel="X", ylabel="Z",
        aspect=DataAspect())
    
    plot_xz_plane!(ax_xz_ode, vf, comparison_settings, Lx=Lx, Lz=Lz, y_slice=0.0)
    plot_dns_xz_plane!(ax_xz_dns, dns, comparison_settings, Lx=Lx, Lz=Lz)
    
    # === XY PLANE ===
    fig_xy = Figure(size=comparison_settings.fig_size)
    ax_xy_ode = Axis(fig_xy[1, 1],
        title="Model (u, v) in xy-plane",
        xlabel="X", ylabel="Y",
        aspect=DataAspect())
    ax_xy_dns = Axis(fig_xy[1, 2],
        title="DNS (u, v) in xy-plane",
        xlabel="X", ylabel="Y",
        aspect=DataAspect())
    
    plot_xy_plane!(ax_xy_ode, vf, comparison_settings, Lx=Lx, z_slice=0.0)
    plot_dns_xy_plane!(ax_xy_dns, dns, comparison_settings, Lx=Lx)
    
    # === YZ PLANE ===
    fig_yz = Figure(size=comparison_settings.fig_size)
    ax_yz_ode = Axis(fig_yz[1, 1],
        title="Model (v, w) in yz-plane",
        xlabel="Z", ylabel="Y",
        aspect=DataAspect())
    ax_yz_dns = Axis(fig_yz[1, 2],
        title="DNS (v, w) in yz-plane",
        xlabel="Z", ylabel="Y",
        aspect=DataAspect())
    
    plot_yz_plane!(ax_yz_ode, vf, comparison_settings, Lz=Lz, x_slice=0.0)
    plot_dns_yz_plane!(ax_yz_dns, dns, comparison_settings, Lz=Lz)
    
    Colorbar(fig_yz[1, 3],
             colormap=comparison_settings.colormap,
             limits=(-1, 1),
             label="Streamwise velocity u"
            )
    # Force the row to match the axis height
    rowsize!(fig_yz.layout, 1, Fixed(600))  # Adjust 600 for desired height...
    
    ny, nx = size(dns.u_xy)
    nz = size(dns.u_xz, 1)
    println("DNS grid: Nx × Ny × Nz = $nx × $ny × $nz")
    
    # Save if directory provided
    if save_path !== nothing
        mkpath(save_path)
        base_name = basename(root)
        CairoMakie.save(joinpath(save_path, "$(base_name)_xz_comparison.png"), fig_xz)
        CairoMakie.save(joinpath(save_path, "$(base_name)_xy_comparison.png"), fig_xy)
        CairoMakie.save(joinpath(save_path, "$(base_name)_yz_comparison.png"), fig_yz)
        println("Saved comparison figures to: $save_path")
    end
    
    return (fig_xz, fig_xy, fig_yz)
end
end