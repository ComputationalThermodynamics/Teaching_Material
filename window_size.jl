#=
window_size.jl  –  open the window with a size that fits on the screen

The figures of the examples are designed for a window of about 1500 × 950.
On a smaller screen such a window does not fit. Before the window is opened,
the examples therefore call

    fit_to_screen!(fig)

which makes the figure smaller if it is larger than the fraction
WINDOW_FRACTION of the screen. Nothing is scaled: text, sliders and buttons
keep their size, the diagrams get the space that is left. To get smaller (or
larger) windows, change WINDOW_FRACTION below.
=#

using GLMakie

const WINDOW_FRACTION = (0.85, 0.80)    # largest window at start-up, as a fraction of the screen (width, height)

"""
    screen_size()

Size (width, height) of the primary screen in the units of `Figure(size = ...)`,
or `nothing` if it cannot be found out.
"""
function screen_size()
    try
        GLFW    = GLMakie.GLFW
        monitor = GLFW.GetPrimaryMonitor()
        monitor.handle == C_NULL && return nothing
        mode = GLFW.GetVideoMode(monitor)
        # macOS and Wayland give the screen size in the units of a window; Windows and X11 give
        # pixels, of which a window has `scale` per unit ("150 %" -> 1.5)
        scale = GLFW.GetPlatform() in (GLFW.PLATFORM_COCOA, GLFW.PLATFORM_WAYLAND) ? 1.0 :
                Float64(minimum(GLFW.GetMonitorContentScale(monitor)))
        return (mode.width / scale, mode.height / scale)
    catch
        return nothing
    end
end

"""
    fit_to_screen!(fig; fraction = WINDOW_FRACTION, screen = screen_size())

Make `fig` smaller if it is wider or higher than `fraction` of the screen
(width and height are limited separately; the figure is never made larger).
Call it before `display(fig)`: the window opens with the size of the figure.
"""
function fit_to_screen!(fig; fraction = WINDOW_FRACTION, screen = screen_size())
    screen === nothing && return fig
    w0, h0 = size(fig.scene)
    w = min(w0, round(Int, fraction[1] * screen[1]))
    h = min(h0, round(Int, fraction[2] * screen[2]))
    (w, h) == (w0, h0) || resize!(fig, w, h)
    return fig
end
