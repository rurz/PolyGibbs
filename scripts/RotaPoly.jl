# Preamble ########################################################################################
using Pkg
Pkg.activate(Base.current_project())

cd("..") # Go to the parent directory
###################################################################################################
# Packages needed #################################################################################
using PolyGibbs
using DelimitedFiles
using ImageCore
using ImageIO
using ImageMagick
using FileIO
using PyPlot
###################################################################################################
# Angle of rotation ###############################################################################

const θ = π/4

###################################################################################################
# Retrieve of data from image #####################################################################

@info "Retrieving the image channels and extracting the two-dimensional fields"
input = Float64.(channelview(RGB.(load("data/origin/wmap/wmap50.png"))))

rchn = input[1, :, :]
gchn = input[2, :, :]
bchn = input[3, :, :]

###################################################################################################
# Identification of the RGB fields dimensions ################################################

j = (size(input)[2] - 1)/2;

N = Integer(2*j + 1);

@info "The size of the field is $N"

###################################################################################################
# Rotation of the image channels ###########################################################################

@info "Performing rotation of R channel at θ = $θ"
routput = Emodes(rchn, θ)

@info "Performing rotation of G channel at θ = $θ"
goutput = Emodes(gchn, θ)

@info "Performing rotation of B channel at θ = $θ"
boutput = Emodes(bchn, θ)

###################################################################################################
# Saving results ##################################################################################

@info "Writting the results as separate data files for each channel"
writedlm("data/target/wmap/wmap50R_rpi4", real(routput))
writedlm("data/target/wmap/wmap50G_rpi4", real(goutput))
writedlm("data/target/wmap/wmap50B_rpi4", real(boutput))

###################################################################################################

begin
    imgrep = permutedims(StackedView(transpose(real(routput)), transpose(real(goutput)), transpose(real(boutput))), [3, 2, 1])
	fig = figure()
	ax = gca()
	ax.imshow(imgrep, cmap = "gray")
	tight_layout()
	savefig("data/target/wmap/wmap50_rpi4", dpi = 300)
end
