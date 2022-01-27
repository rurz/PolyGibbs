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

@info "Retrieving the image and extracting the two-dimensional field"
input = Float64.(channelview(Gray.(load("data/origin/symm/symm50.png"))))

###################################################################################################
# Identification of the grayscale field dimensions ################################################

j = (size(input)[1] - 1)/2;

N = Integer(2*j + 1);

@info "The size of the field is $N"

###################################################################################################
# Rotation of the image ###########################################################################

@info "Performing rotation at θ = $θ"

output = Emodes(input, θ)

###################################################################################################
# Saving results ##################################################################################

@info "Writting the results"
writedlm("data/target/symm/symm50_rpi4.dat", real(output))

###################################################################################################

begin
	fig = figure()
	ax = gca()
	ax.imshow(real(output), cmap = "gray")
	tight_layout()
	savefig("data/target/symm/symm50_rpi4.png", dpi = 300)
end
