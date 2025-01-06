from microsim.cosem import CosemDataset

# can find all the names with CosemDataset.names()
dset_name = "jrc_hela-3"
# fetch the dataset
dataset = CosemDataset.fetch(dset_name)

# get all image names
# here you could skip if it doesn't have the label
print([x.name for x in dataset.images])

# load a single label image
img = dataset.image(name="er-mem_seg")

# get some data, paying attention to bin_mode.
# where bin_mode:
# - 'standard' is what cosem does (based on the mode of instance segmentation)\
# - 'sum' is more like fluorophore counts (more often useful)
# - 'auto' (default) is 'standard' if it's a prediction, 'sum' if it's a segmentation

# get the highest resolution level (level 0) of the image
# because there IS no binning here, the bin_mode is irrelevant
# however to avoid the bug we discussed, I use bin_mode='standard'
data = img.read(level=0, bin_mode="standard")

# show the pixel scale (this is the nanometers per pixel in the data)
print(img.grid_scale)


# level 0 data will have a scale of 0.004 um/px
# to save locally (figure it out)


# THIS level here is what I would probably augment ...
# then downsample to something like 0.032 um/px
# then pass to from_ground_truth
your_augmentation = ...
your_downscaled_augmentation = ...  # use multiples of 3 ... see below


# sim = Simulation.from_ground_truth(
#     your_downscaled_augmentation,
#     # tell the simulation the scale after you've downscaled
#     scale=(0.032, 0.032, 0.032),
#     # after the full optical PSF simulation we need to simulate camera pixelation
#     # 3 is a decent choice because it will yield a 96nm pixel size, which
#     # is approximately nyquist 1.4NA that is the DEFAULT of the simulation
#     # (change only if you change the NA of the objective.)
#     # currently: downscale of 3 will likely fail if `your_downscaled_augmentation` is
#     # not factor of 3 smaller than the original image
#     output_space={"downscale": 3},
#     modality=ms.Confocal(),
#     detector=ms.CameraCCD(qe=0.82, read_noise=2, bit_depth=12, offset=100),
# )
