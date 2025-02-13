import argparse
from pathlib import Path

import tifffile as tf

from microsim import schema as ms

class SplittingDataMode:
    SingleChannel = 0
    MultiChannel = 1

# RATIO = 1.5


def run_simulation(dset: str, outputdir: str | None = None, 
                   splitting_mode:SplittingDataMode = SplittingDataMode.MultiChannel,
                   fluorophore1_str:str='EGFP',
                   fluorophore2_str:str='Venus',
                   label1:str='er-mem_pred',
                     label2:str='mito-mem_pred',
                    concentration1:float=1.0,
                    concentration2:float=1.0,
                   scale_factor:int=1,
                   shape_factor:float=1.0,
                   downscale:int=4,
                   exposure_ms:float=0.5,
                   ) -> None:
    print("Running simulation for", dset)
    assert scale_factor in [1, 2, 4, 8], "scale_factor must be one of 1, 2, 4, 8"
    flf1 = ms.Fluorophore.from_fpbase(fluorophore1_str)
    flf2 = ms.Fluorophore.from_fpbase(fluorophore2_str)
    if splitting_mode == SplittingDataMode.MultiChannel:
        sim = ms.Simulation(
            # note: this is a rather coarse simulation, but it's fast
            # scale should be a one of .004 * 2^n, where n is an integer from 0 to 4
            # space basically determines the field of view.
            truth_space=ms.ShapeScaleSpace(
                shape=(96, int(shape_factor*1400), int(shape_factor*1400)), scale=(0.004*scale_factor, 0.004*scale_factor, 0.004*scale_factor)
            ),
            output_space={"downscale": downscale},
            sample=[
                # pick dataset and layer name from https://openorganelle.janelia.org/datasets
                ms.FluorophoreDistribution(
                    distribution=ms.CosemLabel(dataset=dset, label=label1),
                    fluorophore=flf1,
                    concentration=concentration1
                ),
                ms.FluorophoreDistribution(
                    distribution=ms.CosemLabel(dataset=dset, label=label2),
                    fluorophore=flf2,
                    concentration=concentration2,
                ),
            ],
            channels=["i6WL::Widefield Dual Green", "i6WL::Widefield Triple Yellow"],
            modality=ms.Confocal(pinhole_au=1.5),
            detector=ms.CameraCCD(qe=0.82, read_noise=2),
            settings=ms.Settings(max_psf_radius_aus=2, cache=False),
            exposure_ms=exposure_ms,
            # output_path="bleedout.tif",
        )

    oipf = sim.optical_image_per_fluor()
    with_bleed = sim.digital_image(oipf.sum("f"))
    ch1 = sim.digital_image(oipf.sel(f=flf1))
    ch2 = sim.digital_image(oipf.sel(f=flf2))
    if outputdir is not None:
        print("Writing to files")
        dest = Path(outputdir) / dset
        dest.mkdir(parents=True, exist_ok=True)
        write_to_files(with_BT_data=with_bleed, ch1_data=ch1, ch2_data=ch2, dest=dest, dset=dset, 
                       concentration1=concentration1, concentration2=concentration2, 
                       scale_factor=scale_factor, downscale=downscale, fluor1_str=fluorophore1_str, 
                       fluor2_str=fluorophore2_str, exposure_ms=exposure_ms)
        # write_to_files(dest, dset, ratio, fluorophore1_str, fluorophore2_str, with_bleed, ch1, ch2)
    
    return dest, with_bleed, ch1, ch2

def fnames(dset:str, concentration1:float,concentration2:float, 
           scale_factor:int, downscale:int, fluor1_str:str, fluor2_str:str, exposure_ms:float):
    """
    R3.0_S4_D1_Ex100.0ms
    """
    postfix = f"{fluor1_str}_{fluor2_str}_R{concentration1}-{concentration2}_S{scale_factor}_D{downscale}_Ex{exposure_ms}ms"
    inp = f"{dset}_bleedthrough_{postfix}.tif"
    ch1 = f"{dset}_ch1_{postfix}.tif"
    ch2 = f"{dset}_ch2_{postfix}.tif"
    return inp, ch1, ch2

def write_to_files(with_BT_data, ch1_data, ch2_data, dest:Path, dset:str, 
                   concentration1:float,
                   concentration2:float,
                    scale_factor, downscale:int, fluor1_str:str, fluor2_str:str,
                   exposure_ms:float):
    inp_fname, ch1_fname, ch2_fname = fnames(dset=dset, 
                                             concentration1=concentration1,
                                                concentration2=concentration2, 
                                             scale_factor=scale_factor, 
                                             downscale=downscale, 
                                             fluor1_str=fluor1_str, 
                                             fluor2_str=fluor2_str,
                                             exposure_ms=exposure_ms)
    tf.imwrite(
        dest / inp_fname,
        with_BT_data.transpose("z", "c", "y", "x"),
        imagej=True,
    )
    tf.imwrite(dest / ch1_fname, ch1_data.isel(c=0), imagej=True)
    
    tf.imwrite(
        dest / ch2_fname, ch2_data.isel(c=1), imagej=True
    )



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--dset", type=str, default="jrc_hela-3")
    parser.add_argument("--label1", type=str, default="er-mem_pred")
    parser.add_argument("--label2", type=str, default="mito-mem_pred")
    parser.add_argument("--concentration1", type=float, default=1.0)
    parser.add_argument("--concentration2", type=float, default=1.0)
    parser.add_argument("--ratio", type=float, default=1.5, help="Relative concentration of fluorophore 2 to fluorophore 1")
    parser.add_argument("--fluor1", type=str, default="EGFP")
    parser.add_argument("--fluor2", type=str, default="Venus")
    parser.add_argument("--scale_factor", type=int, default=1)
    parser.add_argument("--downscale", type=int, default=4)
    parser.add_argument("--shape_factor", type=float, default=1.0)
    parser.add_argument("--outputdir", type=str, default="egfp_er_venus_mito")
    parser.add_argument("--exposure_ms", type=float, default=0.5)
    args = parser.parse_args()
    # dsets = CosemDataset.names()
    run_simulation(args.dset, outputdir=args.outputdir,
                   label1=args.label1, label2=args.label2,
                   concentration1=args.concentration1, concentration2=args.concentration2,
                   fluorophore1_str=args.fluor1, fluorophore2_str=args.fluor2,
                   scale_factor=args.scale_factor, downscale=args.downscale,
                   exposure_ms=args.exposure_ms,
                   shape_factor=args.shape_factor)
    # with ThreadPoolExecutor(max_workers=2) as pool:
    #     list(pool.map(run_simulation, dsets))
