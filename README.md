
inverse-mb-supp
===============

This repository contains supplementary code for

>  Inverse-Fourier Non-diffracting Beams for Optical Trapping,
>  Manuel Alberto Martınez-Ruiz, et al. [Journal TBD], [LINK TODO], 2020.

If you find any of this code useful, please consider citing our paper.
Some of this code will be included an upcoming release of the
[optical tweezers toolbox](https://github.com/ilent2/ott).

Requirements
------------

  * Optical Tweezers Toolbox (version 1, should work with 1.5.6)
  * Matlab (tested on R2018a)

Usage
-----

The repository consists of a series of example code files used
to generate different parts of the figures in the paper.
All the `example_` files should be self contained, simply ensure
you have the optical tweezers toolbox installed and run the file.

Other files contains functions used by the `example_` files.
Most of these files will be merged with the optical tweezers toolbox
at some later date.

Files
-----

  * `example_regular_mb_ottv1.m` Generate a regular Mathieu beam
    using the BscBessel class from OTTv1.

  * `example_inverse_mb.m` Generate an inverse Mathieu beam for use
    in OTT using the `BscPmMathieu` class.

  * `example_paraxial_mb.m` Example showing how paraxial Mathieu beams
    can be simulated.

  * `Mathieu.m` Functions used by `example_regular_mb_ottv1.m`, some
    are duplicated in `BscPmMathieu.m`.

  * `BscPmMathieu.m` generates an inverse Mathieu beam using point
    matching in the far-field.  Uses OTTv1 `BscPmParaxial` class.

