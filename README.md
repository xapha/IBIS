# IBIS

This repo. is dedicated to the IBIS super-pixel method.

This method was submitted at the IEEE ICPR 2018 conference.

#### Abstract : 

Superpixels are becoming increasingly popular for use in computer vision applications. However, most state-of-the-art superpixel segmentation algorithms suffer from a high computational cost, which prevents their use for real-time systems. We propose a new approach for superpixel segmentation with a low complexity, decreasing significantly the computation time needed. Instead of evaluating each pixel of an image to assign superpixel labels, this method named IBIS for Iterative Boundaries implicit Identification Segmentation, performs an iterative search for the superpixels boundaries in order to form regular and coherent clusters.
Our method is extremely fast while keeping a high precision of segmentation. Experimental results show that the segmentation can be performed much faster, up to 10 times compared to SLIC. Furthermore, the intrinsic parallelism of this method enables efficient implementations on GPU or any multi-core platforms to be proposed.

![alt text](https://github.com/xapha/IBIS/blob/master/intro.png "intro figure")