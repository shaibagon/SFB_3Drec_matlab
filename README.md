## 3D Reconstruction of a Moving Object

 Re-implementation of

 Simakov, Frolova and Basri *"Dense Shape Reconstruction from a Moving
 Oject under Arbitrary, Unknown Lighting"*, 2003


 This re-imp is restricted to linear lighting model (eq. 10 in the  paper)


### Usage:
````
   z = SimakovFrolovaBasri2003(Io, mo, Ij, Mj, Tj, p)

 Inputs:
   Io  - reference image, 2D matrix of intensity values
   mo  - mask for reference image, 2D binary image of the same size as Io.
         Depth is estimated only for pixels in the mask.
   Ij  - cell array of corresponding images (2D matrices)
   Mj  - cell array of corresponding masks (2D matrices) same size as Ij's
   Tj  - struct array of transformations, each element has
         ().R  - rotation matrix w.r.t reference frame (3x3 matrix)
         ().t  - translation vector w.r.t reference frame (2x1 vector)
         ().scale - scaling scalar w.r.t reference frame
   p   - parameters structure
         ().SFB_or_BC - using SFB consistency measure, or
                        brightness constancy (BC)
         ().z_range   - [min_z max_z] range of labels
         ().NL        - number of labels (100 - 500)
         ().L1_truncate - percent of Z range (.25 - .5)
         ().lambda    - smoothness cost weight (roughly 1e-4 for 'bc', 5e-3 for 'sfb')

 Output:
   z   - depth values for pixels in Io
   e   - energy of solution
````


### Terms of Use

 Copyright (c) Bagon Shai
 Department of Computer Science and Applied Mathmatics
 Wiezmann Institute of Science
 http://www.wisdom.weizmann.ac.il/

 Permission is hereby granted, free of charge, to any person obtaining a copy
 of this software (the "Software"), to deal
 in the Software without restriction, subject to the following conditions:

 The above copyright notice and this permission notice shall be included in
 all copies or substantial portions of the Software.

 The Software is provided "as is", without warranty of any kind.



 If used in an academic research, the follwoing citation must be included
 in any resulting publication:

>   [1] Denis Simakov, Darya Frolova and Ronen Basri
>       "Dense Shape Reconstruction from a Moving Oject under Arbitrary, Unknown Lighting",
>       ICCV, 2003

>   [2] Shai Bagon
>       "Matlab Implementation of Simakov Frolova and Basri 3D Reconstruction",
>       June, 2012


 June. 2012


