# Biomedical Imaging Simulations
Simulating biomedical imaging modalities (X-Ray, Positron Emission Tomography, Magnetic Resonance Imaging, Ultrasound Imaging) using MATLAB.

## Week 2
- Using MATLAB to calculate 2D plane waves (planewaves.mlx). The function takes as input the values of constants kx and ky and returns the 2D plane wave f(x,y).
- Another function (comb_function.mlx) is created which enables the user to create a comb function in discrete-time by specifying the number of impulses and the spacing Δx between them.

## Week 4
- Constructing a focused ultrasound beam by controlling an array of transducer elements. Ability to tune it according to the desired broadening angle and lateral resolution.

## Week 6
- Simulating X-ray planar projection imaging and digital subtraction angiography by using an analytical computer model of the human thorax.

## Week 7
- Emulating the process of Computed Tomography (CT) data generation (sinogram) and implementing various CT image reconstruction approaches such as simple backprojection (BP), filtered backprojection (FBP) and Fourier transform (FT) based reconstruction. Different filter settings in FBP and FT reconstruction and their impact on spatial resolution are also investigated.
- R2U function converts an image from the spatial domain to the spatial frequency domain.
- U2R function converts an image from the spatial frequency domain to the spatial domain.

## Week 8
- Extending the CT simulator to take into consideration realistic experimental parameters including quantum noise.
- Assessing image quality (Signal-to-Noise Ratio) after manipulating experimental parameters and filter settings.
- Studying the impact of radial undersampling and the presence of object motion on image quality.

## Week 9
- Reconstructing the Positron-Emission Tomography/Computed Tomography (PET/CT) in-vivo data.
- Implementing CT-based attenuation correction of PET phantom data by using the analytical computer model of the human thorax which is extended to include an artificial lung tumor accumulating a radio tracer.

## Week 10
- Implementing and studying kinetic modeling and data fitting for quantitative PET data analysis of <sup>18</sup>F-FDG tracer experiments.
