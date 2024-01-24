# Fisher-Kolmogorov-equation-for-neurodegenerative-diseases
Fisher-Kolmogorov equation for neurodegenerative diseases

<script type="text/javascript" async
  src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.7/MathJax.js?config=TeX-MML-AM_CHTML">
</script>

## Parameters

<div>

| **Test Case** | **Time Step <span>$$\Delta t$$</span>** | **Maximum or Final Time (T)** | **Reaction Coefficient <span>$$\alpha$$</span>** | **Diffusion Tensor (D)** | **Comments** |
| --- | --- | --- | --- | --- | --- |
| Test Case 1 (2D Convergence Analysis) | <span>$$10^{-5}$$</span> | <span>$$10^{-3}$$</span> | 1 | <span>$$D = d_{\text{ext}} I, d_{\text{ext}} = 1$$</span> | Short-term, high-precision simulations in 2D. |
| Test Case 2 (Travelling Waves in 2D) | 0.01 | 5 and 10 | 1 | <span>$$D = d_{\text{ext}} I, d_{\text{ext}} = 1$$</span> | Modeling wave propagation in 2D with isotropic diffusion. |
| Test Case 3 (α-Synuclein in 2D Brain Section) | 0.01 years | Not specified | 0.9/year | <span>$$D = d_{\text{ext}} I + d_{\text{axn}} (n \otimes n), d_{\text{ext}} = 8 \, \text{mm}^2/\text{year}, d_{\text{axn}} = 80 \, \text{mm}^2/\text{year}$$</span> | Long-term disease progression with enhanced axonal diffusion. |
| Test Case 4 (3D Convergence Analysis) | <span>$$10^{-5}$$</span> | <span>$$10^{-3}$$</span> | 0.1 | <span>$$D = d_{\text{ext}} I, d_{\text{ext}} = 1$$</span> | 3D simulations with high precision and simplified parameters. |
| Test Case 5 (α-Synuclein Spreading in 3D Brain) | 0.01 years | Not specified | 0.9/year | <span>$$D = d_{\text{ext}} I + d_{\text{axn}} (n \otimes n), d_{\text{ext}} = 8 \, \text{mm}^2/\text{year}, d_{\text{axn}} = 80 \, \text{mm}^2/\text{year}$$</span> | 3D brain model simulation with detailed diffusion modeling. |

</div>

This table provides a comprehensive overview of the simulation parameters for each test case, highlighting the time step, reaction coefficient, and diffusion tensor details, tailored to the specific objectives and scales of the various scenarios.
