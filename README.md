# Fisher-Kolmogorov-equation-for-neurodegenerative-diseases
Fisher-Kolmogorov equation for neurodegenerative diseases

##Parameters
| **Test Case**                                     | **Time Step (\(\Delta t\))** | **Maximum or Final Time (T)** | **Reaction Coefficient (\(\alpha\))** | **Diffusion Tensor (D)**                         | **Comments**                                               |
|---------------------------------------------------|------------------------------|-------------------------------|-------------------------------------|-------------------------------------------------|------------------------------------------------------------|
| Test Case 1 (2D Convergence Analysis)             | \(10^{-5}\)                  | \(10^{-3}\)                   | 1                                   | \(D = d_{\text{ext}} I, d_{\text{ext}} = 1\)    | Short-term, high-precision simulations in 2D.              |
| Test Case 2 (Travelling Waves in 2D)              | 0.01                         | 5 and 10                      | 1                                   | \(D = d_{\text{ext}} I, d_{\text{ext}} = 1\)    | Modeling wave propagation in 2D with isotropic diffusion.  |
| Test Case 3 (α-Synuclein in 2D Brain Section)     | 0.01 years                   | Not specified                 | 0.9/year                            | \(D = d_{\text{ext}} I + d_{\text{axn}} (n \otimes n), d_{\text{ext}} = 8 \, \text{mm}^2/\text{year}, d_{\text{axn}} = 80 \, \text{mm}^2/\text{year}\) | Long-term disease progression with enhanced axonal diffusion. |
| Test Case 4 (3D Convergence Analysis)             | \(10^{-5}\)                  | \(10^{-3}\)                   | 0.1                                 | \(D = d_{\text{ext}} I, d_{\text{ext}} = 1\)    | 3D simulations with high precision and simplified parameters. |
| Test Case 5 (α-Synuclein Spreading in 3D Brain)   | 0.01 years                   | Not specified                 | 0.9/year                            | \(D = d_{\text{ext}} I + d_{\text{axn}} (n \otimes n), d_{\text{ext}} = 8 \, \text{mm}^2/\text{year}, d_{\text{axn}} = 80 \, \text{mm}^2/\text{year}\) | 3D brain model simulation with detailed diffusion modeling. |

This table provides a comprehensive overview of the simulation parameters for each test case, highlighting the time step, reaction coefficient, and diffusion tensor details, tailored to the specific objectives and scales of the various scenarios.
