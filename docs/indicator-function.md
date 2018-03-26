# Actuator disk indicator function

Each wind turbine actuator disk with finite thickness \\(s\\) is represented by
the normalized indicator function

\\[ \mathcal{I}(\mathbf{x}) = V^{-1}\left[ H(\hat{x}+s/2) - H(\hat{x} - s/2)  \right] H(D/2 - \hat{r}) \\]
where \\(V = s \pi D^2/4\\) is the volume of the disk, \\(H(x)\\) is the
Heaviside function. The coordinate system of the actuator disk is denoted by
hats with the \\(\mathbf{\hat{e}_1}\\) unit vector opposite of the thrust force
and \\(\hat{r}^2 = \hat{y}^2 + \hat{z}^2\\). To avoid Gibbs phenomena around
the actuator disk forcing, a smoothed indicator function

\\[ \mathcal{R}(\mathbf{x}) = \int G(\mathbf{x}-\mathbf{x'}) \, \mathcal{I}(\mathbf{x'}) \, d^3\mathbf{x'} \\]

is used to apply the thrust force on the flow and average the velocity field
over the disk. The filtering is applied using a Gaussian kernel

\\[ G(\mathbf{x}) = \left(\frac{6}{\pi \Delta}\right)^{3/2} \exp \left( -\frac{6\lVert\mathbf{x}\rVert^2}{\Delta^2} \right)  \\]

with a filter width \\(\Delta = \alpha h \\) that is based on the grid size
\\(h = \sqrt{dx^2 + dy^2 + dz^2}\\).

The smoothed indicator function can be decomposed
\\[ \mathcal{R}(\mathbf{x}) = V^{-1} \mathcal{R}_1(\hat{x}) \mathcal{R}_2(\hat{r})\\]
into a normal component
\\[ \mathcal{R}_1(\hat{x}) = \left(\frac{6}{\pi \Delta}\right)^{1/2} \int \left[ H(\hat{x}'+s/2) - H(\hat{x}' - s/2) \right]    \exp \left( -6\frac{(\hat{x}-\hat{x}')^2}{\Delta^2} \right) \, d \hat{x}' \\]
and radial component
\\[ \mathcal{R}_2(\hat{r}) = \frac{6}{\pi \Delta} \int \int H(D/2 - \hat{r}')  \exp \left( -6\frac{(\hat{y}-\hat{y}')^2 + (\hat{z}-\hat{z}')^2}{\Delta^2} \right) \, d \hat{y}' \, d \hat{z}'. \\]
LESGO calculates the the smoothed indicator function at each point using the
analytic solution to the normal component

\\[ \mathcal{R}_1(\hat{x}) = \frac{1}{2}\mathrm{erf}\left(\frac{\sqrt{6}}{\Delta}\left(x+\frac{s}{2}\right) \right) - \frac{1}{2}\mathrm{erf}\left(\frac{\sqrt{6}}{\Delta}\left(x-\frac{s}{2}\right) \right)\\]

and numerically computing \\(\mathcal{R}_2(\hat{r})\\) on a grid finer than the simulation grid.

