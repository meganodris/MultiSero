MultiSero
================

This method uses multi-pathogen serological data to simultaneously
estimate the prevalence of antibodies against each pathogen and
between-pathogen cross-reactive relationships. Sources of multi-pathogen
serological data can come from multiplex binding immunoassays and/or
separate immunoassays (e.g. ELISAs) that have been run on the same
samples. To account for cross-reactivity among related pathogens,
antibodies should be measured against the same antigen across pathogens
(e.g. envelope protein for flaviviruses).

## Concept: antibodies as multivariate Gaussian components

We use a simulated example to demonstrate the model concept, where
antibody titers against two pathogens, A and B, have been measured in a
study population. With only two pathogens, we can characterize this data
using two-dimensional Gaussian components. Assuming both pathogens have
transmitted in this population there are four possible infection
statuses (or antibody profiles) that we will characterize using Gaussian
components:

- Negative (grey)
- A positive (red)
- B positive (blue)
- A and B positive (purple)

The scatter plots below show how the characteristics of these four
Gaussian components change in scenarios of varying levels of
cross-reactivity, with log antibody titer measurements against pathogens
A and B shown on the x and y-axes. Here, $\phi$ is the cross-reactivity
parameter defined as the relative titer increase (e.g. binding) induced
by the infecting pathogen against the other pathogen(s). A value of
$\phi=0$ indicates no recognition/binding of the antibodies to the other
pathogen and $\phi=1$ indicates identical recognition/binding of the
antibodies to both pathogen. All other parameters are assumed constant
across scenarios (i)-(iv). The histograms on the top and right axes show
the single-dimension antibody titer distributions for pathogen A and
pathogen B respectively. Dashed lines show the means of the Gaussian
components. We describe each of the scenarios (i-iv) below.

![](ReadMeFiles/ConceptPlot.png)<!-- -->

**Scenario (i): No cross-reactivity**

In this baseline scenario we assume pathogens A and B to be unrelated
(i.e. no cross-reactivity). The mean titer of the A-positive component
(red) against pathogen B is equal to that of the negative component
(grey horizontal line), while the mean titer of the B-positive component
(blue) against pathogen A is also equal to that of the negative
component (grey vertical line).

**Scenario (ii): Cross-reactivity from A to B**

In this scenario we assume that antibodies induced by infection with
pathogen A partially recognize pathogen B ($\phi_{AB}=0.25$), but
antibodies induced from infection by pathogen B do not recognize
pathogen A ($\phi_{BA}=0$). Here, the Gaussian component that is
A-positive (red) has shifted upwards on the y-axis, shown by the
horizontal red line, due to the these antibodies now partially
recognizing/binding pathogen B. The mean titer of these antibodies
against pathogen B is now 0.5 (calculated as shown below), higher than
the mean of the negative component (grey). Here, $\mu_{0,B}$ is the mean
of the negative titers against pathogen B and $\mu_{1,A}$ is the mean
titer increase against pathogen A induced by infection with pathogen A.
The increase in common binding also proportionally increases the
correlation and covariance of the A-positive (red) component.

$$\mu_{0,B} + \phi_{AB} * \mu_{1,A}$$ $$= 0 + 0.25*2$$ $$= 0.5$$

**Scenario (iii): Cross-reactivity from B to A**

We now assume that antibodies induced by infection with pathogen B
recognize pathogen A ($\phi_{BA}=0.5$), while antibodies induced from
infection by pathogen A do not recognize pathogen B ($\phi_{AB}=0$).
Here, the Gaussian component that is B-positive (blue) has shifted right
on the x-axis, shown by the vertical blue line, due to these antibodies
partially recognizing/binding pathogen A. The mean antibody titer of
this component against pathogen A is now 1.25 (calculated as shown
below), higher than the mean of the negative component (grey). Here,
$\mu_{0,A}$ is the mean of the negative titers against pathogen A and
$\mu_{1,B}$ is the mean titer increase against pathogen B induced by
infection with pathogen B. The increase in common binding also
proportionally increases the correlation and covariance of the
B-positive (blue) component.

$$\mu_{0,A} + \phi_{BA} * \mu_{1,B}$$ $$= 0 + 0.5*2.5$$ $$= 1.25$$

**Scenario (iv): Cross-reactivity in both directions**

In this last scenario we assume cross-reactivity in both directions,
where antibodies induced from infection by each pathogen recognize the
other pathogen to varying extents. The A-positive (red) Gaussian
component now shifts up on the y-axis (due to binding against pathogen
B) and the B-positive (blue) Gaussian component shifts right on the
x-axis due to binding against pathogen A. The means of these components
are calculated in the same way as scenarios (ii) and (iii).

## Applying the model

The model is written in Stan and executed using cmdStan. It can be
applied to any number of pathogens, though computational time increases
significantly as we consider more pathogens. As the number of pathogens
is increased, the dimensions of the Gaussian components increases and
the number of possible infection statuses also increases.
