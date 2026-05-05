#import "@preview/charged-ieee:0.1.4": ieee

#show: ieee.with(
  title: [Analysis of an X17-like Search],
  abstract: [
    Simulation and analysis of a simplified nuclear experiment sensitive to a hypothetical light boson with mass near 17 MeV/c2, often referred to as the X17. The purpose of the exercise is not only to study a possible beyond-Standard-Model signal, but also to learn the standard workflow used in high-energy and nuclear physics: Monte Carlo event generation, detector simulation with Geant4, data storage in ROOT format, and offline analysis with a small ROOT-based framework.
  ],
  authors: (
    (
      name: "D. Vázquez Lago",
      department: [Particle & Nuclear Physics],
      organization: [Universidad de Santiago de Compostela],
      location: [Santiago de Compostela, España],
      email: "daniel.vazquez.lago@rai.usc.es"
    ),
  ),
  bibliography: bibliography("refs.bib"),
  figure-supplement: [Fig.],
)

= Introduction

#lorem(650)

= Event generation

#figure(
  image("../plots/angular_distributions_ROOT.pdf", width: 100%),
  caption: [Angular distribution of the $e$ and $e^-$.],
  )

= Detector construction

== Detector in the ATOMKI geometry 

#figure(
  image("../plots/ATOMKI_geometry._0000.pdf", width: 100%),
  caption: [Detector in the ATOMKI geometry, with 20 events of $e^- e^+$ production.],
  )

== Optimization of the geometry 

= Transport and hits

= Analysis
