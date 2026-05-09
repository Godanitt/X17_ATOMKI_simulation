#set document(title: [A Geant4--ROOT Feasibility Study of X17-like Pair Emission with Alternative Silicon--Scintillator Geometries], author: "Daniel Vázquez")
#set page(paper: "a4", margin: (left: 20mm, right: 20mm, top: 18mm, bottom: 20mm))
#set text(font: "New Computer Modern", size: 10.2pt, lang: "en")
#set par(justify: true, leading: 0.58em, first-line-indent: 0.0em)
#set heading(numbering: "1.1")
#set figure.caption(separator: [ -- ])
#set table(stroke: 0.35pt, gutter: 0pt)
#show heading.where(level: 1): it => block(above: 1.15em, below: 0.55em, it)
#show heading.where(level: 2): it => block(above: 0.85em, below: 0.35em, it)
#let MeVc = [MeV/$c^2$]
#let deg = [$degree$]
#let thetaee = [$theta_(e^+e^-)$]
#let root = smallcaps[ROOT]
#let geant = smallcaps[Geant4]

#align(center)[
  #text(size: 17pt, weight: "bold")[A Geant4--ROOT Feasibility Study of X17-like Pair Emission with Alternative Silicon--Scintillator Geometries]

  #v(0.45em)
  #text(size: 10.8pt)[Daniel Vázquez]

  #text(size: 9.2pt)[Master in Physics -- Universidade de Santiago de Compostela]

  #text(size: 9.2pt)[Técnicas de análisis y simulación en física nuclear y de partículas]

  #v(0.35em)
  #text(size: 8.8pt)[#datetime.today().display("[day] [month repr:long] [year]")]
]

#v(0.8em)

#block(stroke: 0.6pt + luma(130), inset: 10pt, radius: 2pt)[
  *Abstract.* This work presents a compact end-to-end feasibility study for an X17-like search in nuclear electron--positron pair emission.  A supplied event table is used as the generator-level signal model, while a smooth internal-pair-conversion-like background is generated with the same kinematic scale and propagated through the same reconstruction chain.  Four detector configurations are compared: an ATOMKI-like six-telescope layout, a downstream $2 pi$ semispherical shell, an approximate $4 pi$ double shell, and a forward pad plane.  The simulation chain is deliberately split into an ideal-transport stage, where #geant records only geometrical crossings, and an analysis-level detector-response stage, where finite efficiency, angular smearing, energy resolution and thresholds are applied before #root template fits are performed.  For equal generated signal and background samples of 10 000 events, the forward pad plane gives the largest reconstructed X17 purity in the present toy model, while the $4 pi$ shell maximizes total coincidence statistics but accepts substantially more smooth background.  The study is not a claim of sensitivity to a real X17 branching ratio; it is a controlled demonstration of how geometry, acceptance, detector effects and background modeling enter an X17-like angular-correlation analysis.

  #v(0.4em)
  *Github Respository*: https://github.com/Godanitt/X17_ATOMKI_simulation

  #v(0.4em)
  *Keywords:* X17, internal pair conversion, Geant4, ROOT, silicon detector, scintillator, template fit, geometrical acceptance.
]

#outline(title: [Contents], depth: 2)
#pagebreak()

= Introduction

Several anomalies in the angular correlation of electron--positron pairs emitted in light-nucleus transitions have motivated searches for a hypothetical boson with a mass close to 17 #MeVc, conventionally denoted X17.  The original #super[$8$]Be result reported an enhancement in the large-aperture-angle region of internal pair creation and interpreted it as a possible light neutral boson @krasznahorkay2016be8.  A widely discussed theoretical interpretation is a light protophobic vector boson, whose effective coupling to protons is suppressed relative to neutrons @feng2016protophobic.  Recent reviews emphasize both the interest of the anomaly and the need for independent confirmation: the X17 interpretation has been discussed for #super[$8$]Be, #super[$4$]He and #super[$12$]C transitions, but detector acceptance, resolution, nuclear-reaction modeling and analysis choices remain central systematic issues @gustavino2024x17 @alves2023shedding @krasznahorkay2024update.

The experimental observable of interest is the opening angle $theta_(e^+e^-)$, or equivalently the reconstructed invariant mass

$ m_(e e)^2 = 2 m_e^2 + 2 (E_+ E_- - |p_+| |p_-| cos theta_(e^+e^-)) . $

For a narrow boson decaying promptly into an $e^+e^-$ pair, the event population is expected to distort the smooth Standard-Model internal-pair-conversion (IPC) distribution, with the largest discriminating power appearing at large opening angle.  In the original ATOMKI-style interpretation the main irreducible background is IPC of virtual photons, which gives a smooth and generally falling aperture-angle spectrum rather than a narrow excess @gustavino2024x17.  More recent independent searches, including the MEG II study in #super[$7$]Li$(p,e^+e^-)$#super[$8$]Be, show that the experimental situation is still open and that robust detector-response studies are essential @megii2025x17.

The purpose of the present work is narrower and methodological.  We implement a reproducible simulation and analysis chain for an X17-like search using #geant @agostinelli2003geant4 @allison2016geant4 and #root @brun1997root.  The chain follows the standard workflow of the course exercise @tasfnp2026brief: event generation, geometry transport, hit storage, detector-effect parametrization, and offline analysis through histograms and fits.  The emphasis is on comparing detector geometries under identical input physics and identical reconstruction.  This is important because an apparent excess in $theta_(e^+e^-)$ is only meaningful after the acceptance and smearing of the spectrometer have been understood.

Four geometries are considered.  The `current` geometry is a simplified ATOMKI-like six-telescope arrangement in the transverse plane.  The `2pi` geometry is a downstream semispherical shell, excluding the beam-axis hole.  The `4pi` geometry mirrors the semisphere upstream and downstream to approximate nearly full angular coverage.  The `padplane` geometry is a forward plane of pads with the central beam pad removed.  The four layouts are not meant to be engineering-final designs; they are controlled alternatives that expose the trade-off between angular coverage and background acceptance.

= Detector model and geometries

== Sensitive volumes and materials

Each telescope is modeled as a thin double-sided silicon layer followed by a plastic scintillator block.  The #geant implementation uses a lithium target at the origin, air as the world material, silicon for the strip plane and vinyl-toluene plastic scintillator for the energy-measuring block.  The relevant dimensions are summarized in @tab:detector.

#figure(
  kind: table,
  table(
    columns: (1.25fr, 1.05fr, 1.05fr, 1.15fr),
    inset: 5pt,
    align: (left, right, right, left),
    [*Element*], [*Half thickness*], [*Face half-size*], [*Role*],
    [Li target], [0.05 mm], [5 mm radius], [interaction vertex],
    [Silicon], [0.25 mm], [4.1 cm $times$ 4.3 cm], [position and coincidence],
    [Plastic scintillator], [1.0 cm], [4.1 cm $times$ 4.3 cm], [energy response proxy],
    [Si--scintillator gap], [2.0 mm], [--], [mechanical separation],
  ),
  caption: [Detector elements used in the simplified geometry model.  The silicon layer is the default volume used for geometrical coincidence efficiencies.]
) <tab:detector>

The physics list is deliberately transport-only.  Electrons and positrons are generated at the target and transported in straight lines; no electromagnetic energy loss, multiple scattering, bremsstrahlung or positron annihilation is attached.  This choice makes the first simulation stage an acceptance calculation rather than a detailed detector simulation.  Detector response is introduced later at analysis level.  This separation is useful because it keeps the raw #geant output interpretable: an event is either geometrically accepted or not before any stochastic detector effect is imposed.

== Layouts under study

@fig:geo shows the four layouts used in the comparison.  The original ATOMKI-inspired layout has six telescope modules around the target.  The semispherical designs use three polar rings per hemisphere: six modules at 28#deg from the beam axis, ten modules at 55#deg and fourteen modules at 80#deg.  The `4pi` option mirrors this shell in the backward hemisphere.  The pad-plane option is a 5 by 5 downstream grid at $z=18$ cm with the central beam cell removed.

#figure(
  grid(
    columns: 2,
    gutter: 8pt,
    image("figures/geometry/ATOMKI._0000.pdf", width: 100%),
    image("figures/geometry/2pi._0000.pdf", width: 100%),
    image("figures/geometry/4pi_0000.pdf", width: 100%),
    image("figures/geometry/pad_plane_0000.pdf", width: 100%),
  ),
  caption: [Visualizations of the four detector configurations: ATOMKI-like six telescopes, downstream $2 pi$ shell, approximate $4 pi$ double shell and forward pad plane.]
) <fig:geo>

The designs intentionally keep a beam-axis hole.  In a real proton-beam experiment, instrumenting the beam direction with a sensitive block would be mechanically and physically problematic: it would intercept the incoming beam and introduce additional material and background.  The omission is therefore not a visualization accident but part of the geometry definition.

== Detector effects

After the ideal hit reconstruction, a compact detector-response model is applied to the coincidence tree.  The response parameters used in the present production are listed in @tab:effects.  They are not tuned to a particular apparatus; they provide a realistic first-order degradation of the ideal acceptance sample.

#figure(
  kind: table,
  table(
    columns: (1.45fr, 0.8fr, 1.9fr),
    inset: 5pt,
    align: (left, right, left),
    [*Parameter*], [*Value*], [*Meaning*],
    [Per-particle efficiency], [0.90], [independent survival of the $e^-$ and $e^+$],
    [Angular smearing], [2#deg], [Gaussian smearing in $theta$ and $phi$],
    [Energy resolution], [5%], [Gaussian relative kinetic-energy resolution],
    [Energy threshold], [1.0 MeV], [event rejected if either lepton is below threshold],
    [Signal region], [120º--180º], [large-angle region used for fit summaries],
  ),
  caption: [Analysis-level detector effects used after ideal geometrical hit reconstruction.]
) <tab:effects>

The implemented effects are the minimal effects needed for a credible first pass: finite particle efficiency, angular resolution, energy resolution and threshold.  Effects not included at this stage are equally important to state.  There is no timing model, so accidental coincidences are not simulated.  A stored event corresponds to an event in which both an electron and a positron are present in the same generated event and both survive the analysis-level detector response.  Electronic noise is not treated as an independent physics background; it is approximated only through smearing and threshold losses.  External pair conversion in passive material and multiple scattering are also absent because the transport stage intentionally uses a process-free physics list.

= Event generation and simulation chain

== Signal model

The X17-like signal is not generated from an analytic closed-form distribution.  Instead, the generator reads `data/data_pair_creation.txt`, whose five columns are

$ theta_(e e), quad theta_(e^-), quad T_(e^-), quad theta_(e^+), quad T_(e^+) . $

The event generator uses the table directly.  Event $i$ takes row $i$ modulo the number of rows, so requesting more events than rows simply wraps the input table.  The absolute azimuth is randomized because the file contains polar angles and energies but not an absolute azimuthal orientation.  The electron and positron are emitted on opposite azimuthal sides of the same plane, preserving the tabulated opening angle.  This choice avoids inventing a new signal distribution and makes the generator traceable to the supplied input data.

The sampled signal distributions are shown in @fig:signal-sampled.  The energy sum is effectively fixed at 17.5 MeV, while the individual lepton kinetic energies and polar angles vary event by event.  The generated invariant-mass distribution is therefore mostly driven by the supplied opening-angle structure and the energy sharing.

#figure(
  grid(
    columns: 2,
    gutter: 8pt,
    image("figures/sampling/signal_angles_sampled.pdf", width: 100%),
    image("figures/sampling/signal_energies_sampled.pdf", width: 100%),
  ),
  caption: [Generator-level signal sample.  Left: opening angle and lepton polar-angle distributions.  Right: total and per-lepton kinetic-energy distributions.]
) <fig:signal-sampled>

== IPC-like background model

The background is a smooth IPC-like toy model, not a tuned IPC calculation.  This distinction is important.  Its purpose is to supply a broad falling angular template that passes through the same detector and analysis chain as the signal.  The model uses a truncated exponential opening-angle distribution,

$ f(theta_(e e)) prop exp(-theta_(e e) / 40 degree), quad 0 <= theta_(e e) <= 180 degree, $

an isotropic absolute orientation, a fixed total kinetic energy of 17.5 MeV and a random energy sharing with both leptons above 0.20 MeV.  This captures the qualitative fact that IPC is smooth in opening angle and preferentially concentrated at smaller angles, without pretending to reproduce a full nuclear matrix-element calculation.  The sampled background distributions are shown in @fig:background-sampled.

#figure(
  grid(
    columns: 2,
    gutter: 8pt,
    image("figures/sampling/background_angles_sampled.pdf", width: 100%),
    image("figures/sampling/background_energies_sampled.pdf", width: 100%),
  ),
  caption: [Generator-level IPC-like toy background.  The opening angle follows a smooth falling law, while the absolute orientation is isotropic and the kinetic-energy sum is fixed to the signal scale.]
) <fig:background-sampled>

@fig:templates-sampled overlays the signal and background templates before geometry.  In the generated sample the signal populates larger opening angles and a narrower invariant-mass region than the smooth toy background.  The analysis problem is therefore not only to reconstruct a large-angle excess, but also to avoid changing the relative signal/background shapes through acceptance in a way that could mimic or hide such an excess.

#figure(
  grid(
    columns: 3,
    gutter: 6pt,
    image("figures/sampling/signal_background_thetaee_sampled.pdf", width: 100%),
    image("figures/sampling/signal_background_mass_sampled.pdf", width: 100%),
    image("figures/sampling/signal_background_energy_sum_sampled.pdf", width: 100%),
  ),
  caption: [Generated signal and background templates before detector acceptance.  The toy study is intentionally shape-based: the relative normalization of signal and background is set by the pseudo-data construction, not by a claimed physical branching ratio.]
) <fig:templates-sampled>

== ROOT data model and analysis flow

The raw #geant output is intentionally minimal.  The `sampled.root` and `background_sampled.root` files contain a `generated` tree and a `hits` tree.  The `generated` tree stores one row per generated pair event.  The `hits` tree stores one row per particle-volume crossing, with both electrons and positrons in the same tree.  Offline reconstruction groups hits by `eventID` and forms a detected event only when one $"pdg"=11$ hit and one $"pdg"=-11$ hit exist in the selected sensitive volume.  This definition is the coincidence criterion used throughout the report.

The analysis then writes `analysis_hits.root`, `background_analysis_hits.root`, detector-effect files, a final signal/background comparison file, and a template-fit file for each geometry.  The #root fit uses a histogram-template decomposition of the reconstructed $theta_(e e)$ distribution.  The pseudo-data tree keeps a truth flag only for validation; a real data analysis would not have access to this label.  The result is therefore an internal closure test of whether the known injected mixture can be recovered by the selected observable and templates.

= Geometrical acceptance

The first performance quantity is the ideal geometrical coincidence efficiency,

$ epsilon_("geom") = N_(e^- "&" e^+ "detected") / N_("generated"). $

The uncertainty in @tab:eff is the binomial standard error

$ sigma_epsilon = sqrt( epsilon_("geom") (1 - epsilon_("geom")) / N_("generated") ). $

The numbers in @tab:eff are produced by the reproducible `memoria/scripts/make_memoria_figures.py` cross-check using the silicon dimensions and geometry definitions.  In a ROOT-enabled environment, the companion macro `memoria/scripts/extract_root_summaries.C` extracts the official `summary` TTrees from `results/<geometry>/analysis_hits.root` and `results/<geometry>/background_analysis_hits.root` and writes replacement tables under `memoria/tables/`.

#figure(
  kind: table,
  include "tables/geometrical_efficiencies.typ",
  caption: [Ideal geometrical pair-coincidence efficiencies for 10 000 generated signal and 10 000 generated background events.  The uncertainty is binomial.]
) <tab:eff>

#figure(
  image("figures/sampling/geometrical_efficiencies.pdf", width: 75%),
  caption: [Comparison of signal and background coincidence acceptance for the four geometries.]
) <fig:eff>

The acceptance table exposes the essential design compromise.  The ATOMKI-like six-telescope geometry has significant acceptance for the isotropic smooth background but essentially no coincidence acceptance for the present forward-biased signal table.  This does not mean that an ATOMKI-like apparatus is intrinsically bad; it means that, with the present coordinate convention and supplied signal table, a transverse six-telescope layout is poorly matched to the generated signal directions.  The downstream $2 pi$ shell accepts roughly one quarter of the signal and a smaller background fraction.  The $4 pi$ shell increases total background acceptance strongly because it views both hemispheres, while the forward pad plane has the best signal-to-background acceptance ratio in this toy model.

= Reconstructed distributions and template fits

After detector effects, the analysis compares reconstructed signal and background distributions in $theta_(e e)$, summed kinetic energy and invariant mass.  The template-fit step uses the reconstructed opening angle as the primary observable because this is the observable in which IPC is expected to be smooth and the X17-like contribution is expected to appear as a localized or large-angle excess.

@fig:root-fits shows representative ROOT template fits for the $2 pi$, $4 pi$ and pad-plane geometries.  The current geometry is not shown because the present signal sample gives no meaningful reconstructed signal in the independent geometrical cross-check.  The three displayed fits should be read as closure plots rather than real-data claims: the fit is asked to recover the mixture of two samples whose labels are known only in the pseudo-data construction.

#figure(
  grid(
    columns: 3,
    gutter: 6pt,
    image("figures/rootfits/2pi/fit_thetaEE.pdf", width: 100%),
    image("figures/rootfits/4pi/fit_thetaEE.pdf", width: 100%),
    image("figures/rootfits/padplane/fit_thetaEE.pdf", width: 100%),
  ),
  caption: [ROOT template fits to the reconstructed $theta_(e e)$ pseudo-data distribution.  From left to right: $2 pi$, $4 pi$ and forward pad plane.]
) <fig:root-fits>

The same fitted fractions can be propagated to energy and invariant-mass projections, as shown in @fig:root-mass-energy for the pad-plane geometry.  These projections are not independent fits in the baseline chain; they are diagnostic checks of whether the mixture inferred from $theta_(e e)$ gives a reasonable description in related observables.

#figure(
  grid(
    columns: 2,
    gutter: 8pt,
    image("figures/rootfits/padplane/fit_energyEE.pdf", width: 100%),
    image("figures/rootfits/padplane/fit_massEE.pdf", width: 100%),
  ),
  caption: [Pad-plane fit projections in summed kinetic energy and reconstructed invariant mass.  The normalization is inherited from the $theta_(e e)$ template fit.]
) <fig:root-mass-energy>

For an equal-size generated signal and background toy sample, the expected reconstructed event counts after geometrical coincidence, per-particle efficiency and thresholds are summarized in @tab:reco.  These estimates are not physical rates; they are a controlled comparison of geometries under the same Monte Carlo inputs.  The final column is the fraction of reconstructed pseudo-data events that are true X17-like signal under that equal-sample convention.

#figure(
  kind: table,
  include "tables/reconstructed_yields.typ",
  caption: [Estimated reconstructed yields for equal generated signal and background samples of 10 000 events each.  The estimates include geometrical coincidence, 90% per-particle efficiency and the 1 MeV threshold.]
) <tab:reco>

The strongest conclusion from @tab:reco is that maximum coverage is not equivalent to maximum purity.  The $4 pi$ shell gives the largest total reconstructed event sample because it accepts many more background events, but its signal fraction is lower than in the $2 pi$ and pad-plane layouts.  The pad plane has the best reconstructed purity in this specific signal model because it is well matched to the forward signal kinematics and rejects much of the isotropic low-opening-angle background.  In a real experiment this conclusion would have to be re-evaluated with a physically normalized IPC rate, realistic beam and target backgrounds, timing, dead material, magnetic fields if present, and systematic uncertainties.

= Background subtraction strategy

The analysis strategy implemented in the project has two complementary forms.

First, the direct signal/background template fit constructs a pseudo-data sample by mixing reconstructed signal and background events, then fits the reconstructed $theta_(e e)$ distribution with normalized signal and background templates.  The fitted signal yield $N_X^"fit"$ and background yield $N_b^"fit"$ are obtained from the fitted fractions multiplied by the total pseudo-data normalization.  This is the cleanest closure test because the template shapes are generated from the same simulation and detector-response chain as the pseudo-data.

Second, the project includes an IPC-like smooth-background-plus-excess fit.  In this approach the smooth component is described by a falling empirical function and the excess by a localized component.  This is closer in spirit to a real search, where the background may be constrained by sidebands or by a smooth theory-driven IPC model while the signal is inferred as an excess.  It is also more fragile: if the background parametrization is too flexible, it can absorb part of the signal; if it is too rigid, it can create a fake excess.  For this reason the template closure test is the more robust baseline, and the smooth-background excess fit should be treated as a diagnostic cross-check rather than the final measurement.

The subtraction logic can be written schematically as

$ N_X(theta) = N_("data")(theta) - alpha N_("bkg")(theta), $

where $alpha$ is either determined from the template fit or from a control-region normalization.  The uncertainty in the subtracted spectrum contains the data counting variance, the background template variance and the uncertainty on $alpha$.  In this simplified study the dominant reported uncertainties are binomial acceptances and template-fit uncertainties from ROOT.  A complete analysis would add nuisance parameters for efficiency, angular resolution, energy scale, target thickness and background-shape choices.

= Discussion

The study makes three points that are useful for the next iteration of the project.

First, detector geometry dominates the result before any sophisticated fitting is attempted.  If the generated signal kinematics are forward-biased, a detector concentrated around the transverse plane can have poor coincidence acceptance.  Conversely, a forward pad plane may perform well for that same input while being less suitable for a signal model with substantial backward or transverse population.  Geometry optimization should therefore be performed against the actual physics distribution, not only against nominal solid angle.

Second, background rejection is not the same problem as signal acceptance.  The $4 pi$ shell gives the broadest coverage and is attractive if the goal is maximum statistics or angular-response completeness.  However, it also accepts an isotropic IPC-like background efficiently.  The pad plane, in the present toy model, wins because it accepts the forward signal more efficiently than the smooth background.  This is an analysis-dependent conclusion, not a universal detector principle.

Third, the current simulation is intentionally honest about its limitations.  The signal model comes from the supplied table; the background is a smooth IPC-like toy model; the detector response is applied after ideal transport; no timing, accidental coincidences, external pair conversion, multiple scattering or detailed scintillation light collection is modeled.  These omissions are acceptable for a first feasibility paper because they isolate geometry and reconstruction.  They would not be acceptable for a claim of discovery, exclusion or physical branching-ratio sensitivity.

= Conclusions

A complete geant4--root workflow has been organized for an X17-like pair-emission study with four detector geometries and a reproducible Typst report.  The workflow separates generator-level physics, ideal geometrical transport, analysis-level detector response and ROOT-based template fitting.  This separation makes the result auditable: the origin of each distortion in the reconstructed $theta_(e e)$ spectrum can be traced to either the input physics, the geometry, the response model or the fit.

Within the present equal-sample toy study, the forward pad plane is the best geometry for the supplied signal table.  It gives the largest estimated reconstructed signal yield and the highest true X17 fraction after detector effects.  The $2 pi$ shell is a reasonable compromise with lower total statistics but a cleaner signal/background balance than the $4 pi$ layout.  The $4 pi$ shell is valuable for acceptance completeness but, for the smooth isotropic background used here, it admits too much background to be the purity-optimal design.  The current ATOMKI-like transverse layout is not well matched to the present input coordinate convention and should either be reinterpreted with the correct physical axes or reserved as a historical/reference geometry.

The most important next improvements are clear: replace the toy IPC background by a physically normalized IPC generator, include passive-material effects and multiple scattering, add a timing model for accidental coincidences, validate the geometry with overlap-free placements, and treat detector-response parameters as nuisance parameters in the fit.  With those additions, the same code structure can evolve from a feasibility demonstration into a quantitatively defensible sensitivity study.

#bibliography("refs.bib", style: "american-physics-society")
