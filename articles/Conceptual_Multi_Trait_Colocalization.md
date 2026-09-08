# Conceptual Framework for Multi-trait Colocalization and ColocBoost

This vignette introduces key conceptual definitions for multi-trait
colocalization, describes their practical implementation in ColocBoost,
and relates them to established concepts in fine-mapping and existing
colocalization methods.

## 1. Multi-trait colocalization problem as a colocalization event decomposition

The multi-trait colocalization problem concerns identifying *shared
causal signals* across $`L \geq 2`$ traits within a genomic region of
interest. To formalize this problem, ColocBoost introduces a new
conceptual definition of a multi-trait ***Colocalization Event*** as a
triplet $`\{ g(s), T(s), CoS_{\alpha}(s) \}`$ for each detected event
$`s`$, comprising

- $`g(s)`$: ***Genomic Region*** underlying the colocalization analysis.
- $`T(s)`$: ***Trait Configuration*** specifying which subset of $`L`$
  traits share the same causal variant for event $`s`$.
- $`CoS_{\alpha}(s)`$: ***Colocalization Confidence Set*** representing
  the smallest set of variants containing the shared causal variant with
  probability at least $`\alpha`$ (default $`\alpha = 0.95`$).

A central challenge in multi-trait colocalization is that multiple
events can arise within the same genomic region, involving overlapping
or distinct sets of traits and multiple causal variants.

- In Bayesian formulations that specify all possible trait
  configurations a *priori*, the hypothesis space grows combinatorially
  and rapidly becomes computationally prohibitive.
- ColocBoost avoids this enumeration by reformulating colocalization as
  a multi-task learning problem optimized through gradient boosting.
  This formulation jointly resolves multiple colocalization events and
  scales to hundreds of traits.

![Illustration of colocalization events in
ColocBoost.](figures/Colocalization_Events.png)

## 2. Conceptual ColocBoost summaries and their analogies to existing methods

ColocBoost characterizes each event using two complementary summaries:
variant-level localization and event-level evidence for sharing across
the corresponding trait configuration. These summaries have direct
structural analogies to established quantities in fine-mapping and
colocalization, while being defined specifically for multi-trait
colocalization events.

- **Variant-level evidence:** The ***colocalization confidence set***
  (CoS) and ***variant colocalization probability*** (VCP) localize the
  variants underlying each colocalization event. They are structurally
  analogous to the credible set (CS) and posterior inclusion probability
  (PIP), respectively, in statistical fine-mapping methods such as
  SuSiE.

- **Colocalization evidence:** The ***normalized probability of
  colocalization*** (NPC) quantifies support for colocalization. It is
  structurally analogous to PP.H4 in COLOC for pairwise colocalization
  and PPFC in HyPrColoc for multi-trait colocalization.

### 2.1. Variant-level evidence

For each event $`s`$, representing a single causal signal shared by the
subset of traits $`T(s)`$,

- $`\xi^s = (\xi_1^s, \ldots, \xi_P^s)`$ denotes the vector of
  *single-effect colocalization probabilities* across $`P`$ variants in
  the region, quantifying the probability that each variant is the
  shared causal variant underlying traits $`T(s)`$.
- *Structural analogy*: single-effect posterior $`\alpha_l`$ in SuSiE
  for single-effect $`l`$ in a single-trait model.

ColocBoost then defines an $`\alpha`$-level *Colocalization Confidence
Set*, $`CoS_{\alpha}(s)`$ (default $`\alpha = 0.95`$), including the
candidate causal variant and its high-LD proxies, with construction
based on $`\xi^s`$:

``` math
  CoS_{\alpha}(s) = \left\{  v_1, v_2, \ldots, v_{p_0} : p_0 = \min (p: \sum_{j=1}^p \xi_{(j)}^s \geq \alpha )  \right\},
```

where $`\xi_{(1)}^s \geq \ldots \geq \xi_{(P)}^s`$ are the sorted
single-effect colocalization probabilities. ColocBoost discards the CoS
with low *purity* (minimum absolute correlation between all pairs of
variants within CoS, default threshold $`purity < 0.5`$).

- *Structural analogy*: $`\alpha`$-level credible set (CS) in SuSiE in a
  single-trait model.

ColocBoost also defines the *variant colocalization probability* (VCP)
for each variant $`j`$ in the region as
``` math
  VCP_j = 1 - \prod_{s=1}^{S} (1 - \xi_j^s).
```
The construction of VCP is based on the assumption that each event $`s`$
is *conditionally* independent.

- *Structural analogy*: posterior inclusion probability (PIP) in SuSiE
  in a single-trait model.
- *Structural analogy*: variant-level posterior probability for H4
  (SNP.PP.H4) in COLOC in a pairwise colocalization analysis.

![Illustration of variant-level analogues in
ColocBoost.](figures/Variant_Level_Analogues.png)

### 2.2. Colocalization evidence

#### Narrative definition of normalized probability of colocalization (NPC)

For each event $`s`$, ColocBoost defines the *normalized probability of
colocalization* (NPC) as an empirical, event-level measure comparing
shared colocalization with trait-specific alternatives.

- Higher NPC values indicate stronger evidence that the event represents
  a causal signal shared by at least two traits.
- Lower NPC values indicate weaker evidence for sharing, while greater
  consistency with the single-trait specific causal signals.

NPC is *structurally analogous*, but not probabilistically equivalent,
to

- PP.H4: posterior probability of H4, two traits shared the same causal
  variant, in COLOC for pairwise colocalization.
- PPFC: posterior probability of full colocalization in HyPrColoc for
  multi-trait colocalization.

NPC is assigned to each detected colocalization event $`s`$, allowing
multiple distinct events within the same genomic region and a
potentially different trait configuration $`T(s)`$.

**Note:** PP.H4 is defined **only** for two traits and denotes posterior
probability for $`H_4`$, the hypothesis that both traits are associated
and share the same causal variant, relative to hypotheses $`H_0`$,
$`H_1`$, $`H_2`$, and $`H_3`$. For $`L > 2`$, the colocalization
hypothesis space must encompass all possible trait configurations across
shared and distinct causal variants (Foley et al., 2021, *Nature
Communications*). For example, $`H_{(L-2,1,1)}`$ represents a
configuration in which $`L-2`$ traits share one causal variant, while
the remaining two traits have distinct causal variants. The full
hypothesis space grows combinatorially, comprising
$`\mathrm{Bell}(L+1)`$ hypotheses.

ColocBoost avoids this combinatorial explosion by discovering supported
configurations $`T(s)`$ and their corresponding $`CoS_{\alpha}(s)`$ in a
data-driven manner, without enumerating all possible trait
configurations (details in ColocBoost paper).

#### Mathematical definition of normalized probability of colocalization (NPC)

For a detected colocalization event $`s`$ with a triplet
$`\{ g(s), T(s), CoS_{\alpha}(s) \}`$, ColocBoost first defines the
trait-level normalized evidence score $`NP_l^s`$ for each trait
$`l \in T(s)`$:

``` math
  NP_l^s = 1 - \exp(-\lambda_l LRT_l^s),
```

where $`\lambda_l`$ is a trait-specific rate that adjusts the scale of
normalization based on a baseline log-likelihood ratio for each trait
$`l`$ (details in Supplementary Note of ColocBoost paper).

Here, $`LRT_l^s`$ is a log-likelihood ratio test statistic between two
models, $`M^s_{0,l}`$ and $`M^s_{1,l}`$.

- $`M^s_{0,l}`$: the null model where all variants have zero effects on
  trait $`l`$ ($`\beta_l=0`$).
- $`M^s_{1,l}`$: the alternative model that variants in
  $`CoS_{\alpha}(s)`$ have non-zero effects ($`\beta^s_l \neq 0`$).

***Rationale:*** NPC is evaluated only after event $`s`$ has been
detected; therefore, at least one trait is expected to provide strong
evidence for the event. NPC provides event-level evidence that the
detected event is jointly supported by at least two traits, rather than
being driven by evidence from only one trait.

**Two-trait example**

For the two-trait case, without loss of generality, let
$`NP_1^s \ge NP_2^s`$, so that trait 1 is the leading trait for event
$`s`$. ColocBoost approximates the evidence for the single-trait,
non-colocalized configuration, in which trait 1 contributes but trait 2
does not, as

``` math
NPUC_s
=
\underbrace{NP_1^s}_{
\substack{\text{evidence supporting}\\
\text{trait 1 contribution}}
}
\times
\underbrace{(1-NP_2^s)}_{
\substack{\text{lack of evidence supporting}\\
\text{trait 2 contribution}}
}.
```

Accordingly, *a large NPUC indicates that event $`s`$ is supported
**primarily** by the leading trait*. ColocBoost subsequently defines the
event-level colocalization evidence NPC as

``` math
  NPC_s = 1 - NPUC_s.
```

NPUC and NPC are normalized colocalization evidence scores rather than
posterior probabilities.

**Multi-trait generalization**

For $`L>2`$, ColocBoost generalizes NPUC as the *leading-trait-only
explanation* across all traits:

``` math
NPUC_s = NP^s_{l_{\max}} \prod_{l\neq l_{\max}} \left(1-NP_{l}^{s}\right), \, \mathrm{and} \, NPC_s = 1 - NPUC_s.
```

Here, $`NP^s_{l_{\max}}`$ represents the normalized evidence from the
leading trait, whereas the product term captures the lack of support
from the remaining traits. Together, these terms quantify the extent to
which event $`s`$ is supported primarily by the leading trait.
Accordingly, a high NPC indicates support from more than one trait.

**Highlight:** ColocBoost provides complementary evidence at two levels:
$`NPC_s`$ evaluates the overall colocalization event, whereas $`NP_l^s`$
quantifies each trait’s support for that event. In our numerical
studies, the lenient thresholds $`NPC_s \geq 0.5`$ and
$`NP_l^s \geq 0.2`$ maintained well-controlled FDR while retaining
reasonable detection power. More stringent thresholds may be applied to
prioritize the strongest colocalization signals.

See more details about filtering colocalization events by relative
strength of evidence using
[‘get_robust_colocalization’](https://statfungen.github.io/colocboost/articles/Interpret_ColocBoost_Output.html#filter-colocalization-events-by-relative-strength-of-evidence)
function.

![Illustration of event-level analogues in
ColocBoost.](figures/Event_Level_Analogues.png)

Concordance between NPC and PP.H4 or PPFC was assessed only for detected
CoS with at least 95% overlap variants between ColocBoost and COLOC or
HyPrColoc, respectively.

## 3. Practical interpretation for `colocboost` output

This section maps the concepts above to the corresponding `colocboost`
output fields. See [Interpret ColocBoost
Output](https://statfungen.github.io/colocboost/articles/Interpret_ColocBoost_Output.html)
for detailed guidance.

After running `res = colocboost()`,

- `res$cos_summary`: a summary of all colocalization events. Each row
  corresponds to one colocalization event $`s`$ and includes columns
  with
  - `colocalized_outcomes`: *Trait Configuration* – $`T(s)`$;
  - `colocalized_variables`: *Colocalization Confidence Set* –
    $`CoS_{\alpha}(s)`$;
  - `colocalized_variables_vcp`: *Variant Colocalization Probability*
    for shared variants in $`CoS_{\alpha}(s)`$;
  - `cos_npc`: *Normalized Probability of Colocalization* – $`NPC_s`$;
  - `purity`: minimum absolute correlation between all pairs of variants
    within $`CoS_{\alpha}(s)`$.
- `res$cos_details`: detailed information about all colocalization
  events, including sublists with
  - `vcp`: a length-$`P`$ vector of *Variant Colocalization Probability*
    for all $`P`$ variants;
  - `cos_vcp`: a list of single-effect colocalization probabilities for
    each event $`s`$;
  - `cos_outcomes_npc`: a list of trait-level evidence ($`NP^s_l`$) for
    each event $`s`$.

### Example: Causal variant structure

The dataset features two causal variants with indices 194 and 589.

- Causal variant 194 is associated with traits 1, 2, 3, and 4.
- Causal variant 589 is associated with traits 2, 3, and 5.

\
[`library`](https://rdrr.io/r/base/library.html)`(`[`colocboost`](https://github.com/StatFunGen/colocboost)`)`\
`# Loading the Dataset`\
[`data`](https://rdrr.io/r/utils/data.html)`(``Ind_5traits``)`\
`# Run colocboost `\
`res`` ``<-`` `[`colocboost`](https://statfungen.github.io/colocboost/reference/colocboost.md)`(``X ``=`` ``Ind_5traits``$``X``, Y ``=`` ``Ind_5traits``$``Y``)`\
`#> Validating input data.`\
`#> Starting gradient boosting algorithm.`\
`#> Gradient boosting for outcome 4 converged after 40 iterations!`\
`#> Gradient boosting for outcome 5 converged after 59 iterations!`\
`#> Gradient boosting for outcome 1 converged after 61 iterations!`\
`#> Gradient boosting for outcome 3 converged after 91 iterations!`\
`#> Gradient boosting for outcome 2 converged after 94 iterations!`\
`#> Performing inference on colocalization events.`\
`#> Extracting colocalization results with pvalue_cutoff = 0.001, cos_npc_cutoff = 0.2, and npc_outcome_cutoff = 0.2.`\
`#> Keep only CoS with cos_npc >= 0.2. For each CoS, keep the outcomes configurations that pvalue of variants for the outcome < 0.001 and npc_outcome >0.2.`

Colocalization events summary:

\
`cos_summary`` ``<-`` ``res``$``cos_summary`\
`cos_summary``[``, `\
`  `[`c`](https://rdrr.io/r/base/c.html)`(`\
`    ``"colocalized_outcomes"``,`\
`    ``"colocalized_variables"``,`\
`    ``"colocalized_variables_vcp"``,`\
`    ``"cos_npc"``,`\
`    ``"purity"`\
`  ``)`\
`]`\
`#>   colocalized_outcomes          colocalized_variables`\
`#> 1       Y1; Y2; Y3; Y4 rs_186; rs_194; rs_168; rs_205`\
`#> 2           Y2; Y3; Y5                 rs_589; rs_593`\
`#>                                                   colocalized_variables_vcp`\
`#> 1 0.283698494935173; 0.235702194620198; 0.23031425057979; 0.224446848624299`\
`#> 2                                      0.816960678269378; 0.182762268369858`\
`#>   cos_npc    purity`\
`#> 1  0.9989 0.9941612`\
`#> 2  0.9974 0.9761542`

Trait-level evidence:

\
`res``$``cos_details``$``cos_outcomes_npc`\
`` #> $`cos1:y1_y2_y3_y4` ``\
`#>    outcomes_index relative_logLR npc_outcome`\
`#> Y3              3      2.2563485   0.9890312`\
`#> Y1              1      1.8287003   0.9742005`\
`#> Y4              4      0.9499365   0.8504124`\
`#> Y2              2      0.6414479   0.7227667`\
`#> `\
`` #> $`cos2:y2_y3_y5` ``\
`#>    outcomes_index relative_logLR npc_outcome`\
`#> Y5              5       1.994252   0.9814726`\
`#> Y3              3       1.494458   0.9496581`\
`#> Y2              2       1.475390   0.9477011`
