# Journal club discussion guide

**Antão, L. H., Weigel, B., et al. (2022). Climate change reshuffles northern species within their niches. Nature Climate Change 12, 587–592.**

DOI: https://doi.org/10.1038/s41558-022-01381-x

Prepared from the local paper, including its Methods and Extended Data. The separate Supplementary Information and analysis code have not been checked. References below distinguish the printed article page from the PDF page when necessary. Suggested timing is flexible: approximately 40 minutes of discussion.

## The central question

**Can substantial ecological change occur even when a region retains broadly the same species, because the climatic conditions those species experience change?**

The paper's contribution is to look beyond changes in regional species lists. It examines whether associations between occurrence and climate change across decades, regions and taxonomic groups. The authors interpret changing response shapes as species occupying different parts of their climatic niches.

The main question for the group is how convincingly the analysis connects those changing associations to movement within niches.

## A short opening you can read aloud

“This paper asks what climate change does to communities when the regional species list changes relatively little. The authors analyse four decades of Finnish monitoring data across seven groups, fitting separate community models for different regions and decades. They find that climate–occurrence associations change, often more strongly in the north, and that the relative importance of different climatic variables also changes. Their interpretation is that species are being reshuffled within their niches. I suggest we discuss both the value of that perspective and whether the models distinguish movement within a stable niche from changes in the associations themselves.”

## Study design in two minutes

| Component | What the authors did |
|---|---|
| Scope | 1,478 species at 6,504 sites across approximately 1,200 km of Finland, spanning 1978–2017 overall. |
| Groups | Birds, mammals, small rodents, butterflies, moths, forest understory vascular plants and freshwater phytoplankton. |
| Response | Converted abundance data to occurrence; fitted presence/absence models. |
| Subsets | Three boreal zones and four decade bins; coverage differs among groups. Plants have three survey periods, with one survey year per relevant decade; butterflies cover only the last two decades. |
| Models | 63 separate HMSC joint species distribution models, each for one group × zone × decade. Probit link; site and year random effects, with spatially structured latent variables. Plants omit year and NAO. |
| Climate | Annual mean temperature, annual precipitation and snow-cover duration, each with quadratic terms; annual North Atlantic Oscillation (NAO) index with a linear term. Local climate data have 10 × 10 km resolution. |
| Inclusion | At least ten occurrences per species within a subset; subsets require at least 20 samples and six species. Some large subsets were randomly reduced to 1,000 records. |
| Main outputs | Supported response-shape categories; variance partitioning scaled by discriminatory power; comparisons of species lists and spatial community dissimilarity. |

Source: Methods, PDF pages 7–9. The full 40-year span does not apply equally to every group.

## Walk through the figures

### 1. Start with Figure 2: what would “within a niche” mean? (5 minutes)

**Find it:** printed p. 589, PDF p. 3.

On a unimodal climate–occurrence curve, a sampled climate interval might cover the increasing side, the peak, or the decreasing side. As the available climate changes, the sampled interval can occupy a different part of that curve.

The authors label these increasing, bell-shaped and decreasing responses as the lower end, optimum and upper end of a niche. “Lower” means lower values of the climatic predictor; it does not mean lower latitude, poorer habitat or lower abundance.

**Ask:** What assumptions allow a fitted occurrence curve to tell us where a species sits within its niche?

**Facilitator notes:** The schematic encourages a stable underlying niche interpretation. The analysis estimates separate curves within Finland, not complete physiological niches across each species' range. Realized occurrence also reflects habitat, dispersal, interactions and observation. A fitted occurrence peak is therefore not automatically a fitness optimum. The authors explicitly acknowledge that their geographic coverage captures only part of many species' ranges.

### 2. Figure 1: how comparable are the data? (4 minutes)

**Find it:** printed p. 588, PDF p. 2.

Panel a shows unequal time coverage; panel b shows the sampling geography; panel c shows changing climate. Temperature and precipitation increased more strongly in the northernmost zone, while snow-cover declines were strong in both northern and southern zones.

**Ask:** Which comparisons are most persuasive: within a group over time, between zones, or between groups?

**Facilitator notes:** Breadth across taxa is a major strength, but sampling designs and time coverage differ. The authors fit groups separately and examine sampling distributions. That helps, while leaving comparability and detection worth discussing. For example, mammals are recorded through snow tracks, making observation conditions particularly relevant to a snow-related analysis.

### 3. Figure 3: what actually changes? (7 minutes)

**Find it:** printed p. 590, PDF p. 4; Extended Data Fig. 2, PDF p. 12.

Figure 3a pools species across groups for each climate variable. Figure 3b separates the temperature results by group. Yellow indicates increasing, teal bell-shaped, and pink decreasing responses. The bars show proportions of all species in the relevant models with supported responses, using a 90% posterior-support threshold. The unfilled fraction is not classified as supported; it does not establish ecological insensitivity. Grey circles indicate missing models.

Useful patterns to identify together:

- Temperature responses are more prevalent in the north, and the northern pooled proportion increases across decades.
- More northern bird and mammal species show positive temperature associations in later decades.
- Precipitation responses generally become less prevalent; NAO responses become more prevalent, often negative.
- Group-specific patterns differ substantially. Southern butterflies increasingly show bell-shaped or decreasing temperature responses.

**Ask:** Do changing category proportions demonstrate that the same individual species changed niche position?

**Facilitator notes:** These panels summarize category frequencies, not paired species transitions. Species eligibility and group representation vary between subsets. The low regional species-list turnover supports continuity, but does not itself establish which species changed category or why. A particularly informative follow-up would track shared species and their posterior response curves through time on a common climate axis.

### 4. Figure 4: what does increasing “importance” mean? (7 minutes)

**Find it:** printed p. 591, PDF p. 5.

The authors scale variance components by Tjur's R² and regress the resulting species-level values on decade. Tjur's R² is the mean predicted occurrence probability at observed presences minus that at observed absences. Scaling is intended to account for differences in model discrimination.

Red indicates an increasing contribution over time and blue a decreasing contribution. **These colours do not indicate positive and negative climate effects on occurrence.** The colour scales also differ among groups.

In the worked northern-bird example, temperature's weighted contribution increases, snow's decreases, and precipitation and NAO have no significant linear trend. Across groups, plants and mammals tend towards increasing contributions, while moths and phytoplankton tend towards decreasing contributions.

**Ask:** Could a predictor's weighted contribution change without a change in the species' underlying climate response?

**Facilitator notes:** Yes, potentially: the sampled climate range, predictor covariance, species prevalence, model discrimination and competing variance components can change. This makes the metric useful descriptively, but it is not a direct measure of causal climate influence. Ask whether the second-stage regressions account for repeated species and propagate uncertainty from the fitted models. The Methods describe ordinary linear regression; checking the code would be needed to establish all implementation details.

### 5. Extended Data Figure 3: distinguish two kinds of turnover (3 minutes)

**Find it:** PDF p. 13.

This figure measures spatial dissimilarity among sites within each zone and decade using the multiple-site Simpson index. The authors standardize the number of sampled sites within each taxonomic group and repeat sampling 1,000 times. Error bars describe the resulting distributions.

**Ask:** Does low change in a regional species list mean that local communities remain unchanged?

**Facilitator notes:** No. Regional persistence can coexist with considerable rearrangement among sites. This figure concerns spatial differentiation and homogenization; it is not a direct estimate of temporal colonization and extinction at repeatedly sampled sites. The species-list comparisons appear separately in the supplementary analyses. Also note the different vertical scales among panels.

## Three deeper questions (choose one or two)

### A. Have the authors identified movement within stable niches?

**Prompt:** “What evidence would distinguish a fixed niche exposed to new climate from a genuinely changing climate–occurrence relationship?”

**Why it matters:** Separate decade models allow associations to vary, but do not directly compare a shared, stable response curve against time-varying alternatives. Changed curve categories could reflect shifting climate exposure, changed realized relationships, sampling differences or uncertainty.

**Useful extension:** Fit a model with shared climate responses across decades and compare it with one allowing temporal changes. Evaluate predictions in later decades and examine curves over overlapping climate ranges. Simulations under fixed niches and observed sampling/climate changes would test whether the analysis recovers the proposed mechanism.

### B. How strong is the attribution to climate change?

**Prompt:** “What can site and year random effects absorb, and what can they not identify?”

**Why it matters:** Spatial latent variables account for residual spatial structure; they do not guarantee that climate estimates are free from confounding by land use or other omitted factors. They also do not constitute an explicit detection model. The authors acknowledge land-use change as a possible contributor.

**Useful extension:** Add measured habitat and land-use changes, examine consistent-site subsets and survey effort, and separate within-site temporal climate variation from between-site climatic differences. This would help distinguish a temporal response to warming from a spatial climate association.

### C. How much ecological meaning rests on classification?

**Prompt:** “Should we put species into response categories, or retain continuous curves and their uncertainty?”

**Why it matters:** A support threshold can convert a gradual change in precision into a change in category frequency. The ten-occurrence filter also changes which species enter each model. Failure to meet the support threshold does not demonstrate no response.

**Author safeguard:** The authors report sensitivity analyses for the gradient cutoffs used to define response shapes (Supplementary Fig. S6). That is valuable, but differs from testing the posterior-support threshold, sampling balance or changing species membership. The separate supplement has not been inspected for this guide.

**Useful extension:** Repeat comparisons for a shared-species subset, vary support thresholds, and propagate posterior uncertainty into category proportions and temporal contrasts.

## What the paper does well

- Uses unusually broad, long-term monitoring across terrestrial and freshwater groups.
- Treats climate as multidimensional, including snow cover and NAO alongside temperature and precipitation.
- Makes spatial and temporal context central instead of assuming a single relationship everywhere and always.
- Combines modelled climate associations with independent descriptive community metrics.
- Reports sampling checks, response-classification sensitivity analyses, and publicly available data/code links.

## A balanced closing assessment

The strongest supported message is that climate–occurrence associations and their relative explanatory contributions differ across time, latitude and taxa, while broad regional species pools can remain comparatively similar. This is a useful reason to look beyond richness and regional species turnover.

Movement within climatic niches is a plausible interpretation, but the analysis does not uniquely identify movement along fixed niches or measure physiological optima. Consequences for fitness, interactions and ecosystem integrity remain ecological implications rather than directly tested outcomes.

**Final question:** “If we could add one analysis to strengthen the title's claim, which would we choose, and what result would change our interpretation?”

## Connection to an atlas/HMSC discussion

Without assuming the local project's implementation matches this paper, three useful transferable questions are:

1. When comparing models across atlas periods, are we comparing equivalent species sets, sampling support and environmental ranges?
2. Does a change in variance partitioning describe changing predictive contribution, or are we interpreting it as changing ecological sensitivity?
3. Would predictions evaluated over a common environmental range, together with tests of transfer between periods, strengthen the ecological interpretation?

## Quick answers to likely questions

**Did they model abundance?** No. Abundance observations were converted to occurrence for these models.

**Did they estimate fundamental niches?** No. They estimated climate–occurrence relationships within the sampled Finnish environments.

**Does a positive temperature association mean warming benefits a species?** It is consistent with that interpretation, but does not establish a temporal fitness benefit or population increase.

**Does low regional turnover mean little biodiversity change?** No. Local redistribution, abundance change and changing climate associations may occur without regional gains or losses.

**Was predictive transfer across decades demonstrated?** The Methods describe model-fit assessment using AUC and Tjur's R². They do not describe a test training in one decade and predicting another; do not present the reported fit as demonstrated temporal transferability.

**What should participants read first?** The abstract and Figure 2 for the idea, Figures 3–4 for the evidence, and the Methods sections on joint modelling and scoring niche position for the assumptions.
