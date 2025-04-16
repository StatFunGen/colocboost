
## Data

  - `data.R`: Contains example datasets for testing and demonstration.

## Source code structure for developer

Implementation of ColocBoost algorithm falls roughly in the structure of ColocBoost paper. 
That is, we introduce a multi-task regression problem for L traits, followed by the dynamic coupling strategy with SEC learner, 
and post assemble and inference. Implementation-wise,


1. `colocboost.R` implements the main interface function that users interact with directly.

2. `colocboost_workhorse.R`: The core interface of dynamic coupling strategy with SEC learner.
    - `colocboost_check_update_jk.R`: The strategy to determine best update variant for the subset of traits.
    - `colocboost_update.R`: The single effect learner/coupler (SEC) for the best update variant and traits.
    - `colocboost_one_causal.R`: The special case of ColocBoost with per-trait-per-causal assumption with/without LD information.
  
3. `colocboost_assemble.R` implements the core interface of post assemble and inference SEC learners from 2.
    - `colocboost_assemble_cos.R`: The function to create 95% CoS of different colocalization events.
    - `colocboost_assemble_ucos.R`: The function to create 95% CS of trait-specific effects.
    - `colocboost_inference.R`: Post inference functions includes modularity hierarchical clustering method, remove spurious signals, definitation of colocalization evidence, et al.
    - `colocboost_utils.R`: Utility functions includes refining colocalization confidence sets from different SEC and other utilities.
    - `colocboost_output.R`: Format and export analysis results
  
4. `colocboost_plot.R` implements various visualization options for visualize colocboost results.

