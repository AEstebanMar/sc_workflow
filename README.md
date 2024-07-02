Total refactorisation of @SergioAlias's sc-cellranger-workflow, sc-neuron-organoids and sc-retina repos. Heavy optimisation, able to exploit supercomputing resources. The three repos now work inside the same workflow. Easy to use, dynamically coordinates integration according to provided experiment design. Will be integrated in the ExpHunterSuite R package once full workflow is fully functional.

Workflow can analyse an imported counts matrix (check config_default for where to specify it). If that is the case, DO NOT run daemon 1, start at 2.
