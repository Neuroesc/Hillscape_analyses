<img width="303" height="102" alt="neuroesc_long" src="https://github.com/user-attachments/assets/e185a933-e27d-4436-ab1f-52a4fe389e38" />

# Hillscape_analyses
Code used in the analysis and modelling of our Hillscape apparatus, to be used in conjunction with the [available dataset](https://doi.org/10.5281/zenodo.17634454). There is a complete description of the summary dataset contents in that repository.

> [!IMPORTANT]
> These functions require the [summary dataset](https://doi.org/10.5281/zenodo.17634454)

The main code which controls the data analysis is `GIT_audit.m`, the file directories at the beginning of this code will need to be modified so as to load the downloaded summary dataset and set the desired output directories.

The code run_fun_for_audit.m is used to either run on individual data files (which will not be possible with the summary dataset) or run analyses on the summary dataset table in clumaa.mat.

Individual codes are used to generate figures, such as PIT_fig_1_v1.m which creates figure 1. These codes are controlled by an optional switch, set these to 1 to run the individual figure codes, or multiple figure codes sequentially.
