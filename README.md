<img width="303" height="102" alt="neuroesc_long" src="https://github.com/user-attachments/assets/e185a933-e27d-4436-ab1f-52a4fe389e38" />

# Hillscape_analyses
Code used in the analysis and modelling of our Hillscape apparatus, to be used in conjunction with the [available dataset](https://doi.org/10.5281/zenodo.17634454). There is a complete description of the summary dataset contents in that repository and below.

> [!IMPORTANT]
> These functions require the [summary dataset](https://doi.org/10.5281/zenodo.17634454)

# Contents
- [Specifying data 'parts'](https://github.com/Neuroesc/klustest/blob/main/README.md#specifying-data-parts)
- [Usage](https://github.com/Neuroesc/klustest/blob/main/README.md#usage)

# Code overview
The main code which controls the data analysis is `GIT_audit.m`, the file directories at the beginning of this code will need to be modified so as to load the downloaded summary dataset and set the desired output directories.

The code run_fun_for_audit.m is used to either run on individual data files (which will not be possible with the summary dataset) or run analyses on the summary dataset table in clumaa.mat.

Individual codes are used to generate figures, such as PIT_fig_1_v1.m which creates figure 1. These codes are controlled by an optional switch, set these to 1 to run the individual figure codes, or multiple figure codes sequentially.

# Examples







# Summary data file column explanation
After downloading the [summary dataset](https://doi.org/10.5281/zenodo.17634454) you will find a Matlab table called `clumaa.mat`.

Each row of this table gives the data for a single cell in a single session, thus, there will be 3 rows for every cell (i.e. Arena 1, Hillscape, Arena 2).

Each column contains a different data type, described here.

> [!NOTE]
> Some values are calclulated on the fly by the functions that create each figure in the paper, so not all of the values reported will be found here and not all of the values found here were used or are the ones actually used in the figures - check the values are the ones actually used in the figures.

1. electrode - string, the tetrode and cluster reference for this unit, e.g. 't2_c1' means tetrode 2, cluster 1.
2. ele - scalar, the tetrode the cell was recorded on
3. clu - scalar, the cluster corresponding to this cell in the feature space after spike sorting
4. session - scalar, ranges from 1-3, corresponds to Arena 1, Hillscape, Arena 2.
5. session_name - string, the name of the recording files, composed of the recording date and environment name, note that the hillscape is called 'hills' here.
6. session_duration_s - double, duration of the session in seconds
7. spike_times_s - double vector, spike times, in seconds, since the start of the recording that day (i.e. since the start of Arena 1), for all spikes recorded in this session.
8. spike_index - double vector, index into position data, for every spike the position data point nearest in time to the spike, saves on computation
9. n_spikes - scalar, total spikes observed in this session
10. f_rate_hz - scalar, firing rate in Hz observed in this session
11. spike_map - NxM matrix, map of spike locations, binned in x,y
12. rate_map - NxM matrix, map of firing rate, binned in x,y
13. grid_autocorr - NxM matrix, spatial autocorrelation of rate_map
14. mean_waveforms - 1x4 cell array of 32x1 vectors, the mean waveform recorded on each recording channel, in microvolts
15. std_waveforms - 1x4 cell array of 32x1 vectors, the standard deviation waveform recorded on each recording channel, in microvolts
16. wave_width - scalar, the width of the average waveform in milliseconds, on the channel with the hightest peak amplitude in the mean waveform
17. wave_amps - scalar, the peak amplitude of the average waveform in microvolts, on the channel with the hightest peak amplitude in the mean waveform
18. hd_spikemap - Nx1 vector, frequency of spikes, binned in head direction angles
19. hd_ratemap - Nx1 vector, firing rate map, binned in head direction angles
20. hd_rayleigh_r - scalar, Rayleigh vector length calculated using the hd_ratemap
21. hd_pfd_deg - scalar, preferred firing direction in degrees, calculated using the hd_ratemap, the bin angle with the peak firing rate
22. hd_spikemap_half - 1x2 cell array of Nx1 vectors, map of spike locations, binned in head direction angles for first and second half of the session (column 1 and column 2 respectively)
23. hd_ratemap_half - 1x2 cell array of Nx1 vectors, firing rate map, binned in head direction angles for first and second half of the session (column 1 and column 2 respectively)
24. hd_half_correlation - scalar, Pearson correlation coefficient between 1st and second half head direction firing rate map
25. speed_ratemap - Nx1 vector, firing rate map, binned by running speed
26. speed_slope - scalar, slope of the linear best fit line to the speed_ratemap
27. speed_y_intercept - scalar, y-intercept of the linear best fit line to the speed_ratemap
28. speed_score - scalar, speed score, correlation between running speed and firing rate
29. ahv_spikemap - 1xN vector, frequency of spikes, binned by angular head velocity
30. ahv_ratemap - 1xN vector, firing rate map, binned by angular head velocity
31. uci - string, unique cell identifier, composed of the rat number, date, electrode and cluster
32. ratemap_planar - NxM matrix, map of firing rate, binned in x,y, earth horizontal projection
33. planar_amap - NxM matrix, spatial autocorrelation of ratemap_planar
34. planar_field_info - scalar, 6 columns, place field information calculated from the central field of planar_amap (see hill_analysis_v2.m), these values represent the average place field: 
- radius,
- major axis length,
- minor axis length,
- height,
- width,
- orientation
35. planar_spatial_info_shuffles - scalar, 6 columns, spatial information content (SI) calculated on ratemap_planar:
- bootstrapped skaggs SI
- SI z-scored relative to shuffles
- bootstrapped grid score
- grid score z-scored relative to shuffles
- bootstrapped Rayleigh vector length
- Rayleigh vector length z-scored relative to shuffles
36. ratemap_surficial - NxM matrix, map of firing rate, binned in x,y, data projected onto the maze surface, then flattened
37. surficial_amap - same as planar_amap but for ratemap_surficial
38. surficial_field_info - same as planar_field_info but for ratemap_surficial
39. surficial_spatial_info_shuffles - same as planar_spatial_info_shuffles but for ratemap_surficial
40. partn - part number, 1-3, corresponds to Arena 1, Hillscape, Arena 2.
41. cell_type - scalar, 2 columns, the automated and manual cell type, respectively (see get_manual_cell_type_v2.m): 1 = noise, 2 = pyramidal, 3 = place cell, 4 = grid cell, 5 = HD cell, 6 = interneuron, 7 = spatial cell. Automated cell typing was carried out algorithmically, the manual cell type was manually curated.
42. repetition_score - scalar, 2 columns, repetition score (see hill_repeat_analysis_v2.m):
- column 1 is the score calculated using the central row of planar_amap
- column 2 uses the horizontally average of planar_amap.
43. amap_cent - 1xM matrix, the central row of planar_amap
44. planar_fields - scalar, number of place fields detected in ratemap_planar (see hill_field_analysis.m).
45. planar_field_data - Mx9 table, place field data for ratemap_planar (see hill_field_analysis.m), one field per row, output of regionprops, for every field:
- area
- centroid
- weighted centroid
- major axis length
- minoraxis length
- orientation
- maxintensity
- boundingbox
- pixelidxlist
46. surficial_fields - same as planar_fields, but for ratemap_surficial
47. surficial_field_data - same as planar_field_data, but for ratemap_surficial
48. hd_3d_ratemap - 65x65 matrix, azimuth x pitch firing rate map (see PIT_pitch_tuning.m)
49. hd_3d_curves - 2x60 matrix:
- azimuth HD tuning curve on row 1
- pitch HD tuning curve on row 2 (see PIT_pitch_tuning.m).
50. hd_3d_info - scalar, 10 columns, (see PIT_pitch_tuning.m):
- azimuth PFD, pitch PFD
- azimuth spatial information content
- azimuth signal to noise ratio
- azimuth spatial information content z-scored relative to shuffle
- azimuth signal to noise ratio z-scored relative to shuffle
- pitch spatial information content
- pitch signal to noise ratio
- pitch spatial information content z-scored relative to shuffle
- pitch signal to noise ratio z-scored relative to shuffle, 
51. pos_idx - row index into posdata
52. between_session_stability - scalar, 12 columns, correlations between session ratemaps (see PIT_correlations.m):
- arena 1 vs hillscape
- arena 1 vs arena 2
- arena 1 vs hillscape (ridge troughs only)
- arena 1 vs arena 2 (ridge troughs only)
- arena 1 vs hillscape (borders only)
- arena 1 vs arena 2 (borders only)
- arena 1 vs hillscape (center only)
- arena 1 vs arena 2 (center only)
- arena 1 vs hillscape (corners only)
- arena 1 vs arena 2 (corners only)
- arena 1 vs hillscape (aligned to left wall)
- arena 1 vs hillscape (baligned to right wall)
53. ratemap_pitch_half - MxN, two columns, same as ratemap_planar but for all the data where the animal's head was pitched below 0 (horizontal) and all the data where the animal's head was pitched above 0, respectively (see PIT_correlations.m)
54. pitch_stability - scalar, the correlation between the maps in ratemap_pitch_half (see PIT_correlations.m)
55. ratemap_planar_half - MxN, two columns, same as ratemap_planar but for the first and second half of each session, respectively (see PIT_correlations.m)
56. ratemap_surficial_half - MxN, two columns, same as ratemap_surficial but for the first and second half of each session, respectively (see PIT_correlations.m)
57. within_session_stability - scalar, 6 columns, correlations between session ratemaps for session 1st and 2nd halves, uses ratemap_planar_half (see PIT_correlations.m):
- half 1 vs half 2
- empty
- half 1 vs half 2 (troughs only)
- half 1 vs half 2 (borders only)
- half 1 vs half 2 (center only)
- half 1 vs half 2 (corners only)
58. field_anisotropy - unused
59. planar_field_aspect_ratio - scalar, calculated from planar_field_info (see PIT_field_anisotropy.m): field height - width / height + width
60. surficial_field_aspect_ratio - same as planar_field_aspect_ratio but for surficial_field_info
61. planar_field_elongation - unused
62. rat - string, rat name
63. ratn - scalar, rat ID represented as a number
64. surficial_field_elongation - unused

In additon to the above, a behavioural table can also be extracted as:
`posdata = clumaa.Properties.CustomProperties.posdata;`
Each row of this table corresponds to a behavioural session, each column contains specific data:
 
1. pos_idx - position data index, or session number, matches pos_idx in clumaa, so this allows the user to match spike data in clumaa with behaviour/position data in posdata
2. session_times - 3x2 matrix, the start and end time, in seconds, since the start of recording, for each session
3. maze_frame - 36x4, x,y vertices of a wireframe maze fitted to the position data
4. pos - Nx25 table containing the position data:
- session - recording session, 1-3 for Arena 1, hillscape, Arena 2
- pot - time, in seconds, from recording start
- pox - position x coordinates in cm
- poy - position y coordinates in cm
- poh - head direction in degrees
- rx - red LED position x coordinates in cm
- ry- red LED position y coordinates in cm
- gx - green LED position x coordinates in cm
- gy - green LED position y coordinates in cm
- bx - blue LED position x coordinates in cm
- by - blue LED position y coordinates in cm
- poz - position z coordinates in cm
- pox_planar - position x coordinates in cm, earth horizontal
- poy_planar- position y coordinates in cm, earth horizontal
- poz_planar - position z coordinates in cm, earth horizontal
- pox_surficial - position x coordinates in cm, maze projection
- poy_surficial - position y coordinates in cm, maze projection
- poz_surficial - position z coordinates in cm, maze projection
- poz_curve - distance from left wall along surface of maze
- rz - red LED position z coordinates in cm
- gz - green LED position z coordinates in cm
- bz - blue LED position z coordinates in cm
- yaw - 3D head orientation, yaw, degrees
- pitch - 3D head orientation, pitch, degrees
- roll - 3D head orientation, roll, degrees
- anisotropy_map - NxM matrix, 3 columns, columns correspond to arena 1, hillscape, arena 2. Behavioural anisotropy, mapped and binned in x,y (see PIT_behaviour.m)
5. anisotropy_score - scalar, 3 columns, columns correspond to arena 1, hillscape, arena 2. Overall anisotropy score: (time_y - time_x) ./ (time_y + time_x)
6. rat - string, rat name
7. ratn - scalar, rat ID represented as a number
8. hd_3d_dwell - NxM matrix, 3 columns, columns correspond to arena 1, hillscape, arena 2. Dwell time maps corresponding to hd_3d_ratemap in clumaa.
9. mov_3d_dwell - NxM matrix, 3 columns, columns correspond to arena 1, hillscape, arena 2. Azimuth x pitch dwell time maps of movement direction rather than head direction.

