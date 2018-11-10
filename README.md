# Dependencies
* Fieldtrip (21092018 version)
* spm12
* EEGLAB (eeglab14_1_2b)
* ICBM152 (T1, 0.5mm)


# Analysis
## 1-preprocessing
**Remarks:** The scripts has three "big" sections: computing the headmodel, electrode allignment with the template, and the actual preprocessing and decomposition of the data. Here are all the steps:

1. Setup path, dependencies, and important variables
    * Use spm12 because Fieldtrip's default version is not compatible with my MATLAB version
2. Create headmodel based on ICBM 152 (unbiased standard magnetic resonance imaging template brain volume for normal population; the most commonly used reference brain for neuroimaging)
    * Get MRI and add fiducials to it
    * Segment the volume (brain, skull, scalp)
    * Prepare mesh (triangulated surface mesh for the volume conduction model)
    * Create headmodel (Boundary Element Method)
3. Electrode allignment with template 
    * Two passes: first allign template fiducials with electrode fiducials then interactively make the fit better
    * Two electrodes of one subject were digitized incorrectly and corrected by taking those same coordinates from a glass head (g64cap.sfp)
4. Load data
5. High pass filter @ 0.5 Hz (Hamming FIR)
6. Load channel location file
7. Trim data from third event to end of data (only analyze data when participants were actually playing)
    * This quantities were obtained from the MIDI files and Reaper sessions
8. Take out line noise using spectral regression
9. Run clean_rawdata
    * Channel rejection: channels with flat signal for over n seconds and (2) channels that are poorly correlated (r < 0.75) with adjacent channels
    * Run ASR, which removes short-time high-amplitude artifacts in the cont data 
10. Spherical interpolation
11. Common Average Reference
12. Output EEG lab file and save it
13. Process baseline using same steps
14. Prepare the EEG file
15. EEG to Fieldtrip structure (min duration for covariance matrix and one with full duration trials)
16. Prepare the leadfield model (forward model for multiple dipole locations)
17. Compute cortical patches basis functions (forward model)
    * We use a coarse version of the AAL atlas (19 regions)
    * Get the leadfield for each patch (n_channels x (3*points in the patch))
    * We assume Gamma = I, white noise
    * Do SVD of leadfield in that patch, and take n singular values that will give you a gamma representation accuracy of 0.85 (this respresentes a trade off between representation error and resolution)
18. LCMV criterion spatial filters for each patch
    * Get basis for each patch
    * take the smallest eigenvalue of Yk to maximize power output (because we do not know the moment of the dipole)
    * compute each patch spatial filter weights
        * covariance matrix:
            * only take out electrodes interpolated in *all* trials
            * trial level covariance matrix (cut trials to min duration of the trials)
            * average cov matrices together
        * minimize the output power subject to a unit response constraint to patch location 
19. Get patches time series
20. Plot average power of each patch on MRI scan
21. Plot patch resolution on anatomical image
    * Basic idea: ideally, this would be 1, but this is impossible to achieve on practice because we must balance representation error against resolution  
    

## 2-STE_graham
**Remarks** This folder includes all the necessary files and scripts to run the STE analysis in a Compute Canada cluster (in this case, Graham). For details of the bundle Serial farm implementation, you can check out this [video](https://www.youtube.com/watch?v=49dC4bmBCic)

1. Copy all the contents of this file to your Graham home directory
    * You can run this command
    
    ```
    scp -r /file/in/local/computer/ $USERNAME@graham.computecanada.ca:/home/username/project_directory
    ```
    
    * P0* directories include three files per subject: the patch time series (19 sources), the patch labels, and the trial labels
2. Once in the server, run bundle_up.sh
    * This is a bash script that basically bundles up 3 "cases" together to be run together in the serial farm (this number is controlled by N_bundle)
    * You will need a table.dat file. This file is created by the Python script "generate_graham_list.py". Each line of this file is the command that represents a "case" (or an individual job) in the serial farm. 
    * As you can see, each job loads Python, then the SciPy Stack, then runs 2-analysis_ste.py with certain parameters
3. The script 2-analysis_ste.py intakes 7 parameters: Pair, Source frequency, target frequency, subject a, subject b, delay. I wrote it this way to make it easier to create the serial farm on Graham. It does as follows:
    1. It get the parameters from the terminal input
    2. It gets the data from the .mat files (you will need to change the directory here!)
    3. Do the baseline time-frequency decomposition to get the normalizing factor
    4. Do trial time-frequency decomposition
    5. Do computation of symbolic transfer entropy and output a csv
4. Once it's done, copy the csv directory to your computer (run this bash code on your computer)   

    ```
    scp -r $USERNAME@graham.computecanada.ca:/home/username/project_directory/csv_grand_matrices /file/in/local/computer/
    ```