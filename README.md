## License  
This file is part of the project hyperMusic. All of hyperMusic code is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. megFingerprinting is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with megFingerprinting. If not, see <https://www.gnu.org/licenses/>.

## Status
On going - Finishing up data acquisition and data analysis 

## Objective
To characterize a mapping of the neural substrates of social interaction using a music paradigm

## Supervisors
[Dr Laurel Trainor](https://trainorlab.mcmaster.ca/people/ljt) at the [McMaster Institute for Music & the Mind](https://mimm.mcmaster.ca/)

## Contents
### _dependencies
* Includes a .mat file that pairs electrode numbers to electrode names (the gTec system records the electrode names as numbers)
* Includes the preprocessing function, hM_preProEEGLAB

### _figs
* Done with Python's seaborn. This folder includes...
    1. Preliminary results of pre and post self reports of how much people like each other (people seem to like each other a bit more after playing music with each other)
    2. High pass filter (frequency and phase response) used for preprocessing steps
    3. Data before preprocessing
    4. Data before preprocessing
    5. Data after preprocessing
    6. Hyperbrain networks from 10 polyphonic trials: with current analysis parameters we see no network between people, only within connections
    7. Hyperbrain networks from 10 homophonic trials: with current analysis parameters we see no network between people, only within connections
    8. Symbolic transfer entropy and self reports at the trial level plotted against each other (Synergy, or how much people thought they were playing with each other)
    9. Symbolic transfer entropy and self reports at the trial level plotted against each other (Syncrhonization, or how people thought they were syncronized during that trial)
    10. 8. Symbolic transfer entropy and self reports at the trial level plotted against each other (Quality, or how good people thought that specific trial was)

### 1-self_reports
* Takes in a CSV report data (from the MAQ) questionnaire and compares pre and post scores of how much they like each other
    * Please note that CSV format should be...
        * Pre_Participant A
        * Post_Participant B
        * Pre_Participant A
        * Post_Participant B 
* Output: Figure
* To do: Stats, get more people

### (old)2-preprocessing_input2decomposition
*  gTec EEG files are saved as .hdf5 files. HDF stands for "Hierarchical Data Format". Every object in an an HDF5 file has a name, and they are arranged in a POSIX-style hierarchy. The folders in the system are called groups. You get access to those with the method h5py.File(). 
*  To access subgroups of data, you can use Python dictionary-style interface using the item-retrieval syntax.
*  Use printname and the visit method to iterate through all the groups folders
* __Important__: The gTec file has the next structure 
    * eegData.hdf5 
    * AsynchronData 
        * AsynchronSignalTypes: General info on the trigger signal used 
        * Time: What samples have them triggers 
        * TypeID: random id assigned every time you start recording. It is linked to the actual trigger value in a different xml field.
        * Value: The values of the square wave 
    * RawData 
        * AcquisitionTaskDescription: it has information of each of the channels properties, such as filtering, if it is bipolar, its number/name... 
        * DAQDeviceCapabilities: very similar information to the last container 
        * DAQDeviceDescription: description of the amplifier information 
        * Samples: actual raw data stored as a 2D matrix (voltage, channels) [either V, mV, or uV]
        * SessionDescription: literally, you can add this metadata when recording 
        * SubjectDescription: literally, add the subject's metadata (name, birthdate, description...) 
    * SavedFeatues 
        * NumberOfFeatures: features saved in the metadata of the recording 
    *  Version # Version: gRecorder's version 
    
    
### 2-preprocessing_input2decomposition.m
* Ended up deciding to go for MATLAB for the pre-processing
* Uses function called hM_preProEEGLAB.m found in dependencies, based on EEGLAB
    1. Load file and delete extra events
    2. High pass filter 1 Hz (Hamming windowed sinc FIR)
    3. Load channel locations, delete extra channels, and assign 10 20 name
    4. Take out line noise using spectral regression
    5. Run artifact subspace reconstruction
    6. Get info regarding interpolated channels
    7. Interpolate bad channels using spherical method
    8. Perform CAR 
    9. Trim data around trigger event
    10. Save output .set file
    11. Write bad channels to txt file 
* Iterates through participants folders and outputs preprocessed files
    

### 3-analysis_NSTE_EEG.py
* Includes the necessary functions to calculate NSTE 
* Outputs STE, NSTE, tau's used, delay's used...
    * The output of the PCA (% variance explained)
    * Matrices (csv format)
    * Matrices (brainstorm file (.mat))

### 4-preproc_movement_NSTE.py
* Includes the processing movement pipeline. The movement data was obtained from the videos ugins this implementation of [Flow Analysis](https://www.cefala.org/FlowAnalyzer/). For future iterations, we will use our own implementation

### 5-movement_ste_stats.r
* Includes the states to be done to the movement data (analysis of variance)
