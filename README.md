# PDS analysis class

<!-- markdown-toc start - Don't edit this section. Run M-x markdown-toc-refresh-toc -->
**Table of Contents**

- [PDS analysis class](#pds-analysis-class)
    - [Introduction](#introduction)
    - [How to use:](#how-to-use)
        - [Most common scripts](#most-common-scripts)
<!-- markdown-toc end -->

## Introduction

Here are my analysis class for the Photon Detectors data.


## How to use:

After cloning the repository, `cd` inside `./codes/` and run:

``` bash
source ./setup_codes.sh
```

This will create `templates` scripts to help you in your analysis. 

After doing so, in the file `load_my_class.sh` you will find instructions on how to setup an automatic loading of the class for "online" analysis. 

It will give you the instructions (**with the correct path**) as following:

``` example
# insert this line to .bashrc
# alias myclass="source ~/{YOUR_PATH}/load_my_class.sh"
# to use the sample (analyzer) class, go to any folder which have "analyzed.root" file and run:
# myclass
# If you want to execute a specific file instead, you can run
# myclass {path_to_file}/your_file.root 
# If you want to change old data, you need to give the number of points per waveform as following:
# myclass {file.root} {# of pts}
# If you want just to load the class without any data, you can run
# myclass no
```


If everything works you can run `myclass` and it will load it.

To check a few methods, you can run:
``` root
s.sample_plot(4,"",16);
```

This will draw event `4` with TV1D denoise filter parameter equal `16`.<br>
If you have two or more channels, you can draw the event (`4`) in the same graph doing:

``` root
s.kch = 1
s.sample_plot(4,"SAME",16)
```

`kch` refer to the logical numbering of the channels. If you want to set it by name instead you can use, for example, `s.setChannel("Ch4.")`.


For this "online" plottings, the methods are located in `analyzer_plots.C` and `analyzer.C`.


This `ANALYZER` class gives a fast and simple way to write scripts and handle the data. 


### Most common scripts

The most common scripts used are:

#### TODO: update

``` example
adc_read_all_data.C # This process data, subtract baseline and save as .root
giveMeSphe.C        # script to get single photo-electron spectrum
fit_sphe_wave1.C    # script to fit the spectrum
```

To use the scripts, copy the templates to the folder you have your data (or to any place different than their original code). Set the parameters as needed.<br>
Check the parameters carefully, if not clear, check their implementation in the codes `readingCodes.C` and `calibrationCodes.C`.







