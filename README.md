# What is this?
This code is for the cross-section measurement of single-top quark processes. This code is for the final state exactly 1 muon and at least 2 jets. 


## The Analysis Codes folder contains the following files-
SingleTopAna.C -> The main C script <br>
SingleTopAna.h -> The associated header file <br>
ana.C -> Driver script <be>


**What does it do? How do I run?**

These SingleTopAna files are the core of the entire framework. The templates for these files are generated by running MakeSelector() for the Events tree from the nanoAOD files. The header file contains the class named AnaScript, and the C file utilizes its features event by event. The C file is a driver script. It decides which code to run over which sample, the names of output files, its location and so on.

**How do I run this on root?** <br>
To compile the code
> [] .L SingleTopAna.C+

To run the code over the various samples
> [] .x ana.C()
