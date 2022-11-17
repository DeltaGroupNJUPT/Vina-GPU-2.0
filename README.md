# Vina-GPU 2.0
Vina-GPU 2.0 accelerates AutoDock Vina and its related commonly derived docking methods, such as QuickVina 2 and QuickVina-W with GPUs.
Vina-GPU 2.0 includes the following three docking methods,including Vina-GPU+,QuickVina 2-GPU and QuickVina-W-GPU.

## Vina-GPU+
Vina-GPU+ further accelerates Vina-GPU and facilitates single receptor-multi-ligand docking.
The compile and run process is described in [here](https://github.com/DeltaGroupNJUPT/Vina-GPU-2.0/tree/main/Vina-GPU%2B).

If you would like to develop and build Vina-GPU+ from source please look at README.md(https://github.com/DeltaGroupNJUPT/Vina-GPU-2.0/blob/main/Vina-GPU%2B/README.md).

## QuickVina 2-GPU
The compile and run process is described in [here](https://github.com/DeltaGroupNJUPT/QuickVina2-GPU).

## QuickVina-W-GPU
The compile and run process is described in [here](https://github.com/DeltaGroupNJUPT/QVina-W-GPU).

## Graphic User Interface (GUI)
A graphic user interface (GUI) is provided for users on **Windows** OS
1. first make sure that  `Vina-GPU2.0.exe` can run on a terminal
2. put the `Vina-GPU2.0.exe` and all `.bin` formatted files in `./Vina-GPU-2.0/GUI/exec` and overwrite the original files
3. run the `Vina-GPU2.0.exe`file within  `./Vina-GPU-2.0/GUI` to start up the Vina-GPU 2.0 GUI
4. select docking methods
5. select the input and output files
6. set the box center, the box size, thread and search_depth
7. click the `start` button to run Vina-GPU 2.o
8. click `score` button to get docking scores

## Citation
* Tang, Shidi et al. “Accelerating AutoDock Vina with GPUs.” Molecules (Basel, Switzerland) vol. 27,9 3041. 9 May. 2022, doi:10.3390/molecules27093041
* Trott, Oleg, and Arthur J. Olson. "AutoDock Vina: improving the speed and accuracy of docking with a new scoring function, efficient optimization, and multithreading." Journal of computational chemistry 31.2 (2010): 455-461.
* Hassan, N. M. , et al. "Protein-Ligand Blind Docking Using QuickVina-W With Inter-Process Spatio-Temporal Integration." Scientific Reports 7.1(2017):15451.
* Amr Alhossary, Stephanus Daniel Handoko, Yuguang Mu, and Chee-Keong Kwoh. "Fast, accurate, and reliable molecular docking with QuickVina 2. " Bioinformatics (2015): 2214–2216.
