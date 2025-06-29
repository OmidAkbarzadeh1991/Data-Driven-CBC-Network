Our simulations were executed in **MATLAB R2023b** on a **MacBook M2** with **32 GB RAM**.  
The values reported in **Table 1** and **Table 2** refer only to the scripts for **Algorithm 1**  
(`xxxxx_Algorithm_1`). In contrast, the network-simulation scripts (`xxxxx_Network_Simulation`)  
may take **more than 30 minutes** to finish, depending on the number of subsystems and the  
available computational resources.  

> **Note** The reported **running times (RT)** and **memory usage (MU)** correspond to the time  
> needed to solve the SOS problem for each subsystem, as detailed in **Algorithm 1**.

## Prerequisites & Installation

### Supported operating systems
Any relatively recent version of **MATLAB** on any OS should work, but for consistent results we recommend the version used in our tests: **MATLAB R2023b**.

### Required optimization toolboxes & solvers  
Add the following packages to your MATLAB path:

- **YALMIP** R20210331  
- **SOSTOOLS** 4.03  
- **MOSEK** 10.1.29  

### Setup steps
1. Download  
   - **YALMIP**    <https://yalmip.github.io/download/>  
   - **MOSEK**    <https://mosek.com/downloads/>  
   - **SOSTOOLS** <https://github.com/aalto-ics/ssv/tree/master/sostools>  
2. Extract each `.zip` archive.  
3. In MATLAB, open **Set Path… → Add with Subfolders…** and add the extracted folders.  
4. Save the updated path and restart MATLAB (optional but recommended).

> **Note** All three packages **must** be on the MATLAB path before running any of the provided scripts.


