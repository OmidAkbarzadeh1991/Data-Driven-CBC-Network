Our simulations were executed in **MATLAB R2023b** on a **MacBook M2** with **32 GB RAM**.  
The values reported in **Table 1** and **Table 2** refer only to the scripts for **Algorithm 1**  
(`xxxxx_Algorithm_1`). In contrast, the network simulation scripts (`xxxxx_Network_Simulation`)  
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
   - **YALMIP**    <https://github.com/yalmip/YALMIP/archive/refs/tags/R20210331.zip>  
   - **MOSEK**    <https://www.mosek.com/downloads/10.1.29/>  
   - **SOSTOOLS** <https://github.com/oxfordcontrol/SOSTOOLS/archive/refs/tags/v4.03.zip>  
2. Extract each `.zip` archive.  
3. In MATLAB, open **Set Path… → Add with Subfolders…** and add the extracted folders.  
4. Save the updated path and restart MATLAB (optional but recommended).

> **Note** All three packages **must** be on the MATLAB path before running any of the provided scripts.
> 
> **Remark&nbsp;**  
> If you are using a fresh installation of **MATLAB**, you must also install the **Symbolic Math Toolbox**.

### How to install the Symbolic Math Toolbox

1. **Open MATLAB**.  
2. Navigate to the **Home** tab.  
3. Select **Add-Ons → Get Add-Ons**.  
4. In the **Add-On Explorer**, search for **“Symbolic Math Toolbox”**.  
5. Click the toolbox entry and choose **Install**.  
6. Follow the on-screen prompts; MATLAB will handle the download and installation.  
7. After installation, run `ver` in the Command Window to verify that the toolbox is listed.

### MATLAB Path Instructions

Follow the steps below **for each toolbox** (e.g., **YALMIP R20210331** and **SOSTOOLS 4.03**):

1. Open **MATLAB**.  
2. Go to the **Home** tab in the MATLAB toolbar.  
3. Click **Set Path** in the **Environment** section.  
4. In the dialog that appears, select **Add with Subfolders…**.  
5. Locate and select the extracted toolbox folder (after unzipping, if applicable).  
6. Click **Save** to confirm the changes.

### Alternative method

1. Place the toolbox folder in the same location as your other `.m` files.  
2. In MATLAB’s **Current Folder** browser, right-click each folder and choose  
   **Add to Path → Selected Folders and Subfolders**.

### MOSEK 10.1.29 setup

1. **Get a license**  
   - Request a free academic license → [MOSEK Academic Licenses](https://www.mosek.com/products/academic-licenses/)  
   - Or request a time-limited trial license → [MOSEK Trial Licenses](https://www.mosek.com/products/trial/)  
   - Watch for an email (check your spam folder) containing the `mosek.lic` file.

2. **Save the license file**  
   - Download the attached `mosek.lic` from the email.

3. **Place `mosek.lic` in the correct folder**

   | OS       | Destination path (create the folder if it doesn’t exist)             |
   |----------|---------------------------------------------------------------------|
   | Linux    | `/home/<YOUR_USER_NAME>/mosek/mosek.lic`                            |
   | macOS    | `/Users/<YOUR_USER_NAME>/mosek/mosek.lic`                           |
   | Windows  | `C:\Users\<YOUR_USER_NAME>\mosek\mosek.lic`                         |

4. **Install MOSEK binaries**  
   Follow the OS-specific installation guide in the MOSEK documentation: <https://docs.mosek.com/10.1/install/installation.html>.


