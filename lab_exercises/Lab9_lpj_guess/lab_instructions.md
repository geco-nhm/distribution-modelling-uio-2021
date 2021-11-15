# BIOS5211 / BIOS9211 - Dynamic Global Vegetation Models
**Lecturers**: Lasse Torben Keetz, Eva Lieungh
:::success
**Contents**
1. <a href="#1-Lecture">Lecture supporting materials</a>
2. <a href="#2-Lab---LPJ-GUESS-Education">LPJ-GUESS Education - model setup</a>
3. <a href="#3-Exercises">Exercises</a>
:::

### **Ice-breaker question :icecream: :**
### <span style="color:red">What is your favorite plant?</span>

- ADD HERE


</br></br>


# 1 Lecture

<a href="https://docs.google.com/presentation/d/1aDKtcyLN4p2jC1Odej0-4CPkVTCdXC7oFSvpeWPSk8Q/edit?usp=sharing">Click here</a> to view the slides!

## Questions and feedback
:::info
***Template***
- Please write your questions in this format, i.e. preceded by a dash
    - Use indented dashes for follow ups / answers
:::

- ***First question here***

</br></br>


<hr>

# 2 Lab - LPJ-GUESS Education
:::success
**Quick navigation**
1. <a href="#21-Installation">Installation</a>
2. <a href="#22-Simulation-checklist">Simulation checklist</a>
3. <a href="#23-Generating-climate-forcing-data">Generating climate forcing data</a>
4. <a href="#24-Defining-the-LPJ-GUESS-parameters">Defining LPJ-GUESS parameters</a>
5. <a href="#25-Starting-the-simulation">Starting the simulation</a>
6. <a href="#26-Monitoring-the-simulation">Monitoring the simulation</a>
7. <a href="#27-Analyzing-the-output">Analyzing the output</a>
:::
### Preparations before the lab
Before we meet on zoom, please follow these steps to ensure we can have a smooth start! You will also find details about the installation procedure in the next section.

:::warning
:rotating_light: **Attention #1**
LPJ-GUESS Education can only be installed on MS Windows. Please contact the teachers if you do not have access to a Windows OS.
:::
:::warning
:rotating_light: **Attention #2**
Make sure you download the correct version! The older (legacy) version 1.0 differs in terms of implemented features and Graphical User Interface.
:::
- [x] [**Download and install LPJ-GUESS Education version 3.0**](http://web.nateko.lu.se/lpj-guess/education/) :link:
- [x] [**Download / update / install R and R-Studio**](https://www.rstudio.com/products/rstudio/download/#download) :link:
- [x] [**Download the lab R scripts**](https://uio.instructure.com/courses/30419/files/folder/Block%203%20-%20Lab%20materials) :link:

:::warning
:rotating_light: **Attention #3**
The R scripts use relative paths for creating analysis plots. Make sure to not change its internal structure (--> do not move or rename folders/files).
:::

## 2.1 Installation

:ballot_box_with_check: **Step 1**:
Download JPJ-GUESS Education version 3.0. Go to http://web.nateko.lu.se/lpj-guess/education/ and **click on this button** (**NOT** the legacy version link below):

<div style="text-align: center;">

![](https://i.imgur.com/MfKfAqc.png)

</div>

:ballot_box_with_check: **Step 2**:
Execute the `LPJGE.msi` installer (requires administrator rights) and follow the instructions.

:ballot_box_with_check: **Step 3**:
Open the "LPJ-GUESS Education 2013" application (e.g. using the Windows search bar).

:ballot_box_with_check: **Step 4**: 
Create a new "workspace". There will either be a prompt when you start the program, or you can change it using the menu:

<div style="text-align: center;">

![](https://i.imgur.com/lpclaA3.png)

:arrow_double_down: 

</div>

<div style="text-align: center;">

![](https://i.imgur.com/l9AIuTC.png)

</div>


Pick a destination and use the suggested name "Guesswork" when prompted. Here, the program will automatically create instruction / initialization files that are required to run simulations. We will also use this workspace to store the simulation output. A "fresh" workspace will contain the following files:
```
cohort.inz
lpjguess_workspace.ini
population.inz
```
More on this in the next section :100:.

<hr>

## 2.2 Simulation checklist
LPJ-GUESS Education is designed to run single-cell simulations. For this purpose, it also provides a user interface to pick pre-defined (simplified) climate forcing data for grid cells (i.e. geographic locations) around the globe. Before you run a simulation, you therefore need to make some choices for your model experiment.
:::info
- Which **grid cell location** do you want to simulate?
- Which **climate scenario** do you want to study?
- Which **model parameter** settings do you want to use? For example:
    - Which **vegetation mode** (population vs. cohort)?
    - How long should the **spin-up**/**static** (uses provided forcing data until vegetation is in a stable state) and **experimental** (aka **transient**, manipulates forcing data in a chosen manner) **phases** be?
    - Which **output variables** are you interested in?
    - etc.
:::
Don't worry if this sounds overwhelming, we will go through the options step by step. Also, refer to the lecture for additional information.


## 2.3 Generating climate forcing data

Open the "**GetClim**" tool here:

![](https://i.imgur.com/6c7Jr5k.png)

:ballot_box_with_check: **Step 1**: Choose a file name.
The climate forcing data will be generated as a `climate[_name].txt` file. Firstly, you need to pick a name and folder destination here:
<div style="text-align: center;">

![](https://i.imgur.com/AUwcNxE.png)

</div>


:::warning
:rotating_light: To keep your model runs organized, make sure you maintain a good folder structure, e.g. using seperate folders for each simulation. Also, your `climate[_name].txt` file should have a **descriptive name** - it should reflect the location you choose and the climatic changes you apply to the transient phase (e.g. CO2x2, RCP scenario, etc.; see below).
:::

:ballot_box_with_check: **Step 2**: Pick a grid cell location.

You can either use the GUI to select a grid cell with your mouse or manually enter geographic coordinates. 

:::warning
:rotating_light: The spatial resolution of the model is 0.5°, hence you need to round your coordinates to values ending in X.25 or X.75.
:::

<div style="text-align: center;">

![](https://i.imgur.com/LqL0NEd.png)

</div>

:::info
:bulb: We will not cover nitrogen deposition in this exercise, you should leave the value at "Medium". See for example <a href="https://doi.org/10.5194/esd-8-1121-2017">this study</a> for additional background information. The plot under "Climatology" gives a preview of the input data that will be generated.
:::

:ballot_box_with_check: **Step 3**: Define the simulation period and, if applicable, climate anomalies.

<div style="text-align: center;">

![](https://i.imgur.com/XlDWZCP.png)

</div>

:::warning
:rotating_light: Be sure to set the tick at "**Advanced mode**" to gain additional control over applying anomalies to the climate forcing data!
:::

For this exercise, you should use the standard values for **static phase** and **transient phase** (**1000** years and **200** years, respectively) and use the **Future Scenario Mode**. See the explanation of the different phases and modes below. The **Climate timeseries** radio button should be set to **replicate base series**.

Next click :arrow_right: **Enter anomalies**:

<div style="text-align: center;">

![](https://i.imgur.com/AY4A54v.png =320x)

</div>

Here, you enter the changes you want to apply to your base-line climate data during the transient phase. The fields should be more or less self-explanatory. Do not hesitate to ask questions if anything is unclear!

:::warning
:rotating_light: Be aware that you always enter relative changes. For example, if you want to study a predefined target $CO_2$ concentration in ppm, you need to subtract the base-line $CO_2$ first (**GetClim standard** is **280 ppm**, more on this in the exercise description).
:::

:::warning
:rotating_light: It is also important to understand that here, LPJ-GUESS is not coupled to a climate model. The values you enter for climatic variable anomalies **will not** influence other climatic variables. I.e., in contrast to "real-life", if you only alter the atmospheric $CO_2$ concentration, the temperature will remain unchanged.
:::

</br>

:::info
### Background: Simulation phases and simulation modes
The following figure illustrates the simulation phases in LPJ-GUESS Education:

<div style="text-align: center;">

![](https://i.imgur.com/mibFVkr.png)

Source: http://web.nateko.lu.se/lpj-guess/education/docs/protocol.html
</div>


In **Future scenario mode**, during the **static phase**, the modern day climate data is repeatedly **cycled** to bring the modelled vegetation into an equilibrium state. To simplify, once the sufficiently long static phase is over, you can assume to have modeled the "modern day" potential vegetation. In the **transient phase**, the anomalies you enter will linearly adapt the climate forcing to eventually reach the values you specified. It can therefore, for example, be used to model future potential vegetation changes under expected emission scenarios.

</br>

> :bulb: We will only use "**Future scenario mode**" in this exercise!


In **Paleo mode**, the anomalies you pick will be applied to the **static phase**. I.e., modern day climate data +/- anomalies are **cycled** until the vegetation reaches equilibrium and the static phase ends. During the **transient phase**, the climate data will linearly approach the initial modern day forcing.

</br>

> :bulb: To better understand what is going on, it will help to inspect the generated `climate.txt` files. A simple R script to do this is included in the lab materials.

:::warning
:rotating_light: <a id="climate-data-simple"></a>Note that the climate data generated with GetClim are **very** simplified. They represent "modern day" conditions for a 10 year period. Anomalies will linearly transform this data as shown in the figure, which is also an oversimplification.
:::


:ballot_box_with_check: **Step 4**: Click on Generate File :tada:

</br>

---
---
---

## 2.4 Defining the LPJ-GUESS parameters
The next step is to set the parameters that define the "general model behaviour". A detailed discussion of all parameters is out of the scope of this exercise, please refer to the <a href="http://web.nateko.lu.se/lpj-guess/education/docs/starthere.html">user guide</a> for additional information.

:ballot_box_with_check: **Step 1**: Pick a vegetation dynamics mode.

This defines the way PFT population dynamics are modelled. More info here: http://web.nateko.lu.se/lpj-guess/education/docs/vegmode.html. In this exercise, we will use the **cohort mode**. Go to your `Guesswork/` folder and open `cohort.inz` using your favorite editor (e.g. Notepad++).
:::info
:bulb: Note that `cohort.inz` and `population.inz` have exactly the same basic file structure. They only use different pre-defined settings for the respective parameter settings for convenience. 
:::

:::warning
:rotating_light: If you want to be able to reproduce your results, you should keep a copy of the final `.inz` file you used for a model run. In this exercise, we will use the same model parameter settings for all simulations and you will not need to store additional copies.
:::

---
:ballot_box_with_check: **Step 2**: Adapt the model output as follows.

Change the values in the following fields (note the line numbers at the left):

![](https://i.imgur.com/RZmmfw2.png)

This will ensure we have **yearly output** in the generated output files, which will make calculating and interpreting statistics a bit easier. Downside: larger file sizes. Also note that in the example, we specify a (manually created) folder for the model output (defined after `outputdirectory`, <span style="color:red">see warning below</span>). For the exercise, we will directly place the output in `Guesswork`, leave the value at `./`.
:::warning
:rotating_light: LPJ-GUESS Education does not automatically create the specified output directory. You have to create it manually before you run your simulation, otherwise it will throw an error!
:::

:::warning
:rotating_light: Be careful to change the output directory before each new model run! Otherwise, old files will be overwritten!
:::

Also note that you can add additional output files that are not produced by default. To do this, navigate to a variable name that currently is followed by an empty string (e.g. `file_agpp ""`) and add a file name of your choice that, by convention, ends in `.out`. For this exercise, alter the following line to produce file output for annual GPP:

![](https://i.imgur.com/IirXEAj.png)

---

:ballot_box_with_check: **Step 3**: Additional settings.

Also adapt the following settings in the specified lines:
```
116: npatch 50    # Reduce number of replicate patches (see lecture)
                  # Less patches = faster simulation (but statistically "less sound")
121: iffire 0     # Disables random wild fires
```

:ballot_box_with_check: **Step 4**: Familiarize yourself with the other parameters and save the settings file :tada:.

</br>

---
---
---

## 2.5 Starting the simulation

Now that we have defined:
- The climate forcing for a chosen grid cell location
- The LPJ-GUESS parameter settings

It is time to get things running! 

:ballot_box_with_check: **Step 1**: Open the "Run model" window.

Click here:

![](https://i.imgur.com/qmYY9PZ.png)

In the pop-up window, specify your `climate[_name].txt` file and the `population.inz` file that we created/edited. Remember that hint about chosing good, descriptive file names?

![](https://i.imgur.com/P3mSTFW.png)

:ballot_box_with_check: **Step 2**: Click "Run" :tada:!

</br>

---
---
---

## 2.6 Monitoring the simulation
While the model runs, you will see that the program displays the current state of a subset of pre-defined variables in real time.

![](https://i.imgur.com/voZCQ36.png)

You will see different output variables plotted against the current simulation years (x-axis), color coded by PFT. Notice the sudden drop in the variables after 100 years? Don't worry about it, this is implemented on purpose to establish baseline nitrogen pools. The behavior can be controlled with the following parameter in `[population_mode].inz`:

`line 139: freenyears 100          ! yrs to spin up without N limitation (needed to build up N pool)`

> [...] To overcome this, we followed the standard protocol, which is to run LPJ-GUESS for 100 years without N limitation but with normal N deposition to build up the N pools. After 100 years there is sufficient N in the pools, but the vegetation is inconsistent with the desired state as it has been growing without N limitation. Therefore, the vegetation is removed (and the C and N put into the litter pools), and the vegetation is allowed to regrow (this time with N limitation enabled) [...]. Source: Forrest et al. 2020 (https://doi.org/10.5194/gmd-13-1285-2020).

:::info
:bulb: Note that you can manually adjust the display of real-time plots under `Options`:

![](https://i.imgur.com/7bl5RRK.png)

For instance, to display a 3-D graph of PFT-growth on your defined number of `patches` (in the `[population_mode].inz` file), set a tick under `Graphics` -> `Vegetation 3D`. It will look similar to this:

![](https://i.imgur.com/ewupgWy.png =450x)


:rotating_light: Be aware that enabling this feature may slightly slow down your simulation!
:::

After the simulation is finished, a typical output plot looks like this:

<div style="text-align: center;">

![](https://i.imgur.com/1EIPcSt.png)

</div>

The grey area represents the **spin-up / static phase**. Note that the vegetation reaches a "more or less" stable state (=equilibrium) only after approx. 600 years. This is why you carefully need to pick the value of the `Static phase years` setting!

The red area shows the response of the vegetation to your entered climatic anomalies during the **transient phase**. In this example, the grid cell location (coords. corresponding to Oslo) projects highest net primary production for the **TeBS** PFT (temperate, broadleaved tree, summer-green, shade tolerant) after reaching equilibrium in "current" climatic conditions. Under the chosen climate scenario, after approximately 150 years the **IBS** PFT (boreal/temperate, broadleaved tree, summer-green, shade intolerant) becomes more productive in terms of NPP.

:::info
:question: **Realistic? Species? Reasons?** **Other PFTs?** **Implications?** This is the fun part!
:::

---


## 2.7 Analyzing the output

In addition to the direct graphic output in the interface, the raw numbers will also be written to output files located in the output directory you specified in <a href="#24-Defining-the-LPJ-GUESS-parameters">step 2.4</a>. Per convention, the delimited text files are named `[variable_name].out`. You can read them in with your favorite statistical tools (Python, R, Excel, etc.) and plot/analyse the output variables you are interested in. An R script that produces simple "30 year average NPP per PFT" plots is included in the course materials.

### Example plots

$CO_2$ **anomaly**

![](https://i.imgur.com/FCTcT44.png =500x)

**ANPP**

![](https://i.imgur.com/XvlcNWR.png =300x)





---
---

# 3 Exercises

You will be devided into 4 groups. Each group will study one ecosystem type, i.e. one given grid cell location. Using LPJ-GUESS Education, the aim is to tackle the following questions:

1. How can the modelled <a href="#Background-Simulation-phases-and-simulation-modes">"modern day" potential vegetation in equilibrium state</a> be characterized in your grid cell?
    1. Do you think the model output makes sense?
    2. In a real study: how would you evaluate these results and what kinds of questions could you address?
2. You will be given two future "climate anomaly" scenarios (see tables below).
    1. How does the modelled vegetation change in your grid cell **until the year 2100** (modern day + 100 years) under the two different scenarios?
    2. Do you think the vegetation is in a new equilibrium state? What could be the implications for the earth system?

:::warning
:rotating_light: Note that you need to execute only two model runs. The baseline output ("modern day equilibrium") will be included in both and should only differ slightly due to stochastic effects (random disturbances).
:::


### Groups:

|  Group 1  |     Group 2      |     Group 3      |    Group 4    |  Group 5   |
|:---------:|:----------------:|:----------------:|:-------------:|:----------:|
|           |                  |                  |               |            |
|           |                  |                  |               |            |
|           |                  |                  |               |            |
|           |                  |                  |               |            |
|           |                  |                  |               |            |
| :snowman: | :evergreen_tree: | :deciduous_tree: | :ear_of_rice: | :umbrella: |

Details about your biome (<a href="https://en.wikipedia.org/wiki/Biome">Wikipedia</a>) below. To check out satellite images of the coordinates on Google Maps, click on the pins:

|                  |        Biome        | Longitude | Latitude |                                 Map                                 |
|:----------------:|:-------------------:|:---------:|:--------:|:-------------------------------------------------------------------:|
|    :snowman:     |       Tundra        |  97.25 E  | 71.75 N  | <a href="https://goo.gl/maps/1dTBXEMpAXL8V4mGA">:round_pushpin:</a> |
| :evergreen_tree: |        Taiga        |  49.25 E  | 63.25 N  | <a href="https://goo.gl/maps/ZvVCNtmaAE3EXhyT6">:round_pushpin:</a> |
| :deciduous_tree: |  Temperate forest   |  9.25 E   | 50.25 N  | <a href="https://goo.gl/maps/gihCdxAXYkgd6cmt8">:round_pushpin:</a> |
|  :ear_of_rice:   |  Temperate steppe   | 103.75 W  | 47.25 N  | <a href="https://goo.gl/maps/H78N8Ks3wdyvTSri7">:round_pushpin:</a> |
|    :umbrella:    | Tropical rainforest |  66.75 W  |  4.75 S  | <a href="https://goo.gl/maps/iMDno967wk6mZRt4A">:round_pushpin:</a> |

### Climate scenarios:

We will use two <a href="https://en.wikipedia.org/wiki/Representative_Concentration_Pathway">IPCC Relative Concentration Pathway</a> scenarios, **RCP2.6** and **RCP8.5**, to have somewhat realistic values for our climate anomalies. Due to the time limitations, feel free to pick approximate values for the projected changes in your grid cell locations from the maps in the figures. We will assume that our generated base-line climate data reflect the same reference starting period (1986-2005) as is used in the figures and table below.

:::warning

:rotating_light: If the model runs take too long to terminate, feel free to split up the two simulations between your group members.

:::

**Temperature:**
![](https://i.imgur.com/VS59cMu.png)

**Precipitation:**
![](https://i.imgur.com/nrvCwq1.png)
:::warning
:rotating_light: **Attention**! Note that the precipitation change is given in % in the figures, but must be given in mm in LPJ-GUESS Education. Dirty quick hack: eyeball a mean precipitation value from the preview climate data plot you see (or calculate it from the climate .txt file) and then use the % value from the map to calculate an absolute value.
:::

**$CO_2$**:
![](https://i.imgur.com/UP4oXFP.png)

:::warning
:rotating_light: For entering the relative $CO_2$ change anomaly, use the standard GetClim value (280 ppm) and subtract it from the target concentration in 2100 (e.g. 985 - 280 = 705 for RCP8.5).
:::

## 3.1 Ancillary figures

**LPJ-GUESS PFTs:**
<div style="text-align: center;">

![](https://i.imgur.com/U2noGZQ.png =350x)

</div>

![](https://i.imgur.com/I8L6mWV.png)

![](https://i.imgur.com/gQCY56F.png)
<a href="https://academic.oup.com/bioscience/article/54/6/547/294347?login=true"> (Running et al., 2004)</a>

**Terminlogy**
- Autotrophic respiration (R$_a$)
- Heterotrophic respiration (R$_h$)
- Ecosystem respiration (R$_{eco}$)
- Gross Primary productivity (GPP)
- Net Primary Productivity (NPP)
    - NPP = GPP - R$_a$
- Net Ecosystem Exchange (NEE) / Net Ecosystem Production (NEP)
    - NEE = NEP = NPP - R$_h$
- <a href="https://en.wikipedia.org/wiki/Leaf_area_index">Leaf Area Index</a> (LAI)

<div style="text-align: center;">

![](https://i.imgur.com/9koNN3j.png =400x)

**Source**: By M. Campioli, Y. Malhi, S. Vicca, S. Luyssaert, D. Papale, J. Peñuelas, M. Reichstein, M. Migliavacca, M. A. Arain, I. A. Janssens - M. Campioli et. al.: “Evaluating the convergence between eddy-covariance and biometric methods for assessing carbon budgets of forests”; Nature Communications doi:10.1038/ncomms13717, CC BY 4.0, https://commons.wikimedia.org/w/index.php?curid=85907219 

</div>



## 3.2 Results

Add the figures that the R scripts generate into the boxes below. Simply open the figures in a viewer, right click, copy, and then paste them in the box using `CTRL+V`. Please also appoint someone to quickly present the results. Feel free to add additional plots as you wish!

:::success
# Tundra :snowman: 

### RCP2.6


**Temperature**: ![](https://i.imgur.com/pYo7dAw.png)

**Precipitation**: ![](https://i.imgur.com/GqYf4cU.png)

**CO$_2$**: ![](https://i.imgur.com/fYFDhUp.png)

## RCP8.5

**Temperature**: ![](https://i.imgur.com/gvxVYvy.png)

**Precipitation**: ![](https://i.imgur.com/kiI6MPA.png)

**CO$_2$**: ![](https://i.imgur.com/2enznBT.png)


## Compare output variables

### Net primary productivity (NPP)

<div style="display:flex; flex-direction: row;">

<div style="width:48%; border: 2px solid grey; padding: 5px;">
<h4 style="margin-top:0;">RCP 2.6</h4>
<div style="text-align-last: center;">

ADD YOUR PLOT HERE (Copy image, delete this text and CTRL+V)!

</div>
</div>
    
<div style="width:48%; border: 2px solid grey; padding: 5px; margin-left: auto;">
<h4  style="margin-top:0;">RCP 8.5</h4>
<div style="text-align-last: center;">

ADD YOUR PLOT HERE (Copy image, delete this text and CTRL+V)!

</div>
</div>
    
</div>

### Leaf area index (LAI) 

<div style="display:flex; flex-direction: row;">

<div style="width:48%; border: 2px solid grey; padding: 5px;">
<h4 style="margin-top:0;">RCP 2.6</h4>
<div style="text-align-last: center;">

![](https://i.imgur.com/4h78vfJ.png)

</div>
</div>
    
<div style="width:48%; border: 2px solid grey; padding: 5px; margin-left: auto;">
<h4  style="margin-top:0;">RCP 8.5</h4>
<div style="text-align-last: center;">

![](https://i.imgur.com/cm8d7v5.png)

</div>
</div>
    
</div>

:::



---

---

:::success
# Taiga :evergreen_tree:

### RCP2.6

**Temperature**:

**Precipitation**: 

**CO$_2$**: 

### RCP8.5

**Temperature**: 


**Precipitation**: 


**CO$_2$**: 


## Compare output variables

### Net primary productivity (NPP)

<div style="display:flex; flex-direction: row;">

<div style="width:48%; border: 2px solid grey; padding: 5px;">
<h4 style="margin-top:0;">RCP 2.6</h4>
<div style="text-align-last: center;">

ADD YOUR PLOT HERE (Copy image, delete this text and CTRL+V)!

</div>
</div>
    
<div style="width:48%; border: 2px solid grey; padding: 5px; margin-left: auto;">
<h4  style="margin-top:0;">RCP 8.5</h4>
<div style="text-align-last: center;">

ADD YOUR PLOT HERE (Copy image, delete this text and CTRL+V)!

</div>
</div>
    
</div>

### Leaf area index (LAI) 

<div style="display:flex; flex-direction: row;">

<div style="width:48%; border: 2px solid grey; padding: 5px;">
<h4 style="margin-top:0;">RCP 2.6</h4>
<div style="text-align-last: center;">

![](https://i.imgur.com/4h78vfJ.png)

</div>
</div>
    
<div style="width:48%; border: 2px solid grey; padding: 5px; margin-left: auto;">
<h4  style="margin-top:0;">RCP 8.5</h4>
<div style="text-align-last: center;">

![](https://i.imgur.com/cm8d7v5.png)

</div>
</div>
    
</div>

:::



---

---

:::success
# Temperate forest :deciduous_tree: 

### RCP2.6

**Temperature**:

**Precipitation**: 

**CO$_2$**: 

### RCP8.5

**Temperature**: 


**Precipitation**: 


**CO$_2$**: 


## Compare output variables

### Net primary productivity (NPP)

<div style="display:flex; flex-direction: row;">

<div style="width:48%; border: 2px solid grey; padding: 5px;">
<h4 style="margin-top:0;">RCP 2.6</h4>
<div style="text-align-last: center;">

ADD YOUR PLOT HERE (Copy image, delete this text and CTRL+V)!

</div>
</div>
    
<div style="width:48%; border: 2px solid grey; padding: 5px; margin-left: auto;">
<h4  style="margin-top:0;">RCP 8.5</h4>
<div style="text-align-last: center;">

ADD YOUR PLOT HERE (Copy image, delete this text and CTRL+V)!

</div>
</div>
    
</div>

### Leaf area index (LAI) 

<div style="display:flex; flex-direction: row;">

<div style="width:48%; border: 2px solid grey; padding: 5px;">
<h4 style="margin-top:0;">RCP 2.6</h4>
<div style="text-align-last: center;">

![](https://i.imgur.com/4h78vfJ.png)

</div>
</div>
    
<div style="width:48%; border: 2px solid grey; padding: 5px; margin-left: auto;">
<h4  style="margin-top:0;">RCP 8.5</h4>
<div style="text-align-last: center;">

![](https://i.imgur.com/cm8d7v5.png)

</div>
</div>
    
</div>

:::



---

---

:::success
# Temperate steppe :ear_of_rice: 

### RCP2.6

**Temperature**:

**Precipitation**: 

**CO$_2$**: 

### RCP8.5

**Temperature**: 


**Precipitation**: 


**CO$_2$**: 


## Compare output variables

### Net primary productivity (NPP)

<div style="display:flex; flex-direction: row;">

<div style="width:48%; border: 2px solid grey; padding: 5px;">
<h4 style="margin-top:0;">RCP 2.6</h4>
<div style="text-align-last: center;">

ADD YOUR PLOT HERE (Copy image, delete this text and CTRL+V)!

</div>
</div>
    
<div style="width:48%; border: 2px solid grey; padding: 5px; margin-left: auto;">
<h4  style="margin-top:0;">RCP 8.5</h4>
<div style="text-align-last: center;">

ADD YOUR PLOT HERE (Copy image, delete this text and CTRL+V)!

</div>
</div>
    
</div>

### Leaf area index (LAI) 

<div style="display:flex; flex-direction: row;">

<div style="width:48%; border: 2px solid grey; padding: 5px;">
<h4 style="margin-top:0;">RCP 2.6</h4>
<div style="text-align-last: center;">

![](https://i.imgur.com/4h78vfJ.png)

</div>
</div>
    
<div style="width:48%; border: 2px solid grey; padding: 5px; margin-left: auto;">
<h4  style="margin-top:0;">RCP 8.5</h4>
<div style="text-align-last: center;">

![](https://i.imgur.com/cm8d7v5.png)

</div>
</div>
    
</div>

:::



---

---

:::success
# Tropical rainforest :umbrella:

### RCP2.6

**Temperature**:

**Precipitation**: 

**CO$_2$**: 

### RCP8.5

**Temperature**: 


**Precipitation**: 


**CO$_2$**: 


## Compare output variables

### Net primary productivity (NPP)

<div style="display:flex; flex-direction: row;">

<div style="width:48%; border: 2px solid grey; padding: 5px;">
<h4 style="margin-top:0;">RCP 2.6</h4>
<div style="text-align-last: center;">

ADD YOUR PLOT HERE (Copy image, delete this text and CTRL+V)!

</div>
</div>
    
<div style="width:48%; border: 2px solid grey; padding: 5px; margin-left: auto;">
<h4  style="margin-top:0;">RCP 8.5</h4>
<div style="text-align-last: center;">

ADD YOUR PLOT HERE (Copy image, delete this text and CTRL+V)!

</div>
</div>
    
</div>

### Leaf area index (LAI) 

<div style="display:flex; flex-direction: row;">

<div style="width:48%; border: 2px solid grey; padding: 5px;">
<h4 style="margin-top:0;">RCP 2.6</h4>
<div style="text-align-last: center;">

![](https://i.imgur.com/4h78vfJ.png)

</div>
</div>
    
<div style="width:48%; border: 2px solid grey; padding: 5px; margin-left: auto;">
<h4  style="margin-top:0;">RCP 8.5</h4>
<div style="text-align-last: center;">

![](https://i.imgur.com/cm8d7v5.png)

</div>
</div>
    
</div>

:::