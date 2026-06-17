# CellModeller Parameters Guide



## 1\. Simulation-Level Parameters

|Parameter|Type|Description|
|-|-|-|
|`max\\\_cells`|int|Maximum number of cells before simulation stops|
|`jitter\\\_z`|bool|`false` = 2D colony, `true` = 3D|
|`gamma`|float|Frictional drag on cell growth. Higher = cells push harder|
|`pickle\\\_steps`|int|Save state every N steps (lower = more frequent)|
|`random\\\_seed`|int|Seed for reproducibility|



## 2\. Cell Type Parameters

|Parameter|Type|Description|
|-|-|-|
|`display\\\_name`|string|Human-readable label|
|`color`|\[R, G, B]|Visualization color, values 0.0–1.0|
|`growth\\\_rate`|float|How fast the cell grows each timestep|
|`division\\\_length`|float|Cell length that triggers division|
|`division\\\_noise`|float|Random variation added to division length|
|`initial\\\_pos`|\[x, y, z]|Starting position of the seed cell|
|`initial\\\_dir`|\[x, y, z]|Starting orientation of the seed cell|



## 3\. Biochemical Kinetics Parameters

|Parameter|Type|Description|
|-|-|-|
|`production\\\_rate`|float|Basal rate of protein production per timestep|
|`max\\\_production\\\_rate`|float|Maximum rate when fully activated|
|`degradation\\\_rate`|float|Fraction of protein degraded per timestep|
|`hill\\\_coefficient`|float|Cooperativity of activation/repression (n in Hill eq.)|
|`activation\\\_threshold`|float|Signal concentration for half-max activation|
|`repression\\\_threshold`|float|Protein level for half-max repression|



## 4\. Chemical Signaling Parameters

|Parameter|Type|Description|
|-|-|-|
|`enabled`|bool|Whether to use the Chemics signaling module|
|`grid\\\_size`|int|Resolution of the diffusion grid|
|`diffusion\\\_rate`|float|How fast each signal spreads (per signal)|
|`signal\\\_degradation\\\_rate`|float|How fast each signal breaks down|
|`boundary\\\_condition`|string|`"periodic"` or `"fixed"`|



## 5\. SBOL Mapping Overrides (optional)

By default, the converter will try to infer mappings automatically.
Use this section to correct or override those inferences.

|Parameter|Type|Description|
|-|-|-|
|`component\_to\_celltype`|dict|Map SBOL component display IDs to cell type indices|
|`signal\_component\_ids`|list|List of SBOL component IDs that are diffusible signals|
|`ignore\_component\_ids`|list|SBOL components to skip entirely (e.g. backbone, chassis)|



