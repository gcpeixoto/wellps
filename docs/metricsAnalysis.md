## `lib/metricsAnalysis`

**Gustavo Oliveira, @LaMEP/UFPB, gustavo.oliveira@ci.ufpb.br** 

Documentation to understand the fields belonging to the `.csv` files generated from the function `lib/metricsAnalysis`.

### A brief comment on hierarchy

`metricsAnalyzer` handles data that follow the following hierarchy as for reservoir modelling: 

**cell**: the underlying entity of a grid model. For popular reservoir simulators, they are each of the corner-point grid _blocks_. 

**cluster**: a set of face-connected grid cells that form an irregular structure. 

**field**: all the cells cells that make up the grid model. 

For instance, if the model has 1000 cells, any file related to `field` must contain information for all these 1000 cells. However a `cluster` is a volume inside the field with less cells, let us say 50.


### File semantics

Below, we explain the semantics used to name files related to `field` and `cluster` hierarchy. These files store information of the parameter `DRT` (_discrete rock typing_), used by `wellps` to discretize the reservoir model.

#### Placeholders 

Placeholders are denoted by this symbol: `<...>`. Currently, the following options are available: 

`<average>`: one of these three strings: _arithmetic_, _geometric_, or _normalized_.  

`<base>`: one of these two strings: _log10_ (base-10 logarithm) or _ln_ (natural logarithm).

`<D>`: a positive integer value.

`<C>`: cluster ID. Note that we may have many clusters for the same `<D>` value.

<!--`<name>`: a specific nomenclature related to the hierarchy which is explained in the sections below.-->


#### `field` related file

**File template:** `DRT_<average>_<base>_<D>_PerformanceTable.csv`

*Description:* overview of all the clusters related to DRT = D.

*Fields:*

`cluster`: _cluster ID_. The indexing is the natural one: 1,2,3,...

`nce`: number of cluster's cells in descending order.

`s`: slope value computed by linear regression.

`R2`: coefficient of determination computed by linear regression.

`performance`: boolean value that indicates: `0`: low performance; `1`: high performance. 


##### Meaning of "performance"

Here, performance is a decision variable used to separate hydraulic flow units (the clusters) that pass in the linear regression analytical test. This test is based on values of `s` and `R2`. For instance, a cluster is considered a "high performance" one if and only if both the conditions below are satisfied: 

1. $1.0 - \epsilon \le s \le 1.0 + \epsilon$  
2. $0.9 \le R^2 \le 1.0$  

where $\epsilon = \mathcal{O}(10^{-t})$ is a tolerance chosen by the user. For oil recovery simulations, we select only cells belonging to high-performance clusters. Hence, if this field is totally filled by 0 values in the `PerformanceTable` file, this means that none of those clusters should be harnessed. 

#### `cluster` related files

**File template:** `DRT_<average>_<base>_<D>_Cluster_<C>_MetricsAnalytics.csv`

*Description:* quantitative information on centrality metrics over the cluster cells. 

*Fields:*

`I`: I-coordinate (logical index) of the cluster cell.

`J`: J-coordinate (logical index) of the cluster cell.

`K`: K-coordinate (logical index) of the cluster cell.

`closeness`: closeness centrality of the graph vertex mapped onto this cell.

`betweeness`: betweeness centrality of the graph vertex mapped onto this cell.

`degree`: degree centrality of the graph vertex mapped onto this cell.

##### Remark on `betweeness`

It was observed that negative values were found for `betweeness` in some cases and they might appear in some .csv file. We stress that this is a anomalous behaviour. The centrality metrics are computed from `SNAP` and betweeness must not be negative. 

While this investigation is to be concluded, we should say that beetweeness is not yet a critical quantity to our tests because it is not used to determine the optimal positions. 

**File template:** `DRT_<average>_<base>_<D>_Cluster_<C>_MetricsAnalyticsMinMax.csv`

*Description:* detailed information on closeness centrality and topology of the cluster `C`. 

*Fields:*

`maxClo`: value of maximum closeness centrality found over the cluster.

`MIglob`: _global_ I-coordinate (logical index) of the maximum closeness cell.

`MJglob`: _global_ J-coordinate (logical index) of the maximum closeness cell.

`MKglob`: _global_ K-coordinate (logical index) of the maximum closeness cell.

`MIloc`: _local_ I-coordinate (logical index) of the maximum closeness cell.

`MJloc`: _local_ J-coordinate (logical index) of the maximum closeness cell.

`MKloc`: _local_ K-coordinate (logical index) of the maximum closeness cell.

`minClo`: value of minimum closeness centrality found over the cluster.

`mIglob`: _global_ I-coordinate (logical index) of the minimum closeness cell.

`mJglob`: _global_ J-coordinate (logical index) of the minimum closeness cell.

`mKglob`: _global_ K-coordinate (logical index) of the minimum closeness cell.

`mIloc`: _local_ I-coordinate (logical index) of the minimum closeness cell.

`mJloc`: _local_ J-coordinate (logical index) of the minimum closeness cell.

`mKloc`: _local_ K-coordinate (logical index) of the minimum closeness cell.

`Imin`: minimum I-coordinate (logical index) of the cluster bounding box.

`Imax`: minimum I-coordinate (logical index) of the cluster bounding box.

`Jmin`: minimum J-coordinate (logical index) of the cluster bounding box.

`Jmax`: minimum J-coordinate (logical index) of the cluster bounding box.

`Kmin`: minimum K-coordinate (logical index) of the cluster bounding box.

`Kmax`: minimum K-coordinate (logical index) of the cluster bounding box.

`s`: slope value computed by linear regression.

`R2`: coefficient of determination computed by linear regression.

`performance`: boolean value that indicates: `0`: low performance; `1`: high performance. 