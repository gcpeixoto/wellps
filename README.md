# WELLPS - The _Well Placement Strategy Generator_

WELLPS is an in-house code used to generate well placement strategies
based on a couple of computational methods. 

Besides its own routines, it uses a few functions from the following 
third-party softwares

- [MRST](http://mrst.no), essentialy to read complex reservoir grids;
- [SNAP/C++](http://snap.stanford.edu), to compute graph-based centrality measures;

and the MATLAB-specific routines

- [Ellipsoid Fit](https://www.mathworks.com/matlabcentral/fileexchange/24693-ellipsoid-fit), by Y. Petrov, to compute 3D fitting operations.
- [Network Components](https://www.mathworks.com/matlabcentral/fileexchange/42040-find-network-components) by D. Larremore, to compute connected components.

Before running WELLPS, it is recommended to set the environment variables
`MRST_DIR` and `SNAP_DIR`, which should point to local paths on your computer
where MRST and SNAP are located. 

Setting `SNAP_DIR` is mandatory for C++ compilation. Firstly, compile SNAP by using `make all`.
Then, compile `graphMetrics` by using the in-code Makefile (remember to clean the `cpp/main/*.o` file). 
Finally, test the examples. Thenceforward, all things happen from the M-file `startup.m`.
If you are to use WELLPS for the first time, you need to recompile it for your platform accordingly.

WELLPS runs satisfactorily on UNIX-based systems but we have no experience with Windows-based execution.

WELLPS is rarely undergoing maintainance for lack of funding, but it is still 
being applied in oil and gas research. For now, it subsides studies on 
carbon storage and injection-related problems.

For research collaboration, get in touch with Dr. Oliveira.