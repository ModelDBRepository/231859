*How to use this code*

The experiments are run using SkinBrian.py, a framework for running experiments of cells in simple topologies. As of yet, only the tubular topology is implemented and only the integrate-and-fire cell model. 

It allows for easy scanning of parameter spaces; the code below <if __name__ == '__main__':> is fairly self-explanatory. Simply create your own set of model parameters, make sure the model object takes that set as input, and press F5! Also note that there are a few global variables which may impact your performance. Debugging works best on just 1 core, then the multiprocessing module is not called.

There is a caveat when it comes to using one output directory more than once: there is a built-in method to help when very long runs with many parameter combinations break down (as they inevitably will), which will continue where the runs broke down. This goes by run index; this may be different when other parameters are set, so make sure to clear out the output directory if the parameter space and model code is NOT _exactly_ the same.

All other files are for analysis and visualization. They require you to set the correct directory for the output to look into. stripeAnalysis.py is the main file for overall quantitative analysis; SkinBrianVectorViz.py is the main file for detailed, single-run analysis. Both are called in MikadoPaperVizScript.py upon which the analyses in the 'Short and Random'-paper are based.

It requires Python 2.7 and Brian2. All other dependencies are also dependencies of Brian2. 