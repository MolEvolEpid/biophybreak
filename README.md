# biophybreak
Code from Lundgren et al "Combining biomarker and virus phylogenetic models improves epidemiological source identification" published in PLOS Computational Biology 2022 Aug 26;18(8):e1009741 DOI: https://doi.org/10.1371/journal.pcbi.1009741

Modified from phybreak, the package implementing the method described in Klinkenberg et al (2016), doi: http://dx.doi.org/10.1101/069195

Outbreak reconstruction with sequence data and biomarkers

Workflow:

* use 'mbm.predict' to find probabilitity distributions for infection times using biomarkers

* enter data and priors by constructing an object of S3-class 'phybreak', with function 'phybreak'

* do mcmc-updates with functions 'burnin_phybreak' and 'sample_phybreak'; remove samples with 'thin.phybreak'

* access the 'phybreak'-object by get_phybreak-functions such as 'get_transtree', 'get_data', 'get_parameters'

* summarize the mcmc-chain with the functions 'ESS', 'transtree', 'infectorsets', 'phylotree', 'get_mcmc', 'get_phylo', 
  'phybreak.infector.posts'

* plotting with 'plot', 'plotTrans', 'plotPhylo', 'plotPhyloTrans', 'phybreak.plot.posteriors'


* it is possible to simulate data with 'sim_phybreak' or 'rtrans.tree' and 'sim.coal.phybreak'
