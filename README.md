# PEP-Benchmark
The Benchmark of Promoter-Enhancer interaction prediction by AlphaGenome
dependencies:
- python == 3.13.13
- alphagenome == 0.6.1
- biopython == 1.8.7
- numpy == 2.4.6
- pandas == 2.3.3
- scikit-learn == 1.9.1
- matplotlib == 3.10.9

# Run step-by-step
1) Clone repository:

   '''
   git clone https://github.com/NuriddinovMA/PEP-Benchmark.git
   cd PEP-Benchmark
   '''
   
3) install dependenciyes.
4) download the hg38 genome and put it to the genomes folder. The genome used for this benchmark is available on [repository](https://genedev.bionet.nsc.ru/ftp/by_Project/alphagenome/genomes/)
5) generate you unique api key for excess to Alphagenome API (see: [Alphagenome documentaion](https://www.alphagenomedocs.com/installation.html) )
6) run alpha.py to generate variant scores for promoter-enhancer pairs
   
   '''
   python alpha.py --api_key [YOU_API_KEY]
   '''
   
8) run auc.py to calculate auc
   
   '''
   python auc.py
   '''
