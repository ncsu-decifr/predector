DeepLoc 1.0
===========

DeepLoc-1.0 predicts the subcellular localization of eukaryotic proteins. It can differentiate between 10 different localizations: Nucleus, Cytoplasm, Extracellular, Mitochondrion, Cell membrane, Endoplasmic reticulum, Chloroplast, Golgi apparatus, Lysosome/Vacuole and Peroxisome.

Publication
------------

Jose Juan Almagro Armenteros, Casper Kaae Soenderby, Soeren Kaae Soenderby, Henrik Nielsen, Ole Winther; DeepLoc: prediction of protein subcellular localization using deep learning, Bioinformatics, Volume 33, Issue 21, 1 November 2017, Pages 3387-3395,

More information about the method can be found at:

        http://www.cbs.dtu.dk/services/DeepLoc-1.0


Pre-installation
----------------

DeepLoc 1.0 will run and has been tested under Linux and OS X. The only prerequisite is to have python2.7 installed.


Installation
------------

You need to install the following dependencies: 
  Numpy
  Scipy
  Theano==1.01
  Lasagne==0.2.dev1

The installation procedure is:

  1. Install dependencies: 
  	pip install -r requirements.txt

  2. Install DeepLoc package:
        python setup.py install
     or locally:
	python setup.py install --user

  3. Test DeepLoc by running:
	deeploc -f test.fasta
     
     the result should look like the file: test_output.txt 


Running
--------

DeepLoc will be installed under the name 'deeploc'. It has 3 possible arguments:

 * -f, --fasta. Input in fasta format of the proteins.
 * -o, --output. Prefix to use for the output files.
 * -a, --attention. Save attention values in a txt file for each individual protein. 

It is also posible to run deeploc on python directly. This is an example code:

 import os
 os.environ['THEANO_FLAGS']='device=cpu,floatX=float32,optimizer=fast_compile'
 import numpy as np
 from DeepLoc.models import *
 from DeepLoc.utils import * 


 # Read fasta file into a numpy matrix
 ids, prot_seqs = parse_fasta('test.fasta')

 # Perform the prediction
 batch_size = 1
 out_ids, out_loc, out_mem = prediction(ids, prot_seqs, batch_size)

The variable out_loc contains a numpy array with the probability of the localization prediction.


Output
-------

The output is a tabular file with the following format:

 * 1st column: Protein ID.
 * 2nd column: Predicted location.
 * 3rd column: Probability of a protein being membrane-bound.
 * 4th-13th column: Probability for each of the individual localizations. 

Problems and questions
----------------------

In case of technical problems (bugs etc.) please contact packages@cbs.dtu.dk.

Questions on the scientific aspects of the SignalP method  should go to Henrik
Nielsen, hnielsen@bioinformatics.dtu.dk.



