# What is EnPACT?

EnPACT is a modeling strategy which fine tunes DNA sequence-based molecular predictors to new datasets. Generally, the steps to train an EnPACT model are as follows: 

1. Assemble a training dataset of molecular features (gene expression, ChIP-seq peaks, etc.) of which you want to be able to predict the genetic component of variation across features. This will be the "Y" or "output" into the EnPACT model during training.
2. Generate "reference epigenome" predictions at molecular feature loci by passing reference genome DNA sequence from each loci into the DNA sequence-based molecular predictor. These "reference epigenome" predictions will then become the "X" or "input" into the EnPACT model during training.
3. Train the EnPACT model of choice across loci in the genome to predict "Y" from "X" at each loci in the training dataset using standard model training practices.

Pretty simple. In practice, the "En" in "EnPACT" refers to "Enformer", the sequence-based predictor explored here. However, there is no reason that this should be the **only** model which can fulfill this role, so this code is intentionally modular to encourage exploration of other sequence-based predictors. The EnPACT model is also a flexible, but important, choice. Linear elastic net and a simple convolutional neural net are two currently implemented options, but exploration is again encouraged with modular code support. 


## Using EnPACT to explore genetic causes of disease

As disease 

###





## To-do

Store only a single copy of Enformer personalized predictions to reduce memory overhead.
Option to delete personalized EnPACT predictions after contruction of linearized models.
