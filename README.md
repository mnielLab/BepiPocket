# BepiPocket

Current deep learning methods, such as AlphaFold and Chai, can often create antibody-antigen (AbAg) structures with high confidence and accuracy.
However, these methods often fail to predict the correct antibody binding site, placing the antibody incorrectly on the antigen and converging on repeatedly predicting the same redundant binding mode.
Varying seeds and diffusion samples does not guarantee that a diverse set of binding modes explored.
**BepiPocket** is a simple approach for integrating B-cell epitope prediction tools, **BepiPred-3.0** and **DiscoTope-3.0**, to guide antibody-epitope restraints during Chai-1 structure prediction.
This vastly increase diversity of explored antibody binding sites as well as allowing Chai-1 to find substantialy more accurate AbAg structures with high confidence.

Preprint link 

### TODOS
* **TODO**: Upload MMseqs2 MSA Code. MSAs can be provide in the same format used as Chai-1 (.pqt). But are not automatically created on inference. Need to add this code
* **TODO**: Upload DiscoTope-3.0 Code. Not sure if inverse folding environement used for DiscoTope-3.0 is compatible with Chai-1 package environment. Need to sort this out.

## License 
BepiPocket is a tool developed by the Health Tech section at Technical University of Denmark (DTU). The code and data can be used freely by academic groups for non-commercial purposes.
If you plan to use these tools for any for-profit application, you are required to obtain a license (contact Morten Nielsen, morni@dtu.dk).
* If you use BepiPred-3.0 (BepiPocket) to guide to create antibody-epitope pockets, please get a BepiPred-3.0 license.
* If you use DiscoTope-3.0 (DiscoPocket) to guide to create antibody-epitope pockets, please get a DiscoTope-3.0 license.
* The tool uses Chai-1 for structural inference of antibody-antigen compexes. Chai-1 is not open-source, and you are therefore also required to get license from Chai Discovery. 

## Graphical Abstract
![Screenshot](GraphicalAbstract.png)

## Installation 

### Create Conda Environment


### Install Pip Packages 

## Usage 

### Inputs

### Outputs

### Example

## Cite
If you found this tool useful in your research, please cite:<br>
[Pocket Restraints Guided by B-Cell Epitope Prediction improves Chai-1 Antibody-Antigen Structure Modeling](https://doi.org/10.1101/2025.09.17.676770)

