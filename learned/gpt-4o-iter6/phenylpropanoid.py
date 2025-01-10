"""
Classifies: CHEBI:26004 phenylpropanoid
"""
from rdkit import Chem

def is_phenylpropanoid(smiles: str):
    """
    Determines if a molecule is a phenylpropanoid based on its SMILES string.
    Phenylpropanoids are organic aromatic compounds with a structure based on a phenylpropane skeleton.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phenylpropanoid, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string into an RDKit molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define key features to identify a phenylpropanoid
    # A basic phenylpropanoid can be represented as "c1ccccc1CCC" - a phenyl ring followed by 3-carbons (propane).
    phenylpropane_pattern = Chem.MolFromSmarts("c1ccccc1CCC")
    
    # Check for the presence of the phenylpropane core
    if not mol.HasSubstructMatch(phenylpropane_pattern):
        phenyl_pattern = Chem.MolFromSmarts("c1ccccc1")
        if not mol.HasSubstructMatch(phenyl_pattern):
            return False, "No phenyl ring found"
        
        # Check for presence of alternative phenyl linked with small carbon structures, e.g., CO-coupled systems
        coupled_pattern = Chem.MolFromSmarts("[cH1]-[CX3](=[OX1])-[O,N]", asQuery=True)
        if not mol.HasSubstructMatch(coupled_pattern):
            return False, "Missing typical phenylpropanoid linkage (e.g., phenolic or carbonyl coupled system)"
    
    return True, "Structure consistent with phenylpropanoid, containing phenyl/group and typical linkages"

# Note: This is a heuristic approach due to the diverse nature of phenylpropanoids.