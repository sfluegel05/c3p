"""
Classifies: CHEBI:60834 beta-carbolines
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_beta_carbolines(smiles: str):
    """
    Determines if a molecule is a beta-carboline based on its SMILES string.
    A beta-carboline is a pyrido[3,4-b]indole or its hydrogenated derivatives, or spirocyclic variants.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a beta-carboline, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the core indole structure using SMARTS.
    indole_pattern = Chem.MolFromSmarts('c1cc2[nH]c(c1)c2')
    if not mol.HasSubstructMatch(indole_pattern):
        return False, "Does not contain the indole core."
    
    # Define the fused 6-membered ring structure (pyridine, pyrimidine, or derivatives) using SMARTS.
    # This pattern allows for saturation in the 6-membered ring.
    pyridine_pyrimidine_pattern = Chem.MolFromSmarts('[c,n]1ccccc1')

    #check for spiro atoms, which have a valence of 4.
    spiro_pattern = Chem.MolFromSmarts('[C;X4]')
    
    # Check if the molecule has both patterns. We are not checking explicit fusion here,
    # so it will likely overclassify.
    # In the previous try, we tried to use a single pattern, but failed due to the numerous variants of beta-carbolines.
    # Now we are explicitly looking for the two ring systems.
    if mol.HasSubstructMatch(pyridine_pyrimidine_pattern):
      
      # If a spiro atom is present then we should also return true.
       if mol.HasSubstructMatch(spiro_pattern):
           return True, "Contains the beta-carboline core structure (including spirocyclic derivatives)"
       return True, "Contains the beta-carboline core structure (not spirocyclic)"


    return False, "Does not contain the beta-carboline core structure"