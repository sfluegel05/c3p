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

    # Define the core beta-carboline structure, allowing for saturation, using SMARTS.
    # The core pattern captures the fused ring system with a nitrogen in the 5-membered ring.
    # The SMARTS string should match both aromatic and hydrogenated versions.
    beta_carboline_pattern = Chem.MolFromSmarts('[cH1,CH2,CH3]1[cH1,CH2,CH3][cH1,CH2,CH3]2[nH][cH0,CH1,CH2,CH3]3[cH1,CH2,CH3]([cH1,CH2,CH3]2[cH1,CH2,CH3]1)[c,n][cH1,CH2,CH3][cH1,CH2,CH3]3')
    
    spiro_pattern = Chem.MolFromSmarts('[C;X4]')

    if mol.HasSubstructMatch(beta_carboline_pattern):
       if mol.HasSubstructMatch(spiro_pattern):
         return True, "Contains the beta-carboline core structure (including spirocyclic derivatives)"
       return True, "Contains the beta-carboline core structure"

    return False, "Does not contain the beta-carboline core structure"