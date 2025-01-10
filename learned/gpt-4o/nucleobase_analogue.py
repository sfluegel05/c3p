"""
Classifies: CHEBI:67142 nucleobase analogue
"""
from rdkit import Chem

def is_nucleobase_analogue(smiles: str):
    """
    Determines if a molecule is a nucleobase analogue based on its SMILES string.
    A nucleobase analogue is a molecule that can substitute for a normal nucleobase in nucleic acids.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nucleobase analogue, False otherwise
        str: Reason for classification
    """
    
    # Convert SMILES to RDKit Molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS for purine and pyrimidine-like structures
    purine_smarts = Chem.MolFromSmarts("c1ncnc2[nH]cnc12")
    pyrimidine_smarts = Chem.MolFromSmarts("c1cncnc1")
    
    # Check for heterocyclic aromatic rings with nitrogen atoms
    if not mol.HasSubstructMatch(pyrimidine_smarts) and not mol.HasSubstructMatch(purine_smarts):
        return False, "No purine or pyrimidine-like structure found"
    
    # Check for typical functional groups in nucleobase analogues
    functional_groups = [
        Chem.MolFromSmarts("C=O"),  # Keto group
        Chem.MolFromSmarts("NH2"),  # Amino group
        Chem.MolFromSmarts("OH")    # Hydroxy group
    ]
    
    # Search for presence of at least one additional functional group
    modifications = sum(mol.HasSubstructMatch(fg) for fg in functional_groups)
    if modifications == 0:
        return False, "No significant modifications found compared to normal nucleobases"
    
    return True, "Contains heterocyclic aromatic ring with nitrogen and additional functional groups typical of nucleobase analogues"