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
        bool: True if the molecule is a nucleobase analogue, False otherwise
        str: Reason for classification
    """
    
    # Convert SMILES to RDKit Molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define improved SMARTS for purine and pyrimidine-like structures
    purine_smarts = Chem.MolFromSmarts("c1ncnc2[nH]ncnc12")
    pyrimidine_smarts = Chem.MolFromSmarts("c1[nH]cncnc1")
    
    # Check for nitrogen-containing heterocyclic rings characteristic of nucleobases
    if purine_smarts and mol.HasSubstructMatch(purine_smarts):
        return True, "Matches purine-like structure"
    elif pyrimidine_smarts and mol.HasSubstructMatch(pyrimidine_smarts):
        return True, "Matches pyrimidine-like structure"
    
    # If the central ring did not match, check for other characteristic ring structures
    # Check for presence of basic functional groups
    functional_groups = [
        Chem.MolFromSmarts("[CX3]=[OX1]"),  # Keto/Carbonyl group
        Chem.MolFromSmarts("[NX3;H2]"),     # Amino group
        Chem.MolFromSmarts("[OX2H]")        # Hydroxy group
    ]
    
    # Look for at least one nucleobase-typical functional group
    if any(mol.HasSubstructMatch(fg) for fg in functional_groups):
        return True, "Contains significant nucleobase-related functional groups"
    
    # None matched
    return False, "No nucleobase analogue characteristics found"