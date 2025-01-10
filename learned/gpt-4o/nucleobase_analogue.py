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
    
    # Define SMARTS patterns for purine-like structures with modifications
    purine_like_smarts = [
        Chem.MolFromSmarts("c1ncnc2[nH]cnc12"),      # Basic purine
        Chem.MolFromSmarts("c1nc[nH]c2c(ncnc12)"),  # Altered purine derivatives
    ]
    
    # Define SMARTS patterns for pyrimidine-like structures with modifications
    pyrimidine_like_smarts = [
        Chem.MolFromSmarts("c1[nH]cncnc1"),           # Basic pyrimidine
        Chem.MolFromSmarts("c1c[nH]cnc[nH]c1"),       # Variants with additional rings
        Chem.MolFromSmarts("c1ncc(=O)[nH]n1"),        # Typical keto and imidazole groups
        Chem.MolFromSmarts("c1nc[nH]cnc1"),           # Wider pyrimidine scaffolding
    ]
    
    # Check for purine-like or pyrimidine-like ring structures
    for pattern in purine_like_smarts:
        if mol.HasSubstructMatch(pattern):
            return True, "Matches purine-like structure"

    for pattern in pyrimidine_like_smarts:
        if mol.HasSubstructMatch(pattern):
            return True, "Matches pyrimidine-like structure"
    
    # Check for functional groups characteristic of nucleobase analogues
    functional_groups = [
        Chem.MolFromSmarts("[CX3]=[OX1]"),  # Keto/Carbonyl group
        Chem.MolFromSmarts("[NX3;H2]"),     # Amino group
        Chem.MolFromSmarts("[OX2H]"),       # Hydroxy group
        Chem.MolFromSmarts("[Cl,Br,I]")     # Halogen group, common in analogues
    ]
    
    # Consider presence of nucleobase-typical functional groups
    if any(mol.HasSubstructMatch(fg) for fg in functional_groups):
        if mol.GetNumAtoms() > 6:  # A reasonable count to filter out very small molecules
            return True, "Contains significant nucleobase-related functional groups"
    
    # None matched
    return False, "No nucleobase analogue characteristics found"