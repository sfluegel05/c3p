"""
Classifies: CHEBI:18303 phosphatidyl-L-serine
"""
from rdkit import Chem

def is_phosphatidyl_L_serine(smiles: str):
    """
    Classifies a molecule as a phosphatidyl-L-serine based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if the molecule is a phosphatidyl-L-serine, False otherwise
        str: Reason for classification
    """
    
    # Convert the SMILES string to an RDKit molecule object.
    mol = Chem.MolFromSmiles(smiles)

    # Return false with a reason if the molecule object cannot be created.
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns for the glycerol backbone linked to serine and the phosphatidyl group.
    serine_phosphate_pattern = Chem.MolFromSmarts("OC[C@H](N)C(=O)O")
    phosphatidyl_linkage_pattern = Chem.MolFromSmarts("O[P](=O)(O)[O][C@H]")

    # Check for phosphatidyl group attachment and serine configuration
    if not (mol.HasSubstructMatch(serine_phosphate_pattern) and mol.HasSubstructMatch(phosphatidyl_linkage_pattern)):
        return False, "No phosphatidyl-L-serine structure found"

    # Ensure at least two ester-linked fatty acid chains are present based on the core structure
    fatty_acid_linkage = Chem.MolFromSmarts("C(=O)OC")
    fatty_acid_chains = mol.GetSubstructMatches(fatty_acid_linkage)
    if len(fatty_acid_chains) < 2:
        return False, f"Requires 2 ester-linked fatty acid chains, found {len(fatty_acid_chains)}"

    return True, "Molecule meets the criteria for phosphatidyl-L-serine"