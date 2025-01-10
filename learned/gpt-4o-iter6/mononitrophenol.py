"""
Classifies: CHEBI:39362 mononitrophenol
"""
from rdkit import Chem

def is_mononitrophenol(smiles: str):
    """
    Determines if a molecule is a mononitrophenol based on its SMILES string.
    A mononitrophenol is defined as a phenolic compound where the aromatic
    phenol ring carries a single nitro substituent at an unspecified position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a mononitrophenol, False otherwise
        tuple: Reason for classification or potential compound details
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"
    
    # Define flexible phenol pattern with allowance for additional substituents
    phenol_pattern = Chem.MolFromSmarts("c1ccc(O)c([cH,$(c[N+](=O)[O-])])c1")  
    # Define strict nitro group in aromatic context
    nitro_aromatic_pattern = Chem.MolFromSmarts("a-[N+](=O)[O-]")
    
    # Check for presence of phenol structure
    phenol_matches = mol.GetSubstructMatches(phenol_pattern)
    if len(phenol_matches) == 0:
        return False, "No phenolic structure with aromatic characteristics found"
    
    # Check for aromatic nitro group
    nitro_matches = mol.GetSubstructMatches(nitro_aromatic_pattern)
    if len(nitro_matches) != 1:
        return False, f"Found {len(nitro_matches)} unique aromatic-attached nitro groups, need exactly 1"

    # Verify nitro in part of aromatic, tied to phenolic context
    for match in phenol_matches:
        # Exclude OH, capturing only the phenol core as reference
        phenolic_aromatic_carbon = match[:-1]  
        nitro_carbon = nitro_matches[0]  # Using index of nitro found

        # Checking nitro linkage within valid phenol context
        for bond in mol.GetBonds():
            begin_idx = bond.GetBeginAtomIdx()
            end_idx = bond.GetEndAtomIdx()
            if ({begin_idx, end_idx} & set(phenolic_aromatic_carbon)) and \
               (begin_idx in nitro_carbon or end_idx in nitro_carbon):
                return True, "Molecule is a mononitrophenol with nitro group functionally tied to phenol casing"
    
    return False, "The nitro is not integrated correctly in context with the phenol"