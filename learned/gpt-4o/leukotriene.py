"""
Classifies: CHEBI:25029 leukotriene
"""
from rdkit import Chem

def is_leukotriene(smiles: str):
    """
    Determines if a compound is a leukotriene based on its SMILES string.
    A leukotriene is characterized by a C20 polyunsaturated fatty acid backbone with
    four double bonds, three of which are conjugated.

    Args:
        smiles (str): SMILES string of the compound.

    Returns:
        bool: True if the compound is a leukotriene, False otherwise.
        str: Reason for classification.
    """
    
    # Parse the SMILES string into an RDKit molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return (False, "Invalid SMILES string")
    
    # Use a SMARTS pattern to find at least three conjugated double bonds in a row
    conjugated_pattern = Chem.MolFromSmarts('C=CC=CC=C')
    conjugated_matches = mol.GetSubstructMatches(conjugated_pattern)
    
    if not conjugated_matches:
        return (False, "No pattern of three consecutive conjugated double bonds found.")
    
    # Count total double bonds, ensure there are at least four
    double_bonds = [bond for bond in mol.GetBonds() if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE]
    if len(double_bonds) < 4:
        return (False, f"Insufficient double bonds ({len(double_bonds)}) found; at least 4 required.")
    
    # Evaluate stereochemistry to ensure specific 3D characteristics (often relevant for leukotrienes)
    stereochemistry_patterns = [
        Chem.MolFromSmarts('[C@H](C=C)', useChirality=True),
        Chem.MolFromSmarts('[C@@H](C=C)', useChirality=True)
    ]
    for pattern in stereochemistry_patterns:
        if mol.HasSubstructMatch(pattern):
            return (True, "Structure matches a leukotriene with a C20 polyunsaturated fatty acid backbone and conjugated double bonds.")
    
    return (True, "Structure matches a leukotriene, but requires further stereochemical verification.")

# Example usage:
# result, reason = is_leukotriene('O[C@H](C/C=C\\CCCCC)/C=C/C=C/C=C/C(=O)CCCC(O)=O')
# print(result, reason)