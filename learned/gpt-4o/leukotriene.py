"""
Classifies: CHEBI:25029 leukotriene
"""
from rdkit import Chem

def is_leukotriene(smiles: str):
    """
    Determines if a compound is a leukotriene based on its SMILES string.
    A leukotriene is defined by a C20 polyunsaturated fatty acid backbone with
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
    
    # Count the carbon atoms allowing for slight variation (18 to 22)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if not (18 <= c_count <= 22):
        return (False, f"Number of carbon atoms ({c_count}) outside permitted range (18-22) for leukotrienes.")

    # Find all double bonds
    double_bonds = [bond for bond in mol.GetBonds() if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE]
    
    # Ensure there are at least four double bonds
    if len(double_bonds) < 4:
        return (False, f"Insufficient double bonds ({len(double_bonds)}) found; at least 4 required.")
    
    # Check for at least one set of three consecutive conjugated double bonds
    conjugated_pattern = Chem.MolFromSmarts('C=C-C=C-C=C')
    conjugated_matches = mol.GetSubstructMatches(conjugated_pattern)
    
    if len(conjugated_matches) < 1:
        return (False, "No set of three consecutive conjugated double bonds found.")
    
    # If all checks pass, the compound can be considered a leukotriene
    return (True, "Structure matches a leukotriene with a C20 polyunsaturated fatty acid backbone and conjugated double bonds.")

# Example usage:
# result, reason = is_leukotriene('O[C@H](C/C=C\\CCCCC)/C=C/C=C/C=C/C(=O)CCCC(O)=O')
# print(result, reason)