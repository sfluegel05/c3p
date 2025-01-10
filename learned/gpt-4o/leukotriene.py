"""
Classifies: CHEBI:25029 leukotriene
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_leukotriene(smiles: str):
    """
    Determines if a compound is a leukotriene based on its SMILES string.
    A leukotriene is a C20 polyunsaturated fatty acid with four double bonds, three of which are conjugated.

    Args:
        smiles (str): SMILES string of the compound

    Returns:
        bool: True if the compound is a leukotriene, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES string into an RDKit molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return (False, "Invalid SMILES string")
    
    # Check for exactly 20 carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count != 20:
        return (False, f"Incorrect number of carbon atoms: {c_count} found, 20 required")

    # Find double bonds and check if there are exactly 4
    double_bonds = [bond for bond in mol.GetBonds() if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE]
    if len(double_bonds) != 4:
        return (False, f"Incorrect number of double bonds: {len(double_bonds)} found, 4 required")
    
    # Check for groups of three conjugated double bonds
    conjugated_pattern = Chem.MolFromSmarts('[C]=[C]-[C]=[C]-[C]=[C]')
    conjugated_matches = mol.GetSubstructMatches(conjugated_pattern)
    if len(conjugated_matches) < 1:
        return (False, "No set of three consecutive conjugated double bonds found")

    # Look for hydroxyl (-OH) and carboxyl (-COOH) groups, common in leukotrienes
    has_hydroxyl = mol.HasSubstructMatch(Chem.MolFromSmarts('[OX2H]'))
    has_carboxyl = mol.HasSubstructMatch(Chem.MolFromSmarts('C(=O)[O;H1,-1]'))

    if not (has_hydroxyl and has_carboxyl):
        return (False, "Missing common leukotriene functional groups (hydroxyl or carboxyl)")

    # If all checks pass, return True
    return (True, "Structure matches a leukotriene with C20 fatty acid backbone and three conjugated double bonds.")

# Example usage:
# result, reason = is_leukotriene('O[C@H](C/C=C\\CCCCC)/C=C/C=C/C=C/C(=O)CCCC(O)=O')
# print(result, reason)