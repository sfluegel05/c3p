"""
Classifies: CHEBI:25029 leukotriene
"""
from rdkit import Chem
from rdkit.Chem import AllChem

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
    
    # Check for at least 20 carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 20:
        return (False, f"Insufficient carbon atoms: {c_count} found, 20 required")

    # Find double bonds
    double_bonds = [bond for bond in mol.GetBonds() if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE]
    if len(double_bonds) != 4:
        return (False, f"Incorrect number of double bonds: {len(double_bonds)} found, 4 required")
    
    # Check for conjugation of double bonds
    conjugated_double_bonds = 0
    for bond in double_bonds:
        # Get atoms involved in the double bond
        start_atom = bond.GetBeginAtom()
        end_atom = bond.GetEndAtom()

        # Check neighboring bonds
        for neighbor in start_atom.GetBonds():
            if neighbor.GetBondType() == Chem.rdchem.BondType.SINGLE:
                # Check if there's a double bond on the other side of the single bond
                for next_neighbor in neighbor.GetOtherAtom(start_atom).GetBonds():
                    if next_neighbor.GetBondType() == Chem.rdchem.BondType.DOUBLE and next_neighbor != bond:
                        conjugated_double_bonds += 1

    if conjugated_double_bonds < 3:
        return (False, f"Insufficient conjugated double bonds: {conjugated_double_bonds} found, 3 required")
    
    # Check for common functional groups
    # We expect to find features such as hydroxyl groups for leukotrienes
    hydroxyl_group = Chem.MolFromSmarts('[OX2H]')
    if not mol.HasSubstructMatch(hydroxyl_group):
        return (False, "Missing hydroxyl group, a common feature of leukotrienes")
    
    # If all checks pass, return True
    return (True, "Structure matches a leukotriene with C20 fatty acid backbone and conjugated double bonds.")

# Example usage:
# result, reason = is_leukotriene('O[C@H](C/C=C\\CCCCC)/C=C/C=C/C=C/C(=O)CCCC(O)=O')
# print(result, reason)