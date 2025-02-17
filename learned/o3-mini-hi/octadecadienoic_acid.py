"""
Classifies: CHEBI:25627 octadecadienoic acid
"""
"""
Classifies: octadecadienoic acid – a straight-chain C18 polyunsaturated fatty acid having exactly 2 C=C double bonds.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_octadecadienoic_acid(smiles: str):
    """
    Determines if a molecule is an octadecadienoic acid based on its SMILES string.
    An octadecadienoic acid is defined as any straight-chain, C18 polyunsaturated fatty acid with exactly two carbon-carbon double bonds,
    and must have a terminal carboxylic acid group.
    
    Steps of validation:
    1. Parse the SMILES string.
    2. Check for a carboxylic acid substructure.
    3. Check that the total number of carbon atoms is 18.
    4. Check that there are exactly 2 C=C (non‐aromatic) double bonds.
    5. Check that the carbon skeleton is “straight-chain” (i.e. no carbon is connected to more than two other carbons).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is octadecadienoic acid, False otherwise.
        str: A message providing the reason for classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # 1. Verify presence of carboxylic acid group.
    # We look for a pattern for a typical carboxylic acid: carbonyl C(=O) bonded to an OH.
    acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    if not mol.HasSubstructMatch(acid_pattern):
        return False, "Missing carboxylic acid functionality"
    
    # 2. Check that the molecule has exactly 18 carbon atoms (atomic number 6).
    c_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    if len(c_atoms) != 18:
        return False, f"Expected 18 carbon atoms but found {len(c_atoms)}"
    
    # 3. Count the number of carbon-carbon double bonds. Only count non-aromatic ones.
    double_bonds = 0
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
            a1 = bond.GetBeginAtom()
            a2 = bond.GetEndAtom()
            if (a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 6) and (not bond.GetIsAromatic()):
                double_bonds += 1
    if double_bonds != 2:
        return False, f"Found {double_bonds} carbon-carbon double bonds; exactly 2 required"
    
    # 4. Check for “straight-chain” via the carbon connectivity.
    # In a straight chain, each carbon (except the chain endpoints) is connected to exactly 2 other carbons.
    # We check all carbon atoms in the molecule.
    for atom in c_atoms:
        # Count number of neighbouring atoms that are carbon.
        n_c_neighbors = sum(1 for neighbor in atom.GetNeighbors() if neighbor.GetAtomicNum() == 6)
        # End carbons in a linear chain should have one carbon neighbor; internal should have exactly 2.
        if n_c_neighbors > 2:
            return False, f"Carbon atom with index {atom.GetIdx()} has {n_c_neighbors} carbon neighbors, indicating branching"
    
    return True, "Molecule is a straight-chain C18 fatty acid with exactly 2 C=C double bonds and a carboxylic acid group"

# Examples for testing:
if __name__ == "__main__":
    test_smiles = [
        "OC(=O)CCCCCCC\\C=C/C=C\\CCCCCC",           # 9Z,11Z-octadecadienoic acid
        "CCCCCC\\C=C(/C\\C=C/CCCCCCCC(O)=O)[N+]([O-])=O",  # 12-Nitro-9Z,12Z-octadecadienoic acid (likely to fail due to extra branching)
        "CCCCC\\C=C\\C\\C=C\\CCCCCCCC(O)=O"         # linoelaidic acid
    ]
    for sm in test_smiles:
        result, reason = is_octadecadienoic_acid(sm)
        print(f"SMILES: {sm}\nResult: {result}, Reason: {reason}\n")