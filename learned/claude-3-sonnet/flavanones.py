"""
Classifies: CHEBI:28863 flavanones
"""
"""
Classifies: CHEBI:24064 flavanone
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_flavanones(smiles: str):
    """
    Determines if a molecule is a flavanone based on its SMILES string.
    Flavanones have a 3,4-dihydro-2-aryl-2H-1-benzopyran-4-one skeleton.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a flavanone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Core flavanone structure SMARTS pattern:
    # - Benzopyran core with ketone at position 4
    # - No double bond between positions 2 and 3
    # - Aryl group at position 2
    flavanone_pattern = Chem.MolFromSmarts(
        "[#6]1~[#6]~[#6]~[#6]2~[#6](=[O:1])~[#6][#6]~[#6](@[#6]3~[#6]~[#6]~[#6]~[#6]~[#6]3)~[O:2]~2~1"
    )
    
    if not mol.HasSubstructMatch(flavanone_pattern):
        return False, "Missing core flavanone skeleton"

    # Get matches to check specific features
    matches = mol.GetSubstructMatches(flavanone_pattern)
    
    # Check for ketone oxygen and pyran oxygen
    ketone_o = matches[0][0]  # First tagged atom (O:1)
    pyran_o = matches[0][1]   # Second tagged atom (O:2)
    
    # Verify ketone oxygen is indeed a ketone (double bonded)
    ketone_atom = mol.GetAtomWithIdx(ketone_o)
    if not any(bond.GetBondType() == Chem.BondType.DOUBLE 
              for bond in ketone_atom.GetBonds()):
        return False, "Missing ketone group at position 4"

    # Verify pyran oxygen is single bonded (part of the ring)
    pyran_atom = mol.GetAtomWithIdx(pyran_o)
    if not all(bond.GetBondType() == Chem.BondType.SINGLE 
              for bond in pyran_atom.GetBonds()):
        return False, "Incorrect pyran oxygen bonding"

    # Check for absence of double bond between positions 2 and 3
    double_bond_pattern = Chem.MolFromSmarts("O1[CH]=[CH]C(=O)")
    if mol.HasSubstructMatch(double_bond_pattern):
        return False, "Has double bond between positions 2 and 3 (not a flavanone)"

    # Additional check for aryl group at position 2
    aryl_pattern = Chem.MolFromSmarts("O1[CH]([CH2]C(=O))c2ccccc2")
    if not mol.HasSubstructMatch(aryl_pattern):
        return False, "Missing required aryl group at position 2"

    return True, "Contains 3,4-dihydro-2-aryl-2H-1-benzopyran-4-one skeleton"