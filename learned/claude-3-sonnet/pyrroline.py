"""
Classifies: CHEBI:23763 pyrroline
"""
"""
Classifies: CHEBI:51262 pyrroline
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_pyrroline(smiles: str):
    """
    Determines if a molecule contains a pyrroline ring based on its SMILES string.
    A pyrroline is a dihydropyrrole - a 5-membered heterocyclic ring containing 
    one nitrogen atom and one double bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule contains a pyrroline ring, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Convert to Kekulized form to properly handle aromatic bonds
    Chem.Kekulize(mol, clearAromaticFlags=True)
    
    # SMARTS pattern for pyrroline core:
    # [N] - nitrogen in 5-membered ring
    # [CH2,CH] - sp3 carbons (can be substituted)
    # = indicates double bond
    # Size 5 ring with exactly one double bond
    pyrroline_pattern = Chem.MolFromSmarts("[NX3R5]1[CH2X4,CHX4][CH2X4,CHX4][CH1X3]=[CH1X3]1")
    pattern2 = Chem.MolFromSmarts("[NX3R5]1[CH1X3]=[CH1X3][CH2X4,CHX4][CH2X4,CHX4]1")
    
    if mol.HasSubstructMatch(pyrroline_pattern) or mol.HasSubstructMatch(pattern2):
        # Additional validation to exclude aromatic rings
        rings = mol.GetRingInfo().AtomRings()
        for ring in rings:
            if len(ring) == 5:  # Check 5-membered rings
                ring_atoms = [mol.GetAtomWithIdx(i) for i in ring]
                n_count = sum(1 for atom in ring_atoms if atom.GetAtomicNum() == 7)
                if n_count == 1:  # One nitrogen
                    # Count double bonds in ring
                    double_bonds = 0
                    for i in range(len(ring)):
                        bond = mol.GetBondBetweenAtoms(ring[i], ring[(i+1)%5])
                        if bond.GetBondType() == Chem.BondType.DOUBLE:
                            double_bonds += 1
                    if double_bonds == 1:  # Exactly one double bond
                        return True, "Contains pyrroline ring (5-membered ring with N and one double bond)"
        
        # If we get here, we found the pattern but additional validation failed
        return False, "Found similar pattern but validation failed (might be aromatic or have wrong number of double bonds)"
    
    return False, "No pyrroline ring found"