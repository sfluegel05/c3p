from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_tetrachlorobiphenyl(smiles: str):
    """
    Determines if a molecule is a tetrachlorobiphenyl (C12H6Cl4).
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a tetrachlorobiphenyl, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check molecular formula using GetMolecularFormula
    formula = rdMolDescriptors.CalcMolFormula(mol)
    if formula != "C12H6Cl4":
        return False, f"Molecular formula {formula} does not match C12H6Cl4"

    # Check for biphenyl core structure
    # Find aromatic rings
    rings = mol.GetRingInfo()
    aromatic_rings = []
    for ring in rings.AtomRings():
        if len(ring) == 6:
            atoms = [mol.GetAtomWithIdx(i) for i in ring]
            if all(atom.GetIsAromatic() for atom in atoms):
                aromatic_rings.append(ring)
    
    if len(aromatic_rings) != 2:
        return False, "Must contain exactly 2 aromatic rings"

    # Check rings are connected
    ring1_atoms = set(aromatic_rings[0])
    ring2_atoms = set(aromatic_rings[1])
    
    # Look for bonds between rings
    bridge_bonds = 0
    for bond in mol.GetBonds():
        if (bond.GetBeginAtomIdx() in ring1_atoms and bond.GetEndAtomIdx() in ring2_atoms) or \
           (bond.GetBeginAtomIdx() in ring2_atoms and bond.GetEndAtomIdx() in ring1_atoms):
            bridge_bonds += 1
            
    if bridge_bonds != 1:
        return False, "Rings must be connected by exactly one bond"

    return True, "Valid tetrachlorobiphenyl structure"
# Pr=1.0
# Recall=1.0