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

    # Make a copy for Kekulization
    mol_kekulized = Chem.Mol(mol)
    try:
        Chem.Kekulize(mol_kekulized, clearAromaticFlags=True)
    except:
        mol_kekulized = None  # If kekulization fails, set to None
    
    # Basic pyrroline patterns
    patterns = [
        # 2,3-dihydro-1H-pyrrole (3,4-double bond)
        "[NX3R5]1[CH2][CH2][CH]=[CH]1",
        # 2,5-dihydro-1H-pyrrole (3,4-double bond)
        "[NX3R5]1[CH2][CH]=[CH][CH2]1",
        # 4,5-dihydro-1H-pyrrole (2,3-double bond)
        "[NX3R5]1[CH]=[CH][CH2][CH2]1",
        
        # More general patterns
        "[NX3R5]1[CR4,CR3][CR4,CR3][CR4,CR3][CR4,CR3]1",  # Any saturated/unsaturated combo
        "[nX3R5]1[cR5,CR4,CR3][cR5,CR4,CR3][cR5,CR4,CR3][cR5,CR4,CR3]1",  # Aromatic version
        
        # Patterns with carbonyls
        "[NX3R5]1[CR4,CR3][C](=[O,S])[CR4,CR3][CR4,CR3]1",
        "[NX3R5]1[C](=[O,S])[CR4,CR3][CR4,CR3][CR4,CR3]1",
        
        # Imine forms
        "[NX2R5]1=[CR4,CR3][CR4,CR3][CR4,CR3][CR4,CR3]1",
        
        # Charged forms
        "[N+X4R5]1[CR4,CR3][CR4,CR3][CR4,CR3][CR4,CR3]1"
    ]

    valid_patterns = []
    for pattern in patterns:
        patt = Chem.MolFromSmarts(pattern)
        if patt is not None:
            valid_patterns.append(patt)
    
    if not valid_patterns:
        return None, "Failed to create valid SMARTS patterns"

    for patt in valid_patterns:
        if mol.HasSubstructMatch(patt):
            # Additional validation
            rings = mol.GetRingInfo().AtomRings()
            for ring in rings:
                if len(ring) == 5:  # Check 5-membered rings
                    ring_atoms = [mol.GetAtomWithIdx(i) for i in ring]
                    n_count = sum(1 for atom in ring_atoms if atom.GetAtomicNum() == 7)
                    if n_count == 1:  # One nitrogen
                        # Count double bonds and aromatic bonds in ring
                        double_bonds = 0
                        aromatic_bonds = 0
                        for i in range(len(ring)):
                            atom1_idx = ring[i]
                            atom2_idx = ring[(i + 1) % 5]
                            bond = mol.GetBondBetweenAtoms(atom1_idx, atom2_idx)
                            if bond.GetBondTypeAsDouble() == 2:
                                double_bonds += 1
                            if bond.GetIsAromatic():
                                aromatic_bonds += 1
                        
                        # Check conditions for pyrroline
                        if aromatic_bonds == 5:  # Aromatic ring is not pyrroline
                            continue
                        if double_bonds == 1 or (double_bonds == 0 and any(a.GetFormalCharge() > 0 for a in ring_atoms)):
                            return True, "Contains pyrroline ring (5-membered ring with N and appropriate unsaturation)"
            
    return False, "No pyrroline ring found"