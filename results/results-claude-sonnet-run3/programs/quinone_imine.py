from rdkit import Chem
from rdkit.Chem import AllChem

def is_quinone_imine(smiles: str):
    """
    Determines if a molecule is a quinone imine.
    A quinone imine is formed from a quinone by replacement of one or more quinonoid oxygens by =NH or =NR.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a quinone imine, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # Look for cyclic conjugated system with alternating double bonds
    rings = mol.GetRingInfo()
    if not rings.NumRings():
        return False, "No rings found"

    # Find 6-membered rings
    six_membered_rings = []
    for ring in rings.AtomRings():
        if len(ring) == 6:
            six_membered_rings.append(ring)
            
    if not six_membered_rings:
        return False, "No 6-membered rings found"

    # For each 6-membered ring, check if it has quinone imine pattern
    for ring_atoms in six_membered_rings:
        # Count carbons, nitrogens, oxygens in ring
        c_count = 0
        n_count = 0
        o_count = 0
        double_bonds = 0
        
        # Get ring bonds
        ring_bonds = []
        for i in range(len(ring_atoms)):
            atom1_idx = ring_atoms[i]
            atom2_idx = ring_atoms[(i+1)%6]
            bond = mol.GetBondBetweenAtoms(atom1_idx, atom2_idx)
            if bond:
                ring_bonds.append(bond)
        
        for atom_idx in ring_atoms:
            atom = mol.GetAtomWithIdx(atom_idx)
            symbol = atom.GetSymbol()
            
            if symbol == 'C':
                c_count += 1
            elif symbol == 'N':
                n_count += 1
                # Check if nitrogen has double bond
                if any(bond.GetBondType() == Chem.BondType.DOUBLE for bond in atom.GetBonds()):
                    double_bonds += 1
            elif symbol == 'O':
                o_count += 1
                # Check if oxygen has double bond (carbonyl)
                if any(bond.GetBondType() == Chem.BondType.DOUBLE for bond in atom.GetBonds()):
                    double_bonds += 1

        # Count double bonds in ring
        ring_double_bonds = sum(1 for bond in ring_bonds if bond.GetBondType() == Chem.BondType.DOUBLE)
        
        # Pattern matching for quinone imine:
        # - Must have at least one imine (=N-)
        # - Must have at least one carbonyl (C=O) or another double bond
        # - Must be 6-membered ring
        # - Ring should have conjugated system
        if n_count >= 1 and double_bonds >= 2 and ring_double_bonds >= 2:
            # Check for alternating bond pattern
            has_conjugation = False
            for i in range(len(ring_bonds)):
                if (ring_bonds[i].GetBondType() == Chem.BondType.DOUBLE and 
                    ring_bonds[(i+2)%6].GetBondType() == Chem.BondType.DOUBLE):
                    has_conjugation = True
                    break
            
            if has_conjugation:
                details = f"Ring contains {n_count} nitrogen(s), {o_count} oxygen(s), and {double_bonds} double bond(s)"
                return True, f"Contains quinone imine pattern: {details}"

    return False, "No quinone imine pattern found"
# Pr=None
# Recall=None