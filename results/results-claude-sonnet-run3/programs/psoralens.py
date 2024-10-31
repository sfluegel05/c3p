from rdkit import Chem
from rdkit.Chem import AllChem

def is_psoralens(smiles: str):
    """
    Determines if a molecule is a psoralen (furanocoumarin with a 7H-furo[3,2-g]chromen-7-one skeleton).
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a psoralen, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS patterns for the core psoralen structure
    # These patterns describe various representations of the fused ring system:
    patterns = [
        'O=C1Oc2cc3occc3cc2C=C1',  # Basic psoralen core
        'O=C1OC2=CC3=C(C=CO3)C=C2C=C1',  # Alternative representation
        'O=C1Oc2cc3c(cc2C=C1)cco3',  # Another representation
        '[O;D2][c;r6]1[c;r6][c;r6][c;r5]2[o;r5][c;r5][c;r5][c;r5]2[c;r6][c;r6]1[C;r6]=[C;r6][C;r6](=O)[O;r6]'  # Detailed pattern
    ]

    # Check if molecule matches any of the core patterns
    matches_core = False
    for pattern in patterns:
        patt = Chem.MolFromSmarts(pattern)
        if patt and mol.HasSubstructMatch(patt):
            matches_core = True
            break

    if not matches_core:
        return False, "Does not contain psoralen core structure"

    # Verify the presence of three fused rings with correct arrangement
    rings = mol.GetRingInfo()
    ring_systems = rings.AtomRings()
    
    if len(ring_systems) < 3:
        return False, "Insufficient number of rings"

    # Count oxygen atoms in ring systems
    ring_atoms = set([atom for ring in ring_systems for atom in ring])
    ring_oxygens = 0
    carbonyl_oxygen = False
    furan_oxygen = False
    
    for atom_idx in ring_atoms:
        atom = mol.GetAtomWithIdx(atom_idx)
        if atom.GetSymbol() == 'O':
            ring_oxygens += 1
            # Check for carbonyl oxygen
            if any(bond.GetBondType() == Chem.BondType.DOUBLE for bond in atom.GetBonds()):
                carbonyl_oxygen = True
            # Check for furan oxygen
            elif len([n for n in atom.GetNeighbors() if n.GetIdx() in ring_atoms]) == 2:
                furan_oxygen = True

    if not (carbonyl_oxygen and furan_oxygen):
        return False, "Missing required oxygen arrangement"

    # Check for substituents
    substituents = []
    for atom in mol.GetAtoms():
        if atom.GetSymbol() not in ['C', 'H', 'O']:
            substituents.append(atom.GetSymbol())
        elif atom.GetSymbol() == 'O' and atom.GetDegree() == 1:  # Hydroxyl groups
            substituents.append('OH')
        elif atom.GetSymbol() == 'O' and atom.GetIsAromatic() == False and atom.GetDegree() == 2:
            # Check if it's part of an ether group outside the core structure
            if not all(n.GetIdx() in ring_atoms for n in atom.GetNeighbors()):
                substituents.append('OR')

    if substituents:
        return True, f"Substituted psoralen with substituents: {', '.join(set(substituents))}"
    return True, "Unsubstituted psoralen"
# Pr=None
# Recall=0.0