from rdkit import Chem
from rdkit.Chem import AllChem

def is_epsilon_lactone(smiles: str):
    """
    Determines if a molecule contains an epsilon-lactone (7-membered lactone ring).
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule contains an epsilon-lactone, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # SMARTS patterns for 7-membered lactone rings
    lactone_patterns = [
        # Basic 7-membered lactone pattern
        "[O;R1]1-[#6;R1]-[#6;R1]-[#6;R1]-[#6;R1]-[#6;R1]-[#6;R1](=[O;!R])-1",
        # Alternative pattern with different carbonyl position
        "[O;R1]1-[#6;R1]-[#6;R1]-[#6;R1]-[#6;R1]-[#6;R1](=[O;!R])-[#6;R1]-1",
        # Pattern for bridged systems
        "[O;R1]1-[#6;R2]-[#6;R2]-[#6;R1]-[#6;R1]-[#6;R1]-[#6;R1](=[O;!R])-1"
    ]

    for pattern in lactone_patterns:
        pat = Chem.MolFromSmarts(pattern)
        if pat and mol.HasSubstructMatch(pat):
            return True, "Contains a 7-membered lactone ring"

    # Additional check for more complex 7-membered lactone rings
    rings = mol.GetRingInfo()
    for ring in rings.AtomRings():
        if len(ring) != 7:
            continue

        # Check if ring contains oxygen
        ring_atoms = set(ring)
        oxygen_in_ring = False
        carbonyl_carbon = None
        carbonyl_oxygen = None

        for atom_idx in ring:
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetSymbol() == 'O':
                oxygen_in_ring = True
                # Check if this oxygen is connected to a carbon with a carbonyl
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetIdx() in ring_atoms:
                        for n2 in neighbor.GetNeighbors():
                            if (n2.GetSymbol() == 'O' and 
                                n2.GetIdx() not in ring_atoms and 
                                mol.GetBondBetweenAtoms(neighbor.GetIdx(), n2.GetIdx()).GetBondType() == Chem.BondType.DOUBLE):
                                carbonyl_carbon = neighbor.GetIdx()
                                carbonyl_oxygen = n2.GetIdx()

        if oxygen_in_ring and carbonyl_carbon is not None:
            # Verify the connectivity forms a lactone
            for atom_idx in ring:
                if mol.GetAtomWithIdx(atom_idx).GetSymbol() == 'O':
                    # Check if there's a path from the ring oxygen to the carbonyl carbon
                    path = Chem.FindAllPathsOfLengthN(mol, 6, atom_idx, carbonyl_carbon)
                    if path:
                        return True, "Contains a 7-membered lactone ring"

    return False, "No lactone functionality found in 7-membered rings"
# Pr=0.6
# Recall=1.0