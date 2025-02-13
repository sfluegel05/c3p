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
    one nitrogen atom with various degrees of unsaturation.

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

    # Basic patterns for pyrroline and related structures
    patterns = [
        # Basic pyrroline patterns (including saturated and unsaturated forms)
        "[NX3,NX2]1[#6]~[#6]~[#6]~[#6]1",  # Generic pattern
        
        # Lactam forms (pyrrolones)
        "[NX3]1[#6][#6](=[O,S])[#6]~[#6]1",
        "[NX3]1[#6](=[O,S])[#6]~[#6]~[#6]1",
        
        # Imine forms
        "[NX2]1=[#6][#6]~[#6]~[#6]1",
        
        # Charged forms
        "[N+X4]1[#6]~[#6]~[#6]~[#6]1",
        
        # Enamine forms
        "[NX3]1[#6]=[#6][#6]~[#6]1",
        
        # Handle different tautomers
        "[NX3]1[#6]~[#6](=[O,S])[#6]~[#6]1",
    ]

    # Convert patterns to SMARTS objects
    valid_patterns = []
    for pattern in patterns:
        patt = Chem.MolFromSmarts(pattern)
        if patt is not None:
            valid_patterns.append(patt)

    if not valid_patterns:
        return None, "Failed to create valid SMARTS patterns"

    # Check for matches
    for patt in valid_patterns:
        if mol.HasSubstructMatch(patt):
            matches = mol.GetSubstructMatches(patt)
            for match in matches:
                # Get the matched ring atoms
                ring_atoms = [mol.GetAtomWithIdx(i) for i in match]
                
                # Verify it's a 5-membered ring
                ring_info = mol.GetRingInfo()
                rings = ring_info.AtomRings()
                is_in_5ring = False
                for ring in rings:
                    if len(ring) == 5 and all(atom.GetIdx() in ring for atom in ring_atoms):
                        is_in_5ring = True
                        break
                
                if not is_in_5ring:
                    continue
                
                # Count nitrogens in the ring
                n_count = sum(1 for atom in ring_atoms if atom.GetAtomicNum() == 7)
                if n_count != 1:
                    continue
                
                # Count total ring unsaturation (double bonds + formal charges)
                double_bonds = 0
                for i in range(len(match)):
                    atom1_idx = match[i]
                    atom2_idx = match[(i + 1) % 5]
                    bond = mol.GetBondBetweenAtoms(atom1_idx, atom2_idx)
                    if bond.GetBondTypeAsDouble() == 2:
                        double_bonds += 1
                
                formal_charges = sum(abs(atom.GetFormalCharge()) for atom in ring_atoms)
                
                # Check for aromatic ring
                aromatic_atoms = sum(1 for atom in ring_atoms if atom.GetIsAromatic())
                
                # If it's fully aromatic, it's not a pyrroline
                if aromatic_atoms == 5:
                    continue
                    
                # Valid pyrroline should have appropriate unsaturation
                if double_bonds > 0 or formal_charges > 0:
                    return True, "Contains pyrroline ring (5-membered ring with N and appropriate unsaturation)"

    return False, "No pyrroline ring found"