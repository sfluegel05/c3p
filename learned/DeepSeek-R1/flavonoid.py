"""
Classifies: CHEBI:47916 flavonoid
"""
"""
Classifies: CHEBI:72010 flavonoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_flavonoid(smiles: str):
    """
    Determines if a molecule is a flavonoid based on its SMILES string.
    Flavonoids have a 1-benzopyran (chromen-4-one) core with an aryl substituent at position 2.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a flavonoid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"
    
    # Define corrected chromen-4-one core pattern (benzopyran-4-one)
    # SMARTS: benzene fused to pyrone with ketone at position 4
    core_pattern = Chem.MolFromSmarts('c1ccc2c(c1)oc(=O)cc2')
    if not mol.HasSubstructMatch(core_pattern):
        return False, "No chromen-4-one core found"
    
    # Get all core matches
    matches = mol.GetSubstructMatches(core_pattern)
    
    for match in matches:
        try:
            # Indices in match correspond to SMARTS atoms:
            # c1 (0)-c(1)-c(2)-c(3)-c4(c5)... where c5 is part of the pyrone
            # The pyrone oxygen is at index 7 (assuming SMARTS breakdown)
            # Wait, let's get the correct indices. Let's break down the SMARTS:
            # c1ccc2c(c1)oc(=O)cc2
            # Atom numbering in SMARTS:
            # 0: c1
            # 1: c
            # 2: c
            # 3: c
            # 4: c2 connected back to c1
            # 5: o
            # 6: c(=O)
            # 7: c
            # 8: c2 (end of ring)
            # The C2 position (position adjacent to ketone) is atom 7 in the SMARTS?
            # Wait, the pyrone ring is atoms 4 (c2),5 (o),6 (c=O),7 (c),8 (c), and 0 (c1)?
            # The ketone is at atom 6. The adjacent carbons in the pyrone are 4,5,6,7,8, and 0?
            # The C2 position is the carbon next to the ketone in the pyrone ring.
            # In the SMARTS pattern, the pyrone is arranged as c2(oc(=O)cc2). So the atoms in the pyrone are 4 (c2),5 (o),6 (c=O),7 (c),8 (c), and back to 4.
            # Wait, perhaps the C2 position is atom 7 in the match (the first 'c' after the ketone).
            # Alternatively, the C2 position is the carbon adjacent to the ketone oxygen in the pyrone ring.
            # To find this, locate the ketone (atom 6 in the match) and get its neighboring carbon in the pyrone.
            ketone_atom_idx = match[6]
            ketone_atom = mol.GetAtomWithIdx(ketone_atom_idx)
            
            # Find adjacent atoms in the pyrone ring (should be two carbons)
            pyrone_carbons = []
            for neighbor in ketone_atom.GetNeighbors():
                if neighbor.GetIdx() in match:  # Ensure it's part of the core
                    bond = mol.GetBondBetweenAtoms(ketone_atom_idx, neighbor.GetIdx())
                    if bond.GetBondType() == Chem.BondType.SINGLE:
                        pyrone_carbons.append(neighbor)
            
            # C2 is the carbon adjacent to ketone in the pyrone (position 2)
            if not pyrone_carbons:
                continue
            
            # Check each possible C2 candidate
            for c2_carbon in pyrone_carbons:
                # Look for substituents on C2 that are aryl groups
                for bond in c2_carbon.GetBonds():
                    neighbor = bond.GetOtherAtom(c2_carbon)
                    # Check if substituent is part of an aromatic ring
                    if neighbor.GetIsAromatic():
                        # Verify it's part of a 6-membered aromatic ring
                        ring_info = mol.GetRingInfo()
                        for ring in ring_info.AtomRings():
                            if neighbor.GetIdx() in ring and len(ring) >= 6:
                                return True, "Chromen-4-one core with aryl substituent at position 2"
                        # Alternatively, check if the substituent itself is an aromatic ring
                        if neighbor.IsInRing() and neighbor.GetIsAromatic():
                            return True, "Chromen-4-one core with aryl substituent at position 2"
        except IndexError:
            continue
    
    return False, "No aryl substituent at position 2 of chromen-4-one core"