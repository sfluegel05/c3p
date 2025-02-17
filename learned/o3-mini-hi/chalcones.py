"""
Classifies: CHEBI:23086 chalcones
"""
#!/usr/bin/env python3
"""
Classifies: Chalcones – a family defined as 1,3-diphenylpropenone (benzylideneacetophenone)
and derivatives. In our implementation a chalcone (or its derivative) is identified by 
either an α,β-unsaturated ketone core (Ar–C(=O)–CH=CH–Ar) or a related dihydrochalcone variant 
(Ar–C(=O)–CH2–CH2–Ar) with both terminal aromatic atoms belonging to benzene rings.
"""

from rdkit import Chem

def is_chalcones(smiles: str):
    """
    Determines if a molecule belongs to the chalcone class based on its SMILES string.
    We look for a core moiety where a carbonyl group (C(=O)) is conjugated to a two‐carbon
    linker which in the unsaturated case is a C=C double bond or in the dihydro case is –CH2–CH2–.
    In both cases, the ketone carbon is attached at one end to an aromatic ring atom and the
    terminal carbon (of the double or single bond sequence) is also part of an aromatic six‐membered
    ring (benzene). Only when such a match is found (and the aromatic ends are indeed in benzene rings)
    is the molecule classified as chalcone (or a chalcone derivative).

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is classified as a chalcone or chalcone derivative, False otherwise.
        str: Reason for the classification decision.
    """
    # Parse the SMILES string into an RDKit molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Helper function: checks whether an atom (by idx) is a member of at least one benzene ring,
    # here defined as a 6-membered ring in which every atom is aromatic carbon.
    def in_benzene_ring(atom_idx, mol):
        ring_info = mol.GetRingInfo()
        for ring in ring_info.AtomRings():
            if len(ring) == 6 and atom_idx in ring:
                # Check that every atom in the ring is aromatic and is a carbon.
                if all(mol.GetAtomWithIdx(idx).GetIsAromatic() and mol.GetAtomWithIdx(idx).GetAtomicNum() == 6 for idx in ring):
                    return True
        return False

    # Define two SMARTS patterns:
    # 1) Unsaturated chalcone: aromatic ring - carbonyl - C=C - aromatic ring.
    unsat_pattern = Chem.MolFromSmarts("[cR]C(=O)[C]=[C][cR]")
    # 2) Dihydrochalcone variant: aromatic ring - carbonyl - CH2 - CH2 - aromatic ring.
    sat_pattern = Chem.MolFromSmarts("[cR]C(=O)CC[cR]")
    
    # Try to find matches for both patterns.
    # Combine the matches in a list.
    matches = []
    if unsat_pattern is not None:
        matches.extend(mol.GetSubstructMatches(unsat_pattern))
    if sat_pattern is not None:
        matches.extend(mol.GetSubstructMatches(sat_pattern))
    
    if not matches:
        return False, "Chalcone core (unsaturated or dihydro variant) not found"
    
    # Check each match to ensure that the first and last atoms (the aromatic ends) belong to benzene rings.
    for match in matches:
        # match[0] is the aromatic atom connected to the carbonyl,
        # match[-1] is the terminal aromatic atom.
        if in_benzene_ring(match[0], mol) and in_benzene_ring(match[-1], mol):
            return True, "Contains chalcone core (α,β-unsaturated or dihydro ketone with aromatic groups on benzene rings)"
    
    return False, "Chalcone core found, but terminal aromatic atoms are not in benzene rings as required"

# For basic testing you can uncomment the lines below:
# test_smiles = "O=C(\\C=C\\c1ccccc1)c1ccccc1"  # trans-chalcone
# result, reason = is_chalcones(test_smiles)
# print(result, reason)