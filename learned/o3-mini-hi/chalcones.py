"""
Classifies: CHEBI:23086 chalcones
"""
#!/usr/bin/env python3
"""
Classifies: Chalcones – defined as compounds containing a 1,3-diphenylpropenone (benzylideneacetophenone)
core (either unsaturated: Ar-C(=O)-CH=CH-Ar, or the dihydrochalcone variant: Ar-C(=O)-CH2-CH2-Ar)
with both terminal aromatic atoms part of benzene rings.
"""

from rdkit import Chem

def is_chalcones(smiles: str):
    """
    Determines if a molecule belongs to the chalcone class (or its derivatives)
    based on its SMILES string. The core motif is considered valid when: 
       (i) either an α,β-unsaturated ketone (Ar-C(=O)-CH=CH-Ar) or a dihydrochalcone variant 
           (Ar-C(=O)-CH2-CH2-Ar) is present;
       (ii) both terminal aromatic atoms belong to benzene rings (6-membered rings with only aromatic carbons);
       (iii) the ketone carbon and the adjacent linking carbon(s) are not embedded in any ring,
            preventing false positives from fused ring systems.
            
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a chalcone (or chalcone derivative), False otherwise.
        str: Reason for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Helper function: checks whether an atom (by idx) is part of at least one benzene ring.
    def in_benzene_ring(atom_idx, mol):
        ring_info = mol.GetRingInfo()
        for ring in ring_info.AtomRings():
            if len(ring) == 6 and atom_idx in ring:
                # Check that every atom in the ring is aromatic and is a carbon.
                if all(mol.GetAtomWithIdx(idx).GetIsAromatic() and mol.GetAtomWithIdx(idx).GetAtomicNum() == 6 for idx in ring):
                    return True
        return False

    # Define SMARTS patterns:
    # Unsaturated chalcone: aromatic ring (atom0) - carbonyl (atom1) - alkene carbon (atom2) - aromatic ring (atom3)
    unsat_pattern = Chem.MolFromSmarts("[cR]C(=O)[C]=[C][cR]")
    # Dihydrochalcone: aromatic ring (atom0) - carbonyl (atom1) - CH2 (atom2) - CH2 (atom3) - aromatic ring (atom4)
    sat_pattern = Chem.MolFromSmarts("[cR]C(=O)CC[cR]")
    
    valid_match_found = False
    reason = ""
    
    # Attempt both patterns:
    all_matches = []
    if unsat_pattern is not None:
        all_matches.extend(mol.GetSubstructMatches(unsat_pattern))
    if sat_pattern is not None:
        all_matches.extend(mol.GetSubstructMatches(sat_pattern))
    
    if not all_matches:
        return False, "Chalcone core (unsaturated or dihydro variant) not found"
    
    for match in all_matches:
        # For the unsaturated pattern we expect a 4-atom match and for the saturated (dihydro) we expect 5.
        if len(match) == 4:
            # For unsaturated chalcone:
            aromatic1 = match[0]
            carbonyl = match[1]
            linker = match[2]
            aromatic2 = match[3]
            
            # Check that both terminal atoms are in benzene rings.
            if not (in_benzene_ring(aromatic1, mol) and in_benzene_ring(aromatic2, mol)):
                continue
            
            # Check that the carbonyl and its adjacent alkene carbon are not in any ring;
            # in chalcones the carbonyl should be in an open chain.
            if mol.GetAtomWithIdx(carbonyl).IsInRing() or mol.GetAtomWithIdx(linker).IsInRing():
                continue
            
            # Additional check: ensure the carbonyl atom is not aromatic.
            if mol.GetAtomWithIdx(carbonyl).GetIsAromatic():
                continue
            
            valid_match_found = True
            reason = "Contains chalcone core (α,β-unsaturated ketone with aromatic groups on benzene rings)"
            break
            
        elif len(match) == 5:
            # For dihydrochalcone variant:
            aromatic1 = match[0]
            carbonyl = match[1]
            linker1 = match[2]
            linker2 = match[3]
            aromatic2 = match[4]
            
            if not (in_benzene_ring(aromatic1, mol) and in_benzene_ring(aromatic2, mol)):
                continue
            
            if mol.GetAtomWithIdx(carbonyl).IsInRing() or mol.GetAtomWithIdx(linker1).IsInRing() or mol.GetAtomWithIdx(linker2).IsInRing():
                continue
            
            if mol.GetAtomWithIdx(carbonyl).GetIsAromatic():
                continue
            
            valid_match_found = True
            reason = "Contains chalcone core (dihydroketone variant with aromatic groups on benzene rings)"
            break
    
    if valid_match_found:
        return True, reason
    else:
        return False, "Chalcone core found, but connectivity constraints (non-cyclic carbonyl/linker or benzene terminal rings) not met"

# For basic testing you can uncomment the lines below:
# test_smiles = "O=C(\\C=C\\c1ccccc1)c1ccccc1"  # trans-chalcone
# result, reason = is_chalcones(test_smiles)
# print(result, reason)