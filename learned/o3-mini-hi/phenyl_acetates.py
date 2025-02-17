"""
Classifies: CHEBI:140310 phenyl acetates
"""
"""
Classifies: Phenyl Acetates
Definition: An acetate ester obtained by formal condensation of the carboxy group of acetic acid
with the hydroxy group of any phenol.
A molecule qualifies if it contains an acyclic acetoxy group –O–C(=O)CH3 directly attached to an aromatic carbon that is a member of a benzene ring.
"""

from rdkit import Chem

def is_phenyl_acetates(smiles: str):
    """
    Determines whether a molecule is a phenyl acetate based on its SMILES string.
    
    A valid phenyl acetate must contain an acetoxy substructure ([c][O][C](=O)C) where:
        - The oxygen-carbonyl bond is not in a ring (to avoid lactones).
        - The aromatic carbon bearing the group is in a benzene ring (six-membered, all aromatic carbons).
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule qualifies as a phenyl acetate, False otherwise.
        str: Explanation for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the SMARTS pattern for acetoxy group attached to an aromatic carbon.
    acetate_pattern = Chem.MolFromSmarts("[c:1][O:2][C:3](=O)[CH3]")
    if acetate_pattern is None:
        return False, "Error creating SMARTS pattern"
    
    matches = mol.GetSubstructMatches(acetate_pattern)
    if not matches:
        return False, "No phenyl acetate substructure ([c][O][C](=O)C) found."
    
    # Retrieve ring information once
    ring_info = mol.GetRingInfo().AtomRings()
    
    # Process each candidate match.
    for match in matches:
        # Ensure we have exactly four atoms in the match: aromatic C, ester O, carbonyl C, methyl C.
        if len(match) != 4:
            continue
        
        aro_idx, o_idx, carbonyl_idx, methyl_idx = match
        
        # Check that the ester bond (O–C(=O)) is not inside a ring (ruling out lactones).
        in_ring = False
        for ring in ring_info:
            if o_idx in ring and carbonyl_idx in ring:
                in_ring = True
                break
        if in_ring:
            continue
        
        # Verify the aromatic carbon (aro_idx) is aromatic.
        aro_atom = mol.GetAtomWithIdx(aro_idx)
        if not aro_atom.GetIsAromatic():
            continue
        
        # Find a benzene ring: a 6-membered ring of aromatic carbons.
        benzene_found = False
        for ring in ring_info:
            if aro_idx in ring and len(ring) == 6:
                # Check that every atom in the ring is a carbon and is aromatic.
                if all(mol.GetAtomWithIdx(idx).GetAtomicNum() == 6 and mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
                    benzene_found = True
                    break
        if not benzene_found:
            continue
        
        # At this point, we found a valid phenyl acetate substructure.
        return True, "Molecule contains an acyclic phenyl acetate moiety attached to a benzene ring."
    
    return False, "Acetoxy group found but did not meet phenyl acetate criteria (either in a ring or not attached to a benzene ring)."

# Example usage:
if __name__ == "__main__":
    # A known phenyl acetate example
    test_smiles_true = "CC(=O)Oc1ccccc1"  # phenyl acetate
    result, reason = is_phenyl_acetates(test_smiles_true)
    print("Test (phenyl acetate):", result, reason)
    
    # Testing one of the provided examples
    test_smiles_example = "COc1cc(ccc1OC(C)=O)C(C)=O"  # 4-acetyl-2-methoxyphenyl acetate
    result, reason = is_phenyl_acetates(test_smiles_example)
    print("Test (4-acetyl-2-methoxyphenyl acetate):", result, reason)