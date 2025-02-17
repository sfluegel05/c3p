"""
Classifies: CHEBI:132155 hydroxynaphthoquinone
"""
#!/usr/bin/env python
"""
Classifies: hydroxynaphthoquinone
Definition: Any naphthoquinone in which the naphthoquinone moiety is substituted by at least one hydroxy group.
A valid hydroxynaphthoquinone is defined here as a molecule that has:
1. A naphthalene substructure (i.e. two fused aromatic rings, detected via a SMARTS pattern)
2. At least two carbonyl (C=O, ketone) groups attached directly to atoms of that substructure
3. At least one hydroxy group (-OH) attached directly to an atom of that substructure.
Note: This approach is heuristic and may not capture every edge‐case.
"""

from rdkit import Chem

def is_hydroxynaphthoquinone(smiles: str):
    """
    Determines if a molecule is a hydroxynaphthoquinone.
    
    The function looks for a naphthalene substructure (two fused benzene rings)
    and then checks that at least two of its atoms have a carbonyl group (C=O)
    attached and at least one of its atoms has a hydroxy substituent (-OH).
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is recognized as a hydroxynaphthoquinone, else False
        str: Reason explaining the classification decision
    """
    
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a SMARTS pattern for a naphthalene ring system (two fused benzene rings).
    # Note: This pattern is a common representation for naphthalene.
    naphthalene_smarts = "c1ccc2ccccc2c1"
    naph_mol = Chem.MolFromSmarts(naphthalene_smarts)
    
    # Find the naphthalene substructure matches
    matches = mol.GetSubstructMatches(naph_mol)
    if not matches:
        return False, "No naphthalene (fused aromatic) moiety found"
    
    # Iterate over each naphthalene substructure match and check for required carbonyl and hydroxy groups.
    for match in matches:
        # match is a tuple of atom indices that form the naphthalene ring
        ring_atom_indices = set(match)
        carbonyl_count = 0
        hydroxy_count = 0
        
        # Loop over each atom in the found ring substructure
        for idx in ring_atom_indices:
            atom = mol.GetAtomWithIdx(idx)
            # Check all bonds from the ring atom
            for bond in atom.GetBonds():
                # Get the neighboring atom and bond type
                neighbor = bond.GetOtherAtom(atom)
                # Only consider substituents that are not part of the naphthalene ring
                if neighbor.GetIdx() in ring_atom_indices:
                    continue
                
                # Check for a carbonyl group: bond must be double and neighbor an oxygen.
                # (It is assumed that the oxygen in a carbonyl is not further connected to H.)
                if bond.GetBondType() == Chem.BondType.DOUBLE and neighbor.GetAtomicNum() == 8:
                    carbonyl_count += 1
                
                # Check for hydroxy group: single-bonded oxygen that has at least one hydrogen.
                # Note: This heuristic uses the total number of hydrogens on the oxygen.
                if bond.GetBondType() == Chem.BondType.SINGLE and neighbor.GetAtomicNum() == 8:
                    # Using GetTotalNumHs() which includes implicit hydrogens
                    if neighbor.GetTotalNumHs() > 0:
                        hydroxy_count += 1
        
        # For this naphthalene substructure, check if the twin criteria are met
        if carbonyl_count >= 2:
            if hydroxy_count >= 1:
                return True, (f"Found naphthalene ring with {carbonyl_count} carbonyl group(s) and "
                              f"{hydroxy_count} hydroxy substituent(s)")
            else:
                # Substructure found but no hydroxy group attached
                return False, "Naphthoquinone moiety found but lacks a hydroxy substituent"
    # If none of the naphthalene rings meet the criteria
    return False, "Naphthalene ring(s) found but not substituted with both required carbonyl and hydroxy groups"

# For testing purposes:
if __name__ == '__main__':
    test_smiles_list = [
        # lawsone, 2-hydroxy-1,4-naphthoquinone
        "OC1=CC(=O)c2ccccc2C1=O",
        # juglone, another hydroxynaphthoquinone
        "Oc1cccc2C(=O)C=CC(=O)c12",
        # naphthazarin, dihydroxy substituted naphthoquinone 
        "Oc1ccc(O)c2C(=O)C=CC(=O)c12",
        # flaviolin (has hydroxy but only one ketone on ring)
        "Oc1cc(O)c2C(=O)C=C(O)C(=O)c2c1",
        # Not a naphthoquinone 
        "CC(=O)OC1=CC=CC=C1C(=O)O" 
    ]
    for smi in test_smiles_list:
        result, reason = is_hydroxynaphthoquinone(smi)
        print(f"SMILES: {smi}\nResult: {result}, Reason: {reason}\n")