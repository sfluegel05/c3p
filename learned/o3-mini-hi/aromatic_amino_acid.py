"""
Classifies: CHEBI:33856 aromatic amino acid
"""
"""
Classifies: aromatic amino acid
Definition: An amino acid whose structure includes an aromatic ring.
"""
from rdkit import Chem

def is_aromatic_amino_acid(smiles: str):
    """
    Determines if a molecule is an aromatic amino acid based on its SMILES string.
    The molecule must contain a carboxylic acid group, a free (non-amidic) amino group,
    and at least one aromatic ring.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule meets criteria for an aromatic amino acid, False otherwise.
        str: Explanation of the classification result.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add hydrogens so that we can count explicit hydrogens on nitrogen atoms.
    mol = Chem.AddHs(mol)
    
    # 1. Check for a carboxylic acid group.
    # We look for either the protonated form (-C(=O)OH) or the deprotonated carboxylate (-C(=O)[O-]).
    acid_smarts_1 = Chem.MolFromSmarts("[CX3](=O)[OX2H]")   # -C(=O)OH
    acid_smarts_2 = Chem.MolFromSmarts("[CX3](=O)[OX1-]")    # -C(=O)[O-]
    has_acid = mol.HasSubstructMatch(acid_smarts_1) or mol.HasSubstructMatch(acid_smarts_2)
    if not has_acid:
        return False, "Missing carboxylic acid group"
    
    # 2. Check for a free amino group.
    # We define a free amino group as a nitrogen atom (atomic number 7) that
    # has at least one hydrogen attached and is not directly bonded to a carbonyl carbon.
    free_amine_found = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 7:  # nitrogen
            # Count the number of attached hydrogens (after adding explicit hydrogens)
            n_h = sum(1 for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 1)
            if n_h < 1:
                continue  # not a free amine if no H's are attached
            is_amide = False
            # Check if the nitrogen is directly bonded to a carbon that has a double bonded oxygen.
            for nbr in atom.GetNeighbors():
                if nbr.GetAtomicNum() == 6:  # carbon
                    for nbr2 in nbr.GetNeighbors():
                        if nbr2.GetAtomicNum() == 8:
                            bond = mol.GetBondBetweenAtoms(nbr.GetIdx(), nbr2.GetIdx())
                            if bond and bond.GetBondTypeAsDouble() == 2.0:
                                is_amide = True
                                break
                    if is_amide:
                        break
            if not is_amide:
                free_amine_found = True
                break
    if not free_amine_found:
        return False, "Missing free amino group (non-amidic nitrogen with attached hydrogens)"
    
    # 3. Check for an aromatic ring.
    aromatic_ring_found = False
    ri = mol.GetRingInfo()
    for ring in ri.AtomRings():
        # ring is a tuple of atom indices
        if len(ring) < 5:
            # usually aromatic rings are 5- or 6-membered (or more)
            continue
        # Check if every atom in the ring is marked as aromatic.
        if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            aromatic_ring_found = True
            break
    if not aromatic_ring_found:
        return False, "No aromatic ring found"
    
    return True, "Contains carboxylic acid group, free amino group, and aromatic ring"

# Example usage (try any of the provided SMILES strings):
if __name__ == "__main__":
    test_smiles = "NC1=C(C=C(Cl)C=C1)C(O)=O"  # 2-amino-5-chlorobenzoic acid
    result, reason = is_aromatic_amino_acid(test_smiles)
    print(result, reason)