"""
Classifies: CHEBI:36092 clavulone
"""
"""
Classifies: Clavulone – a class of esterified prostanoids obtained from marine corals.
This heuristic classification searches for a cyclopentenone-like (5-membered ring containing a carbonyl)
core and at least one ester (–OC(=O)–R) group.
"""

from rdkit import Chem

def is_clavulone(smiles: str):
    """
    Determines if a molecule is a clavulone based on its SMILES string.
    Heuristic criteria:
      1. The molecule contains a five-membered ring having at least one carbonyl group 
         (a rough approximation of a cyclopentenone core).
      2. At least one ester group (–OC(=O)–R) is present in the structure.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is classified as clavulone, False otherwise.
        str : Reason for the classification decision.
    """
    
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Criterion 1: Look for a five-membered ring (cyclopentane core) that contains a carbonyl.
    ring_info = mol.GetRingInfo()
    cyclopentenone_found = False
    for ring in ring_info.AtomRings():
        if len(ring) == 5:
            # Check if the ring has at least one carbon (atomic num 6) with a double-bonded oxygen.
            for idx in ring:
                atom = mol.GetAtomWithIdx(idx)
                if atom.GetAtomicNum() == 6:  # carbon atom
                    for neighbor in atom.GetNeighbors():
                        if neighbor.GetAtomicNum() == 8:  # oxygen atom
                            bond = mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx())
                            if bond is not None and bond.GetBondType() == Chem.BondType.DOUBLE:
                                cyclopentenone_found = True
                                break
                    if cyclopentenone_found:
                        break
            if cyclopentenone_found:
                break
    if not cyclopentenone_found:
        return False, "No cyclopentenone ring (5-membered ring with a carbonyl group) found"
    
    # Criterion 2: Look for an ester group.
    # The SMARTS "[#6][O][C](=O)[#6]" covers an ester group with a carbon attached to the O.
    ester_pattern = Chem.MolFromSmarts("[#6][O][C](=O)[#6]")
    if not mol.HasSubstructMatch(ester_pattern):
        return False, "No ester functional group (-OC(=O)-R) found in the molecule"
    
    # If both criteria are met, we classify the molecule as belonging to the clavulone class.
    return True, "Contains a cyclopentenone ring and at least one ester group – consistent with clavulone"

# Example usage:
if __name__ == "__main__":
    # Example: punaglandin 2
    smiles_example = "ClC=1C(=O)[C@@]([C@@](O)(C/C=C\\CCCCC)C1)([C@@H](OC(=O)C)[C@H](OC(=O)C)[C@@H](OC(=O)C)CCCC(OC)=O)[H]"
    result, reason = is_clavulone(smiles_example)
    print(result, reason)