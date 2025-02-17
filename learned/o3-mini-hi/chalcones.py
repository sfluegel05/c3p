"""
Classifies: CHEBI:23086 chalcones
"""
#!/usr/bin/env python3
"""
Classifies: Chalcones – a family defined as 1,3-diphenylpropenone (benzylideneacetophenone) 
and derivatives, i.e. compounds with an α,β-unsaturated ketone where one end of the 
alkene and the carbonyl group are each connected to an aromatic ring.
"""

from rdkit import Chem

def is_chalcones(smiles: str):
    """
    Determines if a molecule belongs to the chalcone class based on its SMILES string.
    Chalcones must have an α,β-unsaturated ketone where the carbonyl (–C(=O)–) is
    conjugated with a carbon–carbon double bond and both ends of this unsaturated system 
    are attached to aromatic groups (usually benzene rings).

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is classified as a chalcone, False otherwise.
        str: Reason for the classification decision.
    """
    # Parse the SMILES string into an RDKit molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a SMARTS pattern for the chalcone core.
    # This pattern looks for an aromatic carbon (c) attached to a carbonyl carbon (C(=O))
    # which is then attached through a conjugated double bond (/C=C/) to another aromatic carbon (c).
    # This is a minimal representation of the chalcone (1,3-diphenylpropenone) backbone.
    core_pattern = Chem.MolFromSmarts("cC(=O)/C=C/c")
    if not mol.HasSubstructMatch(core_pattern):
        return False, "Chalcone core (α,β-unsaturated carbonyl system with aromatic termini) not found"

    # Count the number of benzene rings (six-membered rings where all atoms are aromatic carbons)
    ring_info = mol.GetRingInfo()
    benzene_ring_count = 0
    for ring in ring_info.AtomRings():
        if len(ring) == 6:
            # Check if all atoms in the ring are aromatic and are carbons (atomic number 6)
            if all(mol.GetAtomWithIdx(idx).GetIsAromatic() and mol.GetAtomWithIdx(idx).GetAtomicNum() == 6 for idx in ring):
                benzene_ring_count += 1

    if benzene_ring_count < 2:
        return False, f"Insufficient benzene rings; expected at least 2 but found {benzene_ring_count}"

    # If both the chalcone core is present and there are at least two benzene rings, classify as chalcone.
    return True, "Contains chalcone core (α,β-unsaturated ketone with aromatic groups at both ends)"

# For basic testing, you can uncomment the following lines:
# test_smiles = "O=C(\\C=C\\c1ccccc1)c1ccccc1"  # trans-chalcone
# result, reason = is_chalcones(test_smiles)
# print(result, reason)