"""
Classifies: CHEBI:134249 alkanesulfonate oxoanion
"""
"""
Classifies: alkanesulfonate oxoanion

Definition: An alkanesulfonate oxoanion is one in which a –S(=O)(=O)[O–] group 
is directly attached to a sp³-hybridized (acyclic) carbon. Furthermore the carbon’s 
immediate neighbors (other than the attached S) should be “alkane‐like” – i.e. not in rings, 
and not sp² (or aromatic). Overall, molecules with a high fraction of aromatic atoms (e.g. dyes)
are treated as false positives.
"""

from rdkit import Chem
from rdkit.Chem import rdchem

def is_alkanesulfonate_oxoanion(smiles: str):
    """
    Determines if a molecule is an alkanesulfonate oxoanion based on its SMILES string.
    
    The molecule must contain a –S(=O)(=O)[O–] group attached directly to a sp³,
    acyclic and non-conjugated carbon. In addition, if the overall molecule is highly aromatic,
    then we suspect it is not a simple alkanesulfonate oxoanion.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as an alkanesulfonate oxoanion, False otherwise.
        str: Explanation detailing the reasoning.
    """
    # Try to parse the SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Global aromatic content check: count heavy atoms and aromatic heavy atoms.
    heavy_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1]
    aromatic_count = sum(1 for atom in heavy_atoms if atom.GetIsAromatic())
    aromatic_ratio = aromatic_count / len(heavy_atoms) if heavy_atoms else 0

    # A molecular aromatic ratio above 0.4 is a red flag (most alkanesulfonates should be largely aliphatic).
    if aromatic_ratio > 0.4:
        global_aromatic = True
    else:
        global_aromatic = False

    # Define a SMARTS pattern for a sulfonate group directly attached to a carbon.
    # Pattern: a non-aromatic carbon ([#6;!a]) bound to a sulfur which in turn is bound to 2 double oxygens and 1 negatively charged oxygen.
    sulfonate_pattern = Chem.MolFromSmarts("[#6;!a]S(=O)(=O)[O-]")
    matches = mol.GetSubstructMatches(sulfonate_pattern)
    if not matches:
        return False, "No sulfonate pattern (C-S(=O)(=O)[O-]) found"

    # Iterate over each match candidate.
    # The SMARTS match returns a tuple of atom indices: (candidate carbon, sulfur, oxygen, oxygen, oxygen)
    for match in matches:
        candidate_idx = match[0]
        sulfur_idx = match[1]
        candidate = mol.GetAtomWithIdx(candidate_idx)
        sulfur = mol.GetAtomWithIdx(sulfur_idx)

        # Check that the candidate (the aliphatic carbon) is sp3 hybridized.
        if candidate.GetHybridization() != rdchem.HybridizationType.SP3:
            continue

        # The candidate carbon must not be in a ring.
        if candidate.IsInRing():
            continue

        # Check that the candidate's other neighbors (ignoring the attached S) are "alkane‐like":
        #   • They should not be in rings.
        #   • They should be sp3 (i.e. not part of a conjugated/unsaturated system).
        valid_environment = True
        carbon_neighbors = []  # count carbon neighbors (aside from S)
        for nbr in candidate.GetNeighbors():
            if nbr.GetIdx() == sulfur_idx:
                continue  # skip the sulfur
            if nbr.IsInRing():
                valid_environment = False
                break
            if nbr.GetHybridization() != rdchem.HybridizationType.SP3:
                valid_environment = False
                break
            if nbr.GetAtomicNum() == 6:
                carbon_neighbors.append(nbr)
        if not valid_environment:
            continue

        # Allow the candidate carbon to have one or two carbon neighbors.
        if len(carbon_neighbors) > 2:
            continue

        # (Optional) If the overall molecule is highly aromatic, then the candidate's immediate neighborhood
        # must not include any aromatic atoms.
        if global_aromatic:
            if any(nbr.GetIsAromatic() for nbr in candidate.GetNeighbors() if nbr.GetIdx() != sulfur_idx):
                continue

        # If a candidate passes all tests, classify it as an alkanesulfonate oxoanion.
        return True, "Contains an alkanesulfonate oxoanion moiety attached to a sp3 (acyclic, non-conjugated) carbon"

    # If no candidate carbon satisfies all criteria:
    return False, "Found sulfonate group, but none attached to a suitable sp3 (acyclic, non-conjugated) carbon"


# (Optional) Main block with tests
if __name__ == "__main__":
    test_smiles = [
        ("C(C(NCCS([O-])(=O)=O)=O)CCCCCCCCCCCCCCCCCC", "N-icosanoyltaurine(1-)"),
        ("CC(=O)NCCCS([O-])(=O)=O", "acamprosate(1-)"),
        ("OCCS([O-])(=O)=O", "isethionate"),
        ("OC[C@H](O)CS([O-])(=O)=O", "(2S)-3-sulfopropanediol(1-)"),
        ("C(CCC)C/C=C\\C/C=C\\C/C=C\\C/C=C\\CCCC(NCCS([O-])(=O)=O)=O", "N-arachidonoyltaurine(1-)"),
        ("[O-]C(=O)C(=O)CS([O-])(=O)=O", "3-sulfonatopyruvate(2-)"),
        ("CCCCS([O-])(=O)=O", "butane-1-sulfonate"),
        # A false positive example: a dye with many aromatic rings
        ("[K+].[K+].[H]C(=CC([H])=CC([H])=C1N(CCCCS([O-])(=O)=O)c2ccc(cc2C1(C)C)C=C([H])C1=[N+](CCCCS([O-])(=O)=O)c2ccc(cc2C1(C)C)S([O-])(=O)=O", "NIR-3 dye"),
    ]
    for sm, name in test_smiles:
        result, explanation = is_alkanesulfonate_oxoanion(sm)
        print(f"SMILES: {sm}\nName: {name}\nResult: {result}, Reason: {explanation}\n")