"""
Classifies: CHEBI:134249 alkanesulfonate oxoanion
"""
"""
Classifies: alkanesulfonate oxoanion

Definition: An alkanesulfonate oxoanion is one in which a –S(=O)(=O)[O–] group
is directly attached to a non‐ring, non‐aromatic carbon (position 1). 
The carbon can be saturated (sp3) or, if unsaturated (sp2), its non-alkyl neighbors are 
allowed if they are carbonyl groups. Furthermore, the candidate’s other neighbors (aside from S)
should be “alkane‐like” – e.g. not in rings and typically sp3 (or, if sp2, part of a carbonyl).
Overall, molecules with a high fraction of aromatic atoms (e.g. dyes) are treated as negatives.
"""

from rdkit import Chem
from rdkit.Chem import rdchem

def is_alkanesulfonate_oxoanion(smiles: str):
    """
    Determines if a molecule is an alkanesulfonate oxoanion based on its SMILES string.

    The molecule must contain a –S(=O)(=O)[O–] group attached directly to a carbon
    that is non-aromatic and acyclic. The candidate carbon may be sp3, or sp2 if the modification
    is due to an adjacent carbonyl group. In addition, if the overall molecule is highly aromatic,
    then the match is rejected (to weed out dyes and highly conjugated structures).

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is classified as an alkanesulfonate oxoanion, False otherwise.
        str: Explanation detailing the reasoning.
    """
    # Parse the SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Global aromatic content check:
    heavy_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1]
    aromatic_count = sum(1 for atom in heavy_atoms if atom.GetIsAromatic())
    aromatic_ratio = aromatic_count / len(heavy_atoms) if heavy_atoms else 0

    if aromatic_ratio > 0.4:
        return False, "Molecule has high aromatic content, suggesting it is not a simple alkanesulfonate oxoanion"

    # Define a SMARTS pattern for a sulfonate group directly attached to a non-aromatic carbon.
    # We require that the carbon is not aromatic and is not in a ring.
    sulfonate_pattern = Chem.MolFromSmarts("[C;!a]S(=O)(=O)[O-]")
    matches = mol.GetSubstructMatches(sulfonate_pattern)
    if not matches:
        return False, "No sulfonate pattern (non-aromatic C-S(=O)(=O)[O-]) found in molecule"

    # Iterate over candidate matches.
    for match in matches:
        # match tuple: (candidate carbon, sulfur, oxygen, oxygen, oxygen)
        candidate_idx = match[0]
        sulfur_idx = match[1]
        candidate = mol.GetAtomWithIdx(candidate_idx)
        sulfur = mol.GetAtomWithIdx(sulfur_idx)

        # Ensure candidate carbon is not in a ring.
        if candidate.IsInRing():
            continue

        # Allow candidate carbon that is sp3.
        candidate_ok = False
        if candidate.GetHybridization() == rdchem.HybridizationType.SP3:
            candidate_ok = True
        # Or allow candidate that is sp2 if it appears to be part of a carbonyl functionality.
        elif candidate.GetHybridization() == rdchem.HybridizationType.SP2:
            # Check if candidate has a double-bond to oxygen.
            for bond in candidate.GetBonds():
                if bond.GetBondType() == Chem.BondType.DOUBLE:
                    other = bond.GetOtherAtom(candidate)
                    if other.GetAtomicNum() == 8:
                        candidate_ok = True
                        break
        if not candidate_ok:
            continue

        # Check neighbors of the candidate (except the attached sulfur)
        valid_environment = True
        carbon_neighbors = 0
        for nbr in candidate.GetNeighbors():
            if nbr.GetIdx() == sulfur_idx:
                continue
            # Must not be in a ring.
            if nbr.IsInRing():
                valid_environment = False
                break
            # For carbon neighbors:
            if nbr.GetAtomicNum() == 6:
                carbon_neighbors += 1
                # If neighbor is sp3, that is acceptable.
                if nbr.GetHybridization() == rdchem.HybridizationType.SP3:
                    continue
                # If neighbor is sp2, allow only if it is part of a carbonyl group.
                elif nbr.GetHybridization() == rdchem.HybridizationType.SP2:
                    has_carbonyl = False
                    for bond in nbr.GetBonds():
                        if bond.GetBondType() == Chem.BondType.DOUBLE:
                            other = bond.GetOtherAtom(nbr)
                            if other.GetAtomicNum() == 8:
                                has_carbonyl = True
                                break
                    if not has_carbonyl:
                        valid_environment = False
                        break
                else:
                    valid_environment = False
                    break
            else:
                # For other heteroatoms, no special restrictions.
                continue

        # Allow candidate with one or two carbon neighbors.
        if carbon_neighbors > 2:
            valid_environment = False

        if not valid_environment:
            continue

        # If this candidate passes all tests then classify as a positive.
        return True, "Contains an alkanesulfonate oxoanion moiety attached to a suitable aliphatic (acyclic, non-conjugated) carbon"

    # If no candidate carbon satisfies all criteria:
    return False, "Found sulfonate group, but none attached to a suitable acyclic non-aromatic carbon environment"

# (Optional) Main block with tests
if __name__ == "__main__":
    test_cases = [
        ("C(C(NCCS([O-])(=O)=O)=O)CCCCCCCCCCCCCCCCCC", "N-icosanoyltaurine(1-)"),
        ("CC(=O)NCCCS([O-])(=O)=O", "acamprosate(1-)"),
        ("OCCS([O-])(=O)=O", "isethionate"),
        ("OC[C@H](O)CS([O-])(=O)=O", "(2S)-3-sulfopropanediol(1-)"),
        ("[O-]S(C[C@H]([C@@H](O)C=O)O)(=O)=O", "4-deoxy-4-sulfo-D-erythrose(1-)"),
        ("C(CCC)C/C=C\\C/C=C\\C/C=C\\C/C=C\\CCCC(NCCS([O-])(=O)=O)=O", "N-arachidonoyltaurine(1-)"),
        ("CCCCS([O-])(=O)=O", "butane-1-sulfonate"),
        ("C(C(NCCS([O-])(=O)=O)=O)CCCCCCCCCCCCCCCCCCCCC", "N-tricosanoyltaurine(1-)"),
        ("O[C@H](CS([O-])(=O)=O)C([O-])=O", "(S)-3-sulfonatolactate(2-)"),
        ("C(C(NCCS([O-])(=O)=O)=O)CCCCCCCCCCCC", "N-tetradecanoyltaurine(1-)"),
        ("C(C(NCCS([O-])(=O)=O)=O)CCCCCCCCCCCCCCCCCCCCCC", "N-tetracosanoyltaurine(1-)"),
        ("OC[C@@H](O)CS([O-])(=O)=O", "(2R)-3-sulfopropanediol(1-)"),
        ("C(=C/[C@H](C/C=C\\CCCCC)OO)\\C=C/C/C=C\\CCCC(=O)NCCS([O-])(=O)=O", "N-[12(S)-hydroperoxy-(5Z,8Z,10E,14Z)-icosatetraenoyl]taurine(1-)"),
        ("[H][N+](CCO)(CCO)CCS([O-])(=O)=O", "2-[bis(2-hydroxyethyl)ammonio]ethanesulfonate"),
        ("C(C(NCCS([O-])(=O)=O)=O)CCCCCCCCCC", "N-dodecanoyltaurine(1-)"),
        ("[O-]S(C[C@@H](C(=O)[H])O)(=O)=O", "D-3-sulfolactaldehyde(1-)"),
        ("[H][N+]([H])(CCS([O-])(=O)=O)CC(N)=O", "2-[(2-amino-2-oxoethyl)ammonio]ethanesulfonate"),
        ("[O-]S(C[C@H](C(=O)[H])O)(=O)=O", "L-3-sulfolactaldehyde(1-)"),
        ("NCCS([O-])(=O)=O", "2-aminoethanesulfonate"),
        ("OCCN(CCO)CCS([O-])(=O)=O", "2-[bis(2-hydroxyethyl)amino]ethanesulfonate"),
        ("[O-]S(C[C@H](C(=O)CO)O)(=O)=O", "4-deoxy-4-sulfo-D-erythrulose(1-)"),
        ("OCC(O)CS([O-])(=O)=O", "3-sulfopropanediol(1-)"),
        ("C(CCC)C[C@@H](/C=C/C=C\\C/C=C\\C/C=C\\CCCC(NCCS([O-])(=O)=O)=O)OO", "N-[15(S)-hydroperoxy-(5Z,8Z,11Z,13E)-icosatetraenoyl]taurine(1-)"),
        ("[Na+].[H]C(=C([H])C([H])=C([H])C1=[N+](CCCCS([O-])(=O)=O)c2ccc3ccccc3c2C1(C)C)C([H])=C([H])C([H])=C1N(CCCCS([O-])(=O)=O)c2ccc3ccccc3c2C1(C)C", "indocyanine green"),
        ("CCCCCCCCCCCCCCCS([O-])(=O)=O", "hexadecane-1-sulfonate"),
        ("C(CS([O-])(=O)=O)NC(C)=O", "acetyltaurine(1-)"),
        ("FC(S([O-])(=O)=O)(F)F", "triflate"),
        ("OC(CS([O-])(=O)=O)C=O", "3-sulfolactaldehyde(1-)"),
        ("CC(C)S([O-])(=O)=O", "propane-2-sulfonate"),
        ("[O-]C(=O)CS([O-])(=O)=O", "sulfonatoacetate"),
        ("CS([O-])(=O)=O", "methanesulfonate"),
        # False positive test examples (should be rejected)
        ("[K+].[K+].[H]C(=CC([H])=CC([H])=C1N(CCCCS([O-])(=O)=O)c2ccc(cc2C1(C)C)C(O)=O)C=C([H])C1=[N+](CCCCS([O-])(=O)=O)c2ccc(cc2C1(C)C)S([O-])(=O)=O", "NIR-3 dye"),
    ]
    for sm, name in test_cases:
        res, expl = is_alkanesulfonate_oxoanion(sm)
        print(f"SMILES: {sm}\nName: {name}\nResult: {res}, Reason: {expl}\n")