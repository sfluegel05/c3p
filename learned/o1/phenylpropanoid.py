"""
Classifies: CHEBI:26004 phenylpropanoid
"""
"""
Classifies: Phenylpropanoid
"""

from rdkit import Chem

def is_phenylpropanoid(smiles: str):
    """
    Determines if a molecule is a phenylpropanoid based on its SMILES string.
    A phenylpropanoid is any organic aromatic compound with a structure based on a phenylpropane skeleton.
    This includes subclasses such as flavonoids, coumarins, lignans, and anthocyanins.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phenylpropanoid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # List of SMARTS patterns for different phenylpropanoid subclasses

    # Phenylpropane backbone: benzene ring connected to a three-carbon chain
    phenylpropane_pattern = Chem.MolFromSmarts("c1ccccc1CC[CH2]")  # Allow for primary and secondary carbons
    if mol.HasSubstructMatch(phenylpropane_pattern):
        return True, "Contains phenylpropane skeleton"

    # Cinnamic acid derivatives: phenyl group connected to propenoic acid
    cinnamic_acid_pattern = Chem.MolFromSmarts("c1ccccc1C=CC(=O)O")
    if mol.HasSubstructMatch(cinnamic_acid_pattern):
        return True, "Contains cinnamic acid structure"

    # Coumarins: benzopyrone structure
    coumarin_pattern = Chem.MolFromSmarts("O=C1C=CC2=CC=CC=C2O1")
    if mol.HasSubstructMatch(coumarin_pattern):
        return True, "Contains coumarin structure"

    # Flavonoids: 15-carbon skeleton with two phenyl rings and a heterocyclic ring
    flavonoid_pattern = Chem.MolFromSmarts("c1cc(-c2ccc3occ(=O)c(-c4ccccc4)c3c2)ccc1")
    if mol.HasSubstructMatch(flavonoid_pattern):
        return True, "Contains flavonoid backbone"

    # Lignans: dimers formed by oxidative coupling of phenylpropanoids
    lignan_pattern = Chem.MolFromSmarts("c1cc(O)ccc1C[C@@H](C)c1ccc(O)cc1")
    if mol.HasSubstructMatch(lignan_pattern):
        return True, "Contains lignan structure"

    # Stilbenoids: compounds with a 1,2-diphenylethylene structure
    stilbenoid_pattern = Chem.MolFromSmarts("c1ccc(cc1)/C=C/c2ccccc2")
    if mol.HasSubstructMatch(stilbenoid_pattern):
        return True, "Contains stilbenoid structure"

    # Anthocyanins: flavylium ion backbone
    anthocyanin_pattern = Chem.MolFromSmarts("[O+]c1cc2ccc(O)cc2oc1")
    if mol.HasSubstructMatch(anthocyanin_pattern):
        return True, "Contains anthocyanin structure"

    # General phenylpropanoid pattern: benzene ring with a three-carbon chain (allowing for unsaturation)
    general_pattern = Chem.MolFromSmarts("c1ccccc1[C;H0,H1,H2]=[C;H0,H1,H2][C;H0,H1,H2]")
    if mol.HasSubstructMatch(general_pattern):
        return True, "Contains generalized phenylpropanoid structure"

    # Look for presence of benzene ring and count carbons in side chain
    aromatic_rings = mol.GetRingInfo().AtomRings()
    benzene = Chem.MolFromSmarts("c1ccccc1")
    if mol.HasSubstructMatch(benzene):
        # Find side chains connected to benzene ring
        for bond in mol.GetBonds():
            atom1 = bond.GetBeginAtom()
            atom2 = bond.GetEndAtom()
            if atom1.IsInRing() != atom2.IsInRing():
                side_chain = None
                if atom1.IsInRing():
                    side_chain = atom2
                else:
                    side_chain = atom1
                # Count carbons in side chain
                carbon_count = 0
                stack = [side_chain]
                visited = set()
                while stack:
                    atom = stack.pop()
                    if atom.GetIdx() in visited or atom.IsInRing():
                        continue
                    visited.add(atom.GetIdx())
                    if atom.GetAtomicNum() == 6:
                        carbon_count += 1
                    for neighbor in atom.GetNeighbors():
                        if neighbor.GetIdx() not in visited:
                            stack.append(neighbor)
                if carbon_count == 3:
                    return True, "Contains benzene ring with three-carbon side chain"

    return False, "Does not contain phenylpropanoid skeleton"